"""Data layer: fetch, cache and parse protein structures.

Supported inputs
----------------
* PDB ID            e.g. ``"4HHB"``, ``"4HHB:B"`` (chain B), ``"4HHB:B:1-100"`` (+ residues 1-100)
* UniProt accession e.g. ``"P69905"`` -> predicted structure from AlphaFold DB
* CATH domain ID    e.g. ``"1cukA01"`` -> a single domain cut out by CATH
* A local file      ``.pdb``, ``.cif``/``.mmcif`` (optionally ``.gz``), e.g. an upload

Everything that is downloaded is cached under ``data/cache/`` (override with the
``PROTEINLENS_CACHE`` environment variable), so the second time you ask for a
structure no network access is needed.

What "cleaning" a structure means here
--------------------------------------
* Only the first model is used (NMR files contain many models).
* One chain is selected (the first protein chain if none is given).
* Water and ligands (HETATM records) are removed; modified amino acids that are
  part of the chain, like selenomethionine (MSE), are kept and mapped to their
  parent amino acid.
* For atoms with alternate locations (alt locs) Biopython keeps the one with the
  highest occupancy; we write out only that one.
* Residues without a C-alpha atom are skipped (with a warning), because the
  structural alignment works on C-alpha atoms.
* Missing residues (gaps in the experimental model) are detected as "chain
  breaks": consecutive C-alpha atoms further apart than 4.5 Å.
"""

from __future__ import annotations

import gzip
import io
import json
import os
import re
import time
import urllib.error
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
from Bio.Data.PDBData import protein_letters_3to1_extended
from Bio.PDB import MMCIFParser, PDBParser

# --------------------------------------------------------------------------------------
# Errors
# --------------------------------------------------------------------------------------


class ProteinLensError(Exception):
    """Base class for errors with a clear, user-facing message."""


class StructureNotFoundError(ProteinLensError):
    """The requested structure could not be found or downloaded."""


class StructureParseError(ProteinLensError):
    """The structure file could not be parsed into a usable protein chain."""


# --------------------------------------------------------------------------------------
# Constants and identifier patterns
# --------------------------------------------------------------------------------------

PDB_ID_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
# Official UniProt accession pattern (from uniprot.org help pages).
UNIPROT_RE = re.compile(r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$")
# CATH domain: 4-char PDB ID + chain letter + 2-digit domain number, e.g. 1cukA01.
CATH_DOMAIN_RE = re.compile(r"^[0-9][A-Za-z0-9]{3}[A-Za-z0-9][0-9]{2}$")

CATH_VERSION = "v4_4_0"
CA_CA_BREAK = 4.5  # Å; consecutive C-alpha atoms are ~3.8 Å apart; one missing residue gives >5 Å
USER_AGENT = "ProteinLens/0.1 (https://github.com/bafood123/proteinlens)"


def cache_dir() -> Path:
    """Root directory of the download cache (created on demand)."""
    root = os.environ.get("PROTEINLENS_CACHE")
    path = Path(root) if root else Path(__file__).resolve().parent.parent / "data" / "cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


# --------------------------------------------------------------------------------------
# The Structure object
# --------------------------------------------------------------------------------------


@dataclass
class Structure:
    """One cleaned protein chain (or domain).

    Per-residue arrays all have length ``len(self)`` and are in chain order.
    Per-atom arrays describe every kept heavy atom, so the structure can be
    written back out (e.g. for the 3D viewer).
    """

    name: str
    source: str  # "pdb", "alphafold", "cath" or "file"
    chain: str
    sequence: str  # one-letter codes, only residues that have a C-alpha atom
    ca_coords: np.ndarray  # (N, 3) float64, Å
    residue_ids: list[str]  # author residue numbers incl. insertion code, e.g. "52", "52A"
    bfactors: np.ndarray  # (N,) C-alpha B-factor; for AlphaFold models this is pLDDT
    atom_names: list[str] = field(repr=False)
    atom_elements: list[str] = field(repr=False)
    atom_coords: np.ndarray = field(repr=False)  # (M, 3)
    atom_residue_index: np.ndarray = field(repr=False)  # (M,) index into the residue arrays
    residue_names: list[str] = field(repr=False)  # three-letter names
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.sequence)

    @property
    def chain_breaks(self) -> list[int]:
        """Indices i where the gap between residue i and i+1 is a chain break."""
        if len(self) < 2:
            return []
        d = np.linalg.norm(np.diff(self.ca_coords, axis=0), axis=1)
        return [int(i) for i in np.where(d > CA_CA_BREAK)[0]]

    def to_pdb(
        self,
        rotation: np.ndarray | None = None,
        translation: np.ndarray | None = None,
        chain_id: str | None = None,
    ) -> str:
        """Write the cleaned structure as PDB text, optionally moved by x' = R·x + t."""
        coords = self.atom_coords
        if rotation is not None:
            coords = coords @ np.asarray(rotation).T
        if translation is not None:
            coords = coords + np.asarray(translation)
        chain = (chain_id or self.chain or "A")[:1]
        lines = []
        for i, (name, elem, xyz, r) in enumerate(
            zip(self.atom_names, self.atom_elements, coords, self.atom_residue_index, strict=True),
            start=1,
        ):
            resid = self.residue_ids[r]
            num = int(re.match(r"-?\d+", resid).group())  # type: ignore[union-attr]
            icode = resid[len(str(num)) :][:1] or " "
            atom_field = f" {name:<3}" if len(name) < 4 else name[:4]
            lines.append(
                f"ATOM  {i % 100000:5d} {atom_field} {self.residue_names[r]:>3} {chain}"
                f"{num % 10000:4d}{icode}   {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
                f"{1.0:6.2f}{float(self.bfactors[r]):6.2f}          {elem:>2}"
            )
        lines.append("TER")
        lines.append("END")
        return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------------------
# Downloading (with cache)
# --------------------------------------------------------------------------------------


def _download(url: str, dest: Path, retries: int = 3, timeout: float = 60) -> Path:
    """Download ``url`` to ``dest`` unless it is already cached. Raises on failure."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    last_err: Exception | None = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = resp.read()
            tmp = dest.with_suffix(dest.suffix + ".part")
            tmp.write_bytes(data)
            tmp.replace(dest)  # atomic: never leave a half-written file in the cache
            return dest
        except urllib.error.HTTPError as err:
            if err.code == 404:
                raise StructureNotFoundError(f"Not found (HTTP 404): {url}") from err
            last_err = err
        except (urllib.error.URLError, TimeoutError, OSError) as err:
            last_err = err
        time.sleep(2**attempt)
    raise StructureNotFoundError(f"Could not download {url} after {retries} tries: {last_err}")


def _get_json(url: str, dest: Path) -> Any:
    return json.loads(_download(url, dest).read_text())


def fetch_pdb_file(pdb_id: str) -> Path:
    """Download (or reuse from cache) a PDB entry. Tries PDB format, then mmCIF."""
    pdb_id = pdb_id.upper()
    folder = cache_dir() / "pdb"
    for ext in ("pdb", "cif", "pdb.gz", "cif.gz"):
        cached = folder / f"{pdb_id}.{ext}"
        if cached.exists():
            return cached
    try:
        return _download(f"https://files.rcsb.org/download/{pdb_id}.pdb", folder / f"{pdb_id}.pdb")
    except StructureNotFoundError:
        # Very large entries are only distributed as mmCIF.
        try:
            return _download(
                f"https://files.rcsb.org/download/{pdb_id}.cif", folder / f"{pdb_id}.cif"
            )
        except StructureNotFoundError as err:
            raise StructureNotFoundError(
                f"PDB entry '{pdb_id}' was not found at RCSB (files.rcsb.org). "
                "Check the ID, or upload the file instead."
            ) from err


def fetch_pdb_metadata(pdb_id: str) -> dict[str, Any]:
    """Title, method, resolution, ligands and per-chain names/organisms from the RCSB API."""
    pdb_id = pdb_id.upper()
    folder = cache_dir() / "meta"
    entry = _get_json(
        f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}", folder / f"{pdb_id}.json"
    )
    info = entry.get("rcsb_entry_info", {})
    resolution = (info.get("resolution_combined") or [None])[0]
    meta: dict[str, Any] = {
        "pdb_id": pdb_id,
        "title": entry.get("struct", {}).get("title"),
        "method": ", ".join(e.get("method", "") for e in entry.get("exptl", [])) or None,
        "resolution": resolution,
        "ligands": info.get("nonpolymer_bound_components") or [],
        "chains": {},
    }
    ids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids") or []
    for ent_id in ids:
        ent = _get_json(
            f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{ent_id}",
            folder / f"{pdb_id}_entity{ent_id}.json",
        )
        org = (ent.get("rcsb_entity_source_organism") or [{}])[0].get("ncbi_scientific_name")
        desc = ent.get("rcsb_polymer_entity", {}).get("pdbx_description")
        for ch in ent.get("rcsb_polymer_entity_container_identifiers", {}).get("auth_asym_ids", []):
            meta["chains"][ch] = {"description": desc, "organism": org}
    return meta


def fetch_alphafold(uniprot: str) -> tuple[Path, dict[str, Any]]:
    """Download the AlphaFold DB model for a UniProt accession (and its summary info)."""
    uniprot = uniprot.upper()
    folder = cache_dir() / "alphafold"
    try:
        entries = _get_json(
            f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}", folder / f"{uniprot}.json"
        )
    except StructureNotFoundError as err:
        raise StructureNotFoundError(
            f"No AlphaFold DB model found for UniProt accession '{uniprot}'."
        ) from err
    if not entries:
        raise StructureNotFoundError(f"AlphaFold DB returned no models for '{uniprot}'.")
    e = entries[0]
    gz = folder / f"{e['entryId']}.pdb.gz"
    path = gz if gz.exists() else _download(e["pdbUrl"], folder / f"{e['entryId']}.pdb")
    meta = {
        "uniprot": uniprot,
        "title": f"AlphaFold model {e['entryId']}: {e.get('uniprotDescription')}",
        "method": "AlphaFold prediction",
        "model_id": e["entryId"],
        "mean_plddt": e.get("globalMetricValue"),
        "organism": e.get("organismScientificName"),
        "gene": e.get("gene"),
        "description": e.get("uniprotDescription"),
        "resolution": None,
        "ligands": [],
    }
    return path, meta


def fetch_cath_domain(domain_id: str, version: str = CATH_VERSION) -> Path:
    """Download one CATH domain (already cut out of its chain by CATH) as a PDB file."""
    url = f"https://www.cathdb.info/version/{version}/api/rest/id/{domain_id}.pdb"
    gz = cache_dir() / "cath" / f"{domain_id}.pdb.gz"
    if gz.exists():
        return gz
    try:
        return _download(url, cache_dir() / "cath" / f"{domain_id}.pdb")
    except StructureNotFoundError as err:
        raise StructureNotFoundError(f"CATH domain '{domain_id}' could not be fetched.") from err


# --------------------------------------------------------------------------------------
# Parsing
# --------------------------------------------------------------------------------------


def three_to_one(resname: str) -> str | None:
    """Map a residue name to a one-letter code (MSE -> M etc.); None if not an amino acid."""
    return protein_letters_3to1_extended.get(resname.upper())


def parse_residue_range(text: str | None) -> tuple[int, int] | None:
    """Parse ``"10-120"`` into (10, 120). Negative numbers are allowed (``"-5-40"``)."""
    if text is None or text == "":
        return None
    m = re.fullmatch(r"\s*(-?\d+)\s*-\s*(-?\d+)\s*", text)
    if not m:
        raise ProteinLensError(f"Residue range '{text}' should look like '10-120'.")
    start, end = int(m.group(1)), int(m.group(2))
    if end < start:
        raise ProteinLensError(f"Residue range '{text}': end is before start.")
    return start, end


def _read_text(path: Path) -> str:
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as fh:
            return fh.read()
    return path.read_text()


def _looks_like_cif(text: str, filename: str = "") -> bool:
    name = filename.lower().removesuffix(".gz")
    if name.endswith((".cif", ".mmcif")):
        return True
    if name.endswith((".pdb", ".ent")):
        return False
    return text.lstrip().startswith("data_")


def parse_structure(
    source: str | Path | bytes,
    chain: str | None = None,
    residue_range: tuple[int, int] | None = None,
    name: str | None = None,
    origin: str = "file",
    filename: str = "",
) -> Structure:
    """Parse PDB/mmCIF content into a cleaned single-chain :class:`Structure`.

    ``source`` may be a path, raw bytes (e.g. from an upload) or text.
    """
    if isinstance(source, Path) or (isinstance(source, str) and "\n" not in source):
        path = Path(source)
        if not path.exists():
            raise StructureNotFoundError(f"File not found: {path}")
        text = _read_text(path)
        filename = filename or path.name
        name = name or path.name.split(".")[0]
    elif isinstance(source, bytes):
        text = gzip.decompress(source).decode() if source[:2] == b"\x1f\x8b" else source.decode()
    else:
        text = source
    name = name or "structure"

    parser = MMCIFParser(QUIET=True) if _looks_like_cif(text, filename) else PDBParser(QUIET=True)
    try:
        model_set = parser.get_structure(name, io.StringIO(text))
        model = next(iter(model_set))
    except StopIteration as err:
        raise StructureParseError(f"'{name}' contains no models/atoms.") from err
    except Exception as err:  # Biopython raises many exception types on malformed input
        raise StructureParseError(f"Could not parse '{name}': {err}") from err

    protein_chains = [
        c.id for c in model if any(three_to_one(r.get_resname()) and "CA" in r for r in c)
    ]
    if not protein_chains:
        raise StructureParseError(f"'{name}' contains no protein chain with C-alpha atoms.")
    if chain is None:
        chain = protein_chains[0]
    elif chain not in protein_chains:
        raise StructureParseError(
            f"Chain '{chain}' not found in '{name}'. Protein chains available: "
            f"{', '.join(protein_chains)}"
        )

    warnings: list[str] = []
    seq, ca, resids, bfac, resnames = [], [], [], [], []
    a_names, a_elems, a_coords, a_res = [], [], [], []
    ligands: set[str] = set()
    n_no_ca = 0
    n_altloc = 0
    for res in model[chain]:
        het, num, icode = res.id
        resname = res.get_resname()
        if het == "W":
            continue
        one = three_to_one(resname)
        if one is None:  # ligand, ion, sugar, nucleotide ...
            ligands.add(resname)
            continue
        if residue_range and not (residue_range[0] <= num <= residue_range[1]):
            continue
        if "CA" not in res:
            n_no_ca += 1
            continue
        idx = len(seq)
        seq.append(one if len(one) == 1 else "X")
        ca.append(res["CA"].get_coord())
        resids.append(f"{num}{icode.strip()}")
        bfac.append(res["CA"].get_bfactor())
        resnames.append(resname)
        for atom in res:  # for alt locs Biopython yields the highest-occupancy copy
            if atom.is_disordered():
                n_altloc += 1
            if atom.element == "H":
                continue
            a_names.append(atom.get_name())
            a_elems.append(atom.element or atom.get_name()[0])
            a_coords.append(atom.get_coord())
            a_res.append(idx)

    if not seq:
        where = f" in residues {residue_range[0]}-{residue_range[1]}" if residue_range else ""
        raise StructureParseError(f"Chain '{chain}' of '{name}' has no amino acids{where}.")
    if n_no_ca:
        warnings.append(f"{n_no_ca} residue(s) without a C-alpha atom were skipped.")
    if n_altloc:
        warnings.append(f"{n_altloc} atom(s) had alternate locations; kept the highest occupancy.")

    s = Structure(
        name=name,
        source=origin,
        chain=chain,
        sequence="".join(seq),
        ca_coords=np.asarray(ca, dtype=float),
        residue_ids=resids,
        bfactors=np.asarray(bfac, dtype=float),
        atom_names=a_names,
        atom_elements=a_elems,
        atom_coords=np.asarray(a_coords, dtype=float),
        atom_residue_index=np.asarray(a_res, dtype=int),
        residue_names=resnames,
        metadata={"ligands_in_chain": sorted(ligands), "protein_chains": protein_chains},
        warnings=warnings,
    )
    if s.chain_breaks:
        s.warnings.append(
            f"{len(s.chain_breaks)} chain break(s) (missing residues) detected; "
            "alignment uses only the residues present."
        )
    s.metadata["n_chain_breaks"] = len(s.chain_breaks)
    if residue_range:
        s.metadata["residue_range"] = f"{residue_range[0]}-{residue_range[1]}"
    return s


# --------------------------------------------------------------------------------------
# One entry point for everything
# --------------------------------------------------------------------------------------


def split_identifier(identifier: str) -> tuple[str, str | None, tuple[int, int] | None]:
    """``"4HHB:B:1-100"`` -> ("4HHB", "B", (1, 100)). Chain and range are optional."""
    parts = identifier.strip().split(":")
    base = parts[0]
    chain = parts[1] if len(parts) > 1 and parts[1] else None
    rng = parse_residue_range(parts[2]) if len(parts) > 2 else None
    return base, chain, rng


def load_structure(
    identifier: str | Path,
    chain: str | None = None,
    residue_range: str | tuple[int, int] | None = None,
    with_metadata: bool = True,
) -> Structure:
    """Load a structure from a PDB ID, UniProt accession, CATH domain ID or file path.

    Explicit ``chain``/``residue_range`` arguments override anything in the identifier.
    """
    rng = parse_residue_range(residue_range) if isinstance(residue_range, str) else residue_range
    path = Path(identifier)
    if isinstance(identifier, Path) or path.suffix.lower() in {
        ".pdb",
        ".cif",
        ".mmcif",
        ".ent",
        ".gz",
    }:
        s = parse_structure(path, chain=chain, residue_range=rng, origin="file")
        s.metadata.setdefault("title", path.name)
        return s

    base, id_chain, id_rng = split_identifier(str(identifier))
    chain = chain or id_chain
    rng = rng or id_rng

    if CATH_DOMAIN_RE.match(base) and not UNIPROT_RE.match(base.upper()):
        s = parse_structure(
            fetch_cath_domain(base), chain=chain, residue_range=rng, name=base, origin="cath"
        )
        s.metadata["title"] = f"CATH domain {base}"
        return s

    if PDB_ID_RE.match(base):
        pdb_id = base.upper()
        s = parse_structure(
            fetch_pdb_file(pdb_id), chain=chain, residue_range=rng, name=pdb_id, origin="pdb"
        )
        s.name = f"{pdb_id}:{s.chain}"
        if with_metadata:
            try:
                meta = fetch_pdb_metadata(pdb_id)
                chain_meta = meta.pop("chains").get(s.chain, {})
                s.metadata.update(meta)
                s.metadata.update(chain_meta)
            except ProteinLensError as err:
                s.warnings.append(f"Metadata unavailable: {err}")
        return s

    if UNIPROT_RE.match(base.upper()):
        af_path, meta = fetch_alphafold(base)
        s = parse_structure(
            af_path, chain=chain, residue_range=rng, name=meta["model_id"], origin="alphafold"
        )
        s.metadata.update(meta)
        # In AlphaFold files the B-factor column holds pLDDT (0-100 confidence).
        s.metadata["mean_plddt_chain"] = round(float(s.bfactors.mean()), 2)
        return s

    raise ProteinLensError(
        f"'{identifier}' is not a file, a PDB ID (e.g. 4HHB or 4HHB:A), a UniProt "
        "accession (e.g. P69905) or a CATH domain ID (e.g. 1cukA01)."
    )
