"""Rebuild the small offline test fixtures from a populated download cache.

Usage:  PROTEINLENS_CACHE=/path/to/cache python tests/fixtures/make_fixtures.py
(first run e.g. ``proteinlens info 1UBQ`` etc. with that cache so the files exist).
Coordinate files are stripped to ATOM/HETATM/MODEL/ENDMDL/TER/END records and gzipped;
metadata JSON is trimmed to the fields ProteinLens reads.
"""

import gzip
import json
import os
import sys
from pathlib import Path

SRC = Path(os.environ["PROTEINLENS_CACHE"])
DST = Path(__file__).parent / "cache"
KEEP = ("ATOM", "HETATM", "MODEL", "ENDMDL", "TER", "END")
ENTRY_KEYS = ("struct", "exptl", "rcsb_entry_info", "rcsb_entry_container_identifiers")
ENTITY_KEYS = (
    "rcsb_entity_source_organism",
    "rcsb_polymer_entity",
    "rcsb_polymer_entity_container_identifiers",
)


def strip_coords(src: Path, dst: Path) -> None:
    lines = [ln for ln in src.read_text().splitlines() if ln.startswith(KEEP)]
    dst.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(dst, "wt") as fh:
        fh.write("\n".join(lines) + "\n")


def trim_json(src: Path, dst: Path, keys: tuple[str, ...] | None) -> None:
    data = json.loads(src.read_text())
    if keys is not None:
        data = {k: data[k] for k in keys if k in data}
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(json.dumps(data, indent=1))


def main() -> None:
    for sub in ("pdb", "cath", "alphafold"):
        for f in (SRC / sub).glob("*.pdb"):
            strip_coords(f, DST / sub / (f.name + ".gz"))
    for f in (SRC / "meta").glob("*.json"):
        trim_json(f, DST / "meta" / f.name, ENTITY_KEYS if "_entity" in f.name else ENTRY_KEYS)
    for f in (SRC / "alphafold").glob("*.json"):
        trim_json(f, DST / "alphafold" / f.name, None)
    print("fixtures written to", DST, file=sys.stderr)


if __name__ == "__main__":
    main()
