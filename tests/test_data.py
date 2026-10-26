"""Tests for the data layer (fetching, caching, parsing, cleaning)."""

import gzip
import io

import numpy as np
import pytest
from Bio.PDB import MMCIFIO, PDBParser

from proteinlens.data import (
    ProteinLensError,
    StructureParseError,
    load_structure,
    parse_residue_range,
    parse_structure,
    split_identifier,
)

UBQ_START = "MQIFVKTLTGK"


def test_three_known_structures_offline(offline_cache):
    ubq = load_structure("1UBQ")
    assert len(ubq) == 76 and ubq.sequence.startswith(UBQ_START)
    assert ubq.metadata["resolution"] == pytest.approx(1.8)
    assert ubq.metadata["organism"] == "Homo sapiens"

    hbb = load_structure("4HHB:B")
    assert hbb.chain == "B" and len(hbb) == 146
    assert "HEM" in hbb.metadata["ligands"]
    assert hbb.metadata["description"] == "Hemoglobin subunit beta"

    tim = load_structure("8TIM")
    assert tim.chain == "A" and len(tim) == 247  # first protein chain by default


def test_alphafold_and_cath_offline(offline_cache):
    af = load_structure("P69905")
    assert af.source == "alphafold" and len(af) == 142
    assert af.metadata["mean_plddt"] > 90  # pLDDT, a 0-100 confidence score
    assert 0 <= af.bfactors.min() and af.bfactors.max() <= 100

    dom = load_structure("1cukA01")
    assert dom.source == "cath" and len(dom) > 50


def test_hetatm_ligands_and_water_are_stripped(offline_cache):
    hba = load_structure("4HHB:A")
    assert "HEM" in hba.metadata["ligands_in_chain"]  # it was there ...
    assert "HEM" not in hba.residue_names  # ... but is not part of the protein
    assert set(hba.atom_residue_index) == set(range(len(hba)))


def test_chain_and_range_selection(offline_cache):
    part = load_structure("1UBQ:A:1-50")
    assert len(part) == 50 and part.residue_ids[0] == "1" and part.residue_ids[-1] == "50"
    same = load_structure("1UBQ", chain="A", residue_range="1-50")
    assert same.sequence == part.sequence


def test_bad_chain_gives_clear_error(offline_cache):
    with pytest.raises(StructureParseError, match="Protein chains available: A, B, C, D"):
        load_structure("4HHB:Z")


def test_unknown_identifier(offline_cache):
    with pytest.raises(ProteinLensError, match="not a file, a PDB ID"):
        load_structure("hello-world")


def test_identifier_and_range_parsing():
    assert split_identifier("4HHB:B:1-100") == ("4HHB", "B", (1, 100))
    assert split_identifier("1UBQ") == ("1UBQ", None, None)
    assert parse_residue_range("-5-40") == (-5, 40)
    with pytest.raises(ProteinLensError):
        parse_residue_range("40-5")


def _ubq_text(offline_cache) -> str:
    with gzip.open(offline_cache / "pdb" / "1UBQ.pdb.gz", "rt") as fh:
        return fh.read()


def test_upload_bytes_and_missing_residues(offline_cache):
    text = _ubq_text(offline_cache)
    # Delete residues 30-34 to simulate unresolved (missing) residues.
    kept = [
        ln for ln in text.splitlines() if not (ln.startswith("ATOM") and 30 <= int(ln[22:26]) <= 34)
    ]
    s = parse_structure("\n".join(kept).encode(), name="upload")
    assert len(s) == 71
    assert s.metadata["n_chain_breaks"] == 1
    assert any("chain break" in w for w in s.warnings)


def test_altlocs_keep_one_copy():
    pdb = (
        "ATOM      1  N   GLY A   1       0.000   0.000   0.000  1.00 10.00           N\n"
        "ATOM      2  CA  GLY A   1       1.450   0.000   0.000  1.00 10.00           C\n"
        "ATOM      3  N   SER A   2       3.000   1.000   0.000  1.00 10.00           N\n"
        "ATOM      4  CA ASER A   2       3.800   1.000   0.000  0.60 10.00           C\n"
        "ATOM      5  CA BSER A   2       3.900   1.200   0.000  0.40 10.00           C\n"
        "END\n"
    )
    s = parse_structure(pdb, name="alt")
    assert s.sequence == "GS"
    assert sum(1 for n in s.atom_names if n == "CA") == 2  # not 3
    assert np.allclose(s.ca_coords[1], [3.8, 1.0, 0.0])  # highest occupancy (A, 0.60)
    assert any("alternate" in w for w in s.warnings)


def test_mmcif_and_pdb_roundtrip(offline_cache, tmp_path):
    ubq = load_structure("1UBQ")
    # PDB writer round trip
    again = parse_structure(ubq.to_pdb(), name="rt")
    assert again.sequence == ubq.sequence
    assert np.allclose(again.ca_coords, ubq.ca_coords, atol=1e-3)
    # mmCIF parsing
    st = PDBParser(QUIET=True).get_structure("u", io.StringIO(ubq.to_pdb()))
    io_ = MMCIFIO()
    io_.set_structure(st)
    cif = tmp_path / "ubq.cif"
    io_.save(str(cif))
    from_cif = load_structure(cif)
    assert from_cif.sequence == ubq.sequence


def test_selenomethionine_maps_to_met():
    pdb = (
        "HETATM    1  N   MSE A   1       0.000   0.000   0.000  1.00 10.00           N\n"
        "HETATM    2  CA  MSE A   1       1.450   0.000   0.000  1.00 10.00           C\n"
        "ATOM      3  N   GLY A   2       3.000   1.000   0.000  1.00 10.00           N\n"
        "ATOM      4  CA  GLY A   2       3.800   1.000   0.000  1.00 10.00           C\n"
        "END\n"
    )
    assert parse_structure(pdb, name="mse").sequence == "MG"


@pytest.mark.network
def test_three_known_structures_online(online_cache):
    assert load_structure("1UBQ").sequence.startswith(UBQ_START)
    assert len(load_structure("P69905")) == 142
    assert load_structure("1cukA01").source == "cath"
    assert (online_cache / "pdb" / "1UBQ.pdb").exists()  # cached for next time
