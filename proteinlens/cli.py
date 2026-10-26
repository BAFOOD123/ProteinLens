"""Command-line interface: ``proteinlens compare | info | classify | report | benchmark``."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import typer

from proteinlens.align import SAME_FOLD_TM
from proteinlens.align import compare as compare_structures
from proteinlens.data import ProteinLensError, Structure, load_structure

app = typer.Typer(
    add_completion=False, help="ProteinLens: compare and classify protein structures."
)
REPO = Path(__file__).resolve().parent.parent


def _load(identifier: str, chain: str | None, rng: str | None) -> Structure:
    try:
        return load_structure(identifier, chain=chain, residue_range=rng)
    except ProteinLensError as err:
        typer.secho(f"Error: {err}", fg=typer.colors.RED, err=True)
        raise typer.Exit(code=1) from err


@app.command()
def info(
    identifier: str = typer.Argument(
        ..., help="PDB ID (4HHB, 4HHB:B), UniProt ID, CATH domain or file"
    ),
    chain: str = typer.Option(None, help="Chain ID"),
    range_: str = typer.Option(None, "--range", help="Residue range, e.g. 10-120"),
) -> None:
    """Download (or load from cache) one structure and print what was found."""
    s = _load(identifier, chain, range_)
    typer.echo(f"{s.name}  ({s.source}, chain {s.chain}, {len(s)} residues)")
    for key in (
        "title",
        "description",
        "organism",
        "method",
        "resolution",
        "mean_plddt",
        "ligands",
    ):
        if s.metadata.get(key) not in (None, [], ""):
            typer.echo(f"  {key:12s} {s.metadata[key]}")
    typer.echo(f"  sequence     {s.sequence}")
    for w in s.warnings:
        typer.secho(f"  warning: {w}", fg=typer.colors.YELLOW)


@app.command()
def compare(
    a: str = typer.Argument(..., help="First structure (PDB ID, UniProt ID, CATH domain or file)"),
    b: str = typer.Argument(..., help="Second structure"),
    chain_a: str = typer.Option(None, help="Chain of A"),
    chain_b: str = typer.Option(None, help="Chain of B"),
    range_a: str = typer.Option(None, help="Residue range of A, e.g. 1-100"),
    range_b: str = typer.Option(None, help="Residue range of B"),
    as_json: bool = typer.Option(False, "--json", help="Print the full result as JSON"),
) -> None:
    """Structurally align two proteins (TM-align) and print the scores."""
    sa, sb = _load(a, chain_a, range_a), _load(b, chain_b, range_b)
    r = compare_structures(sa, sb)
    if as_json:
        typer.echo(json.dumps(r.to_dict(), indent=2))
        return
    typer.echo(f"A: {sa.name} ({len(sa)} res)   B: {sb.name} ({len(sb)} res)")
    typer.echo(f"TM-score (normalized by A): {r.tm_score_a:.3f}")
    typer.echo(f"TM-score (normalized by B): {r.tm_score_b:.3f}")
    typer.echo(f"RMSD: {r.rmsd:.2f} Å over {r.aligned_length} aligned residue pairs")
    typer.echo(f"Sequence identity (aligned pairs): {100 * r.seq_identity:.1f}%")
    verdict = "likely the SAME fold" if r.likely_same_fold else "likely DIFFERENT folds"
    typer.echo(f"Verdict: {verdict} (rule: max TM-score > {SAME_FOLD_TM})")
    if r.top_regions:
        typer.echo("Largest deviations after superposition:")
        for g in r.top_regions:
            typer.echo(
                f"  A {g.residues_a} <-> B {g.residues_b}: mean {g.mean_distance:.2f} Å, "
                f"max {g.max_distance:.2f} Å ({g.length} pairs)"
            )
    else:
        typer.echo("No stretch of 3+ aligned residues deviates by more than 3 Å.")
    typer.echo("")
    typer.echo(r.aligned_seq_a)
    typer.echo(r.match_line)
    typer.echo(r.aligned_seq_b)


@app.command()
def benchmark(smoke: bool = typer.Option(False, help="Small, fast subset")) -> None:
    """Build the CATH benchmark set and run the classification benchmark."""
    extra = ["--smoke"] if smoke else []
    for script in ("build_dataset.py", "run_benchmark.py", "analyze.py"):
        args = [sys.executable, str(REPO / "benchmarks" / script)]
        args += extra if script != "build_dataset.py" else []
        subprocess.run(args, check=True)


def main() -> None:  # pragma: no cover
    app()


if __name__ == "__main__":  # pragma: no cover
    main()
