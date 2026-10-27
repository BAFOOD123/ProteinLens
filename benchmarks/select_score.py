"""Choose the TM-score variant for 1-NN using ONLY the development sets (main + hard).

Rule (fixed before the held-out test set was built): among the candidate normalizations in
`summary.json -> sensitivity_tm_normalization`, pick the one with the highest mean fold accuracy
over the two development sets; ties go to the earlier candidate. The held-out test set
(`make benchmark-test`) then checks whether the choice holds up on data it was not chosen on.
"""

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DEV_SETS = ["results", "results_hard"]


def main() -> None:
    sens = [
        json.loads((ROOT / d / "summary.json").read_text())["sensitivity_tm_normalization"]
        for d in DEV_SETS
    ]
    candidates = list(sens[0])
    table = {
        c: {
            **{d: s[c]["fold"] for d, s in zip(DEV_SETS, sens, strict=True)},
            "mean_fold": round(sum(s[c]["fold"] for s in sens) / len(sens), 4),
        }
        for c in candidates
    }
    best = max(candidates, key=lambda c: (table[c]["mean_fold"], -candidates.index(c)))
    out = {
        "rule": "highest mean fold accuracy over the development sets",
        "dev_sets": DEV_SETS,
        "candidates": table,
        "selected": best,
    }
    (ROOT / "score_selection.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
