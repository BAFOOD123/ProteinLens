"""Faithfulness checker: verify that a written report only states computed facts.

Given the report text and the facts object it was supposed to be based on, the
checker extracts and verifies:

1. **Every number.** Each number in the text must match some number in the facts,
   allowing for rounding to the precision that was printed: "0.90" matches 0.904
   (tolerance ±0.005), "1.4 Å" matches 1.41 (±0.05), "43%" matches 43.2 (±0.5).
   Fractions in the facts (e.g. a TM-score of 0.904) may also be written as
   percentages (90.4%).
2. **Every identifier** (tokens mixing letters and digits such as PDB IDs ``4HHB``,
   UniProt accessions, CATH domain IDs) and every **CATH code** (``3.10.20``) must
   appear in the facts.
3. **Fold names.** A known fold name (e.g. "TIM barrel", "Rossmann", "globin")
   may only be mentioned if it appears in the facts.
4. **The fold verdict.** If the text says the proteins share a fold (or do not),
   that must agree with the computed verdict (larger TM-score > 0.5).
5. **Outside knowledge about function.** Words about biological function
   ("binds", "catalyzes", "oxygen", ...) are flagged unless they occur in the
   facts, because nothing computed here supports them.

What it does NOT guarantee: that the wording is the clearest possible, that
claims phrased without numbers or listed keywords are correct, or that the
facts themselves are biologically meaningful. See docs/UNDERSTANDING.md.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Any

# --------------------------------------------------------------------------------------
# Vocabulary
# --------------------------------------------------------------------------------------

# Fold / architecture names a writer might mention from memory. Mentioning one is a
# structural claim; it is allowed only if the term appears in the facts.
FOLD_TERMS = [
    "tim barrel",
    "tim-barrel",
    "rossmann",
    "globin",
    "immunoglobulin",
    "ig-like",
    "ubiquitin",
    "beta-grasp",
    "β-grasp",
    "ob fold",
    "ob-fold",
    "jelly roll",
    "jellyroll",
    "beta-propeller",
    "β-propeller",
    "beta propeller",
    "four-helix bundle",
    "4-helix bundle",
    "helix-turn-helix",
    "winged helix",
    "ferredoxin",
    "sh3",
    "sh2",
    "pdz",
    "zinc finger",
    "leucine-rich repeat",
    "ankyrin",
    "armadillo",
    "heat repeat",
    "tpr",
    "lipocalin",
    "cyclophilin",
    "beta-trefoil",
    "β-trefoil",
    "coiled coil",
    "coiled-coil",
    "flavodoxin",
    "p-loop",
    "kinase",
    "thioredoxin",
    "periplasmic binding",
]

# Words about biological function / mechanism that nothing in ProteinLens computes.
FUNCTION_TERMS = [
    "bind",
    "binds",
    "binding",
    "catalyz",
    "catalys",
    "enzym",
    "oxygen",
    "transport",
    "regulat",
    "signal",
    "inhibit",
    "activat",
    "metaboli",
    "evolutionar",
    "ancestor",
    "ancestral",
    "diverged",
    "homolog",
    "mutation",
    "disease",
    "function",
]

SAME_FOLD_PATTERNS = [
    r"same fold",
    r"shared fold",
    r"common fold",
    r"similar fold",
    r"same overall fold",
    r"share (?:a|the|their|one)?\s*(?:common |same |similar |overall )*fold",
    r"structurally similar",
    r"same structural class",
]
DIFFERENT_FOLD_PATTERNS = [
    r"different folds?",
    r"distinct folds?",
    r"unrelated folds?",
    r"dissimilar folds?",
    r"structurally unrelated",
    r"structurally dissimilar",
    r"structurally different",
]
NEGATION = re.compile(r"\b(not|no|never|neither|nor|unlikely|cannot|n't)\b", re.I)

CATH_CODE_RE = re.compile(r"\b\d+\.\d+\.\d+(?:\.\d+)?\b")
ID_RE = re.compile(r"\b(?=[A-Za-z0-9]*\d)(?=[A-Za-z0-9]*[A-Za-z])[A-Za-z0-9]{4,}\b")
NUMBER_RE = re.compile(r"(?<![A-Za-z0-9.])(\d+(?:\.\d+)?)(\s*%)?")
HARMLESS_TOKENS = {"3d", "2d", "1nn"}  # tokens that look like IDs but are vocabulary


# --------------------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------------------


@dataclass
class Issue:
    kind: str  # "number", "identifier", "fold_name", "verdict", "outside_knowledge"
    text: str  # the offending piece of text
    detail: str


@dataclass
class CheckResult:
    passed: bool
    issues: list[Issue] = field(default_factory=list)
    numbers_checked: int = 0
    identifiers_checked: int = 0
    verdict_claims: int = 0

    def summary(self) -> str:
        if self.passed:
            return (
                f"verified: {self.numbers_checked} numbers, {self.identifiers_checked} "
                f"identifiers, {self.verdict_claims} fold-verdict claims"
            )
        return "; ".join(f"[{i.kind}] {i.text}: {i.detail}" for i in self.issues)


# --------------------------------------------------------------------------------------
# Collect what the facts allow
# --------------------------------------------------------------------------------------


def _walk(obj: Any, key: str = ""):
    """Yield (key, value) for every leaf of a nested dict/list."""
    if isinstance(obj, dict):
        for k, v in obj.items():
            yield from _walk(v, k)
    elif isinstance(obj, list | tuple):
        for v in obj:
            yield from _walk(v, key)
    else:
        yield key, obj


@dataclass
class Allowed:
    numbers: list[float]
    fractions: list[float]  # values in [0, 1] that may be written as percentages
    identifiers: set[str]
    cath_codes: set[str]
    text_blob: str  # all string facts, lower-cased, for term lookups
    likely_same_fold: bool | None


def allowed_from_facts(facts: dict[str, Any]) -> Allowed:
    numbers: list[float] = []
    fractions: list[float] = []
    strings: list[str] = []
    for key, v in _walk(facts):
        if isinstance(v, bool) or v is None:
            continue
        if isinstance(v, int | float):
            numbers.append(float(v))
            if 0 <= v <= 1 and "percent" not in key:
                fractions.append(float(v))
        elif isinstance(v, str):
            strings.append(v)
            # Numbers inside strings (titles, residue ranges like "45-52") are facts too.
            numbers.extend(float(m.group(1)) for m in NUMBER_RE.finditer(v))
    blob = " ".join(strings)
    idents = {t.lower() for t in ID_RE.findall(blob)}
    # Names like "4HHB:A" or "AF-P69905-F1" split into parts that are each allowed.
    for s in strings:
        idents.update(p.lower() for p in re.split(r"[^A-Za-z0-9]+", s) if p)
    codes = set(CATH_CODE_RE.findall(blob))
    verdict = facts.get("interpretation", {}).get("likely_same_fold")
    return Allowed(numbers, fractions, idents, codes, blob.lower(), verdict)


def _matches(value: float, decimals: int, is_percent: bool, allowed: Allowed) -> bool:
    tol = 0.5 * 10 ** (-decimals) + 1e-9
    candidates = list(allowed.numbers)
    if is_percent:
        candidates += [f * 100 for f in allowed.fractions]
    return any(abs(value - c) <= tol for c in candidates)


# --------------------------------------------------------------------------------------
# The checker
# --------------------------------------------------------------------------------------


def _sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?;])\s+|\n+", text) if s.strip()]


def _negated(low: str, start: int) -> bool:
    """Is the phrase starting at ``start`` negated by one of the 4 words before it
    (within the same clause)? e.g. "do not share the same fold"."""
    clause = re.split(r"[,;:()]", low[:start])[-1]
    return bool(NEGATION.search(" ".join(clause.split()[-4:])))


def fold_verdict_claims(text: str) -> list[tuple[str, bool]]:
    """Find sentences that claim same/different fold. Returns (sentence, claims_same)."""
    claims = []
    for sent in _sentences(text):
        low = sent.lower()
        found: list[bool] = []  # True = asserts same fold
        for pats, is_same in ((SAME_FOLD_PATTERNS, True), (DIFFERENT_FOLD_PATTERNS, False)):
            for p in pats:
                for m in re.finditer(p, low):
                    found.append(is_same != _negated(low, m.start()))
        if found and all(f == found[0] for f in found):
            claims.append((sent.strip(), found[0]))
        # Sentences making both kinds of claim (e.g. contrasting cases) are too ambiguous
        # to score automatically and are skipped.
    return claims


def check_report(
    text: str, facts: dict[str, Any], claimed_same_fold: bool | None = None
) -> CheckResult:
    """Verify ``text`` against ``facts``. ``claimed_same_fold`` is the structured verdict
    (if the writer returned one); it must agree with the computed verdict."""
    allowed = allowed_from_facts(facts)
    issues: list[Issue] = []
    masked = text

    # 1. CATH codes (checked first, then masked so their digits are not read as numbers).
    for code in CATH_CODE_RE.findall(text):
        if code not in allowed.cath_codes:
            issues.append(Issue("identifier", code, "CATH code not in the computed facts"))
    masked = CATH_CODE_RE.sub(" ", masked)

    # 2. Identifiers (PDB IDs, accessions, domain IDs ...).
    n_ids = 0
    for tok in ID_RE.findall(masked):
        if tok.lower() in HARMLESS_TOKENS:
            continue
        n_ids += 1
        if tok.lower() not in allowed.identifiers:
            issues.append(Issue("identifier", tok, "identifier not in the computed facts"))
    masked = ID_RE.sub(" ", masked)

    # 3. Numbers.
    n_nums = 0
    for m in NUMBER_RE.finditer(masked):
        raw, pct = m.group(1), bool(m.group(2))
        decimals = len(raw.split(".")[1]) if "." in raw else 0
        n_nums += 1
        if not _matches(float(raw), decimals, pct, allowed):
            shown = raw + ("%" if pct else "")
            issues.append(Issue("number", shown, "does not match any computed value"))

    # 4. Fold names from memory.
    low = text.lower()
    for term in FOLD_TERMS:
        if re.search(rf"(?<![a-z]){re.escape(term)}", low) and term not in allowed.text_blob:
            issues.append(Issue("fold_name", term, "fold name not supported by the facts"))

    # 5. Functional / evolutionary claims.
    for term in FUNCTION_TERMS:
        for m in re.finditer(rf"\b{re.escape(term)}\w*", low):
            if m.group(0) not in allowed.text_blob:
                issues.append(
                    Issue("outside_knowledge", m.group(0), "functional claim not computed")
                )
                break

    # 6. Same-fold verdict.
    claims = fold_verdict_claims(text)
    if allowed.likely_same_fold is not None:
        for sent, says_same in claims:
            if says_same != allowed.likely_same_fold:
                issues.append(Issue("verdict", sent, "contradicts the computed fold verdict"))
        if claimed_same_fold is not None and claimed_same_fold != allowed.likely_same_fold:
            issues.append(
                Issue("verdict", f"same_fold={claimed_same_fold}", "structured verdict is wrong")
            )

    return CheckResult(
        passed=not issues,
        issues=issues,
        numbers_checked=n_nums,
        identifiers_checked=n_ids,
        verdict_claims=len(claims) + (claimed_same_fold is not None),
    )
