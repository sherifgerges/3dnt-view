"""Parse variant files into per-residue case/control counts.

Accepts either notation, auto-detected per line:
  * one-letter   e.g.  A985V   or  p.A985V
  * three-letter e.g.  Ala985Val  or  p.Ala985Val
Case-insensitive; an optional 'p.' prefix is allowed on both.
"""
import re
import pandas as pd

# one-letter: a single letter, a number, a single letter (optional p.)
ONE_RE = re.compile(r"^(?:p\.)?([A-Za-z])(\d+)([A-Za-z])$", re.IGNORECASE)
# three-letter: three letters, a number, three letters (optional p.)
THREE_RE = re.compile(r"^(?:p\.)?([A-Za-z]{3})(\d+)([A-Za-z]{3})$", re.IGNORECASE)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}
ONE_LETTERS = set(THREE_TO_ONE.values())


def parse_one(line):
    """Parse a single variant string in either notation.

    Returns (ref, pos, alt, fmt) where fmt is '1-letter' or '3-letter', or None
    if the line isn't a valid missense variant in either format.
    """
    s = line.strip()
    if not s:
        return None
    # three-letter first (its \d is preceded by 3 letters, so it can't be confused
    # with one-letter, but checking it first avoids any ambiguity)
    m = THREE_RE.match(s)
    if m:
        r, a = m.group(1).upper(), m.group(3).upper()
        if r in THREE_TO_ONE and a in THREE_TO_ONE:
            return THREE_TO_ONE[r], int(m.group(2)), THREE_TO_ONE[a], "3-letter"
        return None
    m = ONE_RE.match(s)
    if m:
        r, a = m.group(1).upper(), m.group(3).upper()
        if r in ONE_LETTERS and a in ONE_LETTERS:
            return r, int(m.group(2)), a, "1-letter"
        return None
    return None


def parse_record(line):
    """Parse one line that may be plain (`A985V`) or delimited (.tsv/.csv).

    Splits on tab/comma if present, scans the fields for a variant token, and
    treats a lone integer field (if any) as an allele count for that variant.
    Returns (ref, pos, alt, fmt, count) or None.
    """
    s = line.strip()
    if not s:
        return None
    fields = re.split(r"[\t,]", s) if ("\t" in s or "," in s) else [s]
    fields = [f.strip() for f in fields if f.strip()]
    variant, count = None, None
    for f in fields:
        if variant is None:
            p = parse_one(f)
            if p:
                variant = p
                continue
        if count is None and f.isdigit():
            count = int(f)
    if variant is None:
        return None
    ref, pos, alt, fmt = variant
    return ref, pos, alt, fmt, (count if count is not None else 1)


def parse_variants(text):
    """Return (variants, errors, formats).

    Accepts plain .txt (one variant per line) and tab/comma-delimited .tsv/.csv
    (variant in some column, optional integer count column; header rows are
    skipped automatically since they contain no valid variant token).

    variants: list of {ref, pos, alt, fmt, count} for each well-formed line.
    errors:   list of (line_number, raw_line) with no recognizable variant.
    formats:  dict counting how many lines used each notation.
    Blank lines are ignored.
    """
    variants, errors, formats = [], [], {}
    for i, raw in enumerate(text.splitlines(), start=1):
        if not raw.strip():
            continue
        parsed = parse_record(raw)
        if parsed is None:
            errors.append((i, raw))
            continue
        ref, pos, alt, fmt, count = parsed
        variants.append({"ref": ref, "pos": pos, "alt": alt, "fmt": fmt, "count": count})
        formats[fmt] = formats.get(fmt, 0) + 1
    return variants, errors, formats


def validate_against_sequence(variants, sequence):
    """Return list of mismatch warnings where ref AA != sequence[pos-1]."""
    warnings = []
    n = len(sequence)
    for v in variants:
        p = v["pos"]
        if p < 1 or p > n:
            warnings.append(f"{v['ref']}{p}{v['alt']}: position {p} is outside the protein (length {n})")
        elif sequence[p - 1] != v["ref"]:
            warnings.append(
                f"{v['ref']}{p}{v['alt']}: expected {sequence[p-1]} at position {p} "
                f"(possible isoform/numbering mismatch)")
    return warnings


def build_counts(case_variants, control_variants, uniprot_id):
    """Aggregate per-residue allele counts. Each variant line = one allele.

    Returns a DataFrame with columns uniprot_id, aa_pos, ac_case, ac_control
    (one row per residue that carries at least one case or control variant).
    """
    counts = {}
    for v in case_variants:
        counts.setdefault(v["pos"], [0, 0])[0] += v.get("count", 1)
    for v in control_variants:
        counts.setdefault(v["pos"], [0, 0])[1] += v.get("count", 1)
    rows = [
        {"uniprot_id": uniprot_id, "aa_pos": pos, "ac_case": c[0], "ac_control": c[1]}
        for pos, c in sorted(counts.items())
    ]
    return pd.DataFrame(rows)
