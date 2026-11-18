"""
variants.py - calculation logic for the Genomic Variant Frequency Calculator.

Pure Python (no heavy dependencies) so it is easy to test and read.

Main functions
--------------
genotype_frequencies(n_hom_ref, n_het, n_hom_alt)  -> dict
allele_frequencies(n_hom_ref, n_het, n_hom_alt)    -> dict
hardy_weinberg(n_hom_ref, n_het, n_hom_alt)        -> dict
analyze_counts(n_hom_ref, n_het, n_hom_alt)        -> dict   (everything at once)
parse_vcf(file_like)                               -> list of variant records
analyze_vcf(file_like, maf_threshold)              -> pandas DataFrame
"""

import math
import re

import pandas as pd


# ---------------------------------------------------------------------------
# Core math (works on genotype counts)
# ---------------------------------------------------------------------------

def _validate(n_hom_ref, n_het, n_hom_alt):
    for name, n in (("n_hom_ref", n_hom_ref), ("n_het", n_het), ("n_hom_alt", n_hom_alt)):
        if n < 0:
            raise ValueError(f"{name} cannot be negative")
    if n_hom_ref + n_het + n_hom_alt == 0:
        raise ValueError("Total number of individuals must be greater than 0")


def genotype_frequencies(n_hom_ref, n_het, n_hom_alt):
    """Proportion of individuals with each genotype (AA, Aa, aa)."""
    _validate(n_hom_ref, n_het, n_hom_alt)
    n = n_hom_ref + n_het + n_hom_alt
    return {
        "hom_ref": n_hom_ref / n,
        "het": n_het / n,
        "hom_alt": n_hom_alt / n,
    }


def allele_frequencies(n_hom_ref, n_het, n_hom_alt):
    """
    Frequency of the reference allele (p) and alternate allele (q).
    Each person carries 2 alleles, so there are 2N alleles in total.
    """
    _validate(n_hom_ref, n_het, n_hom_alt)
    total_alleles = 2 * (n_hom_ref + n_het + n_hom_alt)
    ref_count = 2 * n_hom_ref + n_het
    alt_count = 2 * n_hom_alt + n_het
    p = ref_count / total_alleles
    q = alt_count / total_alleles
    return {
        "p_ref": p,
        "q_alt": q,
        "maf": min(p, q),  # Minor Allele Frequency
        "total_alleles": total_alleles,
    }


def hardy_weinberg(n_hom_ref, n_het, n_hom_alt):
    """
    Chi-square test (1 degree of freedom) comparing observed genotype counts
    to the counts expected under Hardy-Weinberg Equilibrium.

    Expected: AA = p^2 * N, Aa = 2pq * N, aa = q^2 * N
    A small p-value (< 0.05) suggests the variant deviates from HWE, which can
    point to genotyping error, selection, or population structure.
    """
    _validate(n_hom_ref, n_het, n_hom_alt)
    n = n_hom_ref + n_het + n_hom_alt
    freqs = allele_frequencies(n_hom_ref, n_het, n_hom_alt)
    p, q = freqs["p_ref"], freqs["q_alt"]

    expected = [p * p * n, 2 * p * q * n, q * q * n]
    observed = [n_hom_ref, n_het, n_hom_alt]

    # If the variant is monomorphic (p or q = 0) there is nothing to test.
    if p == 0 or q == 0:
        return {
            "expected": expected,
            "observed": observed,
            "chi_square": 0.0,
            "p_value": 1.0,
            "in_hwe": True,
        }

    chi2 = sum((o - e) ** 2 / e for o, e in zip(observed, expected))
    # For 1 degree of freedom: p-value = erfc(sqrt(chi2 / 2))
    p_value = math.erfc(math.sqrt(chi2 / 2))
    return {
        "expected": expected,
        "observed": observed,
        "chi_square": chi2,
        "p_value": p_value,
        "in_hwe": p_value >= 0.05,
    }


def analyze_counts(n_hom_ref, n_het, n_hom_alt):
    """Run every calculation and return one combined dictionary."""
    result = {}
    result.update(genotype_frequencies(n_hom_ref, n_het, n_hom_alt))
    result.update(allele_frequencies(n_hom_ref, n_het, n_hom_alt))
    hwe = hardy_weinberg(n_hom_ref, n_het, n_hom_alt)
    result["hwe_chi_square"] = hwe["chi_square"]
    result["hwe_p_value"] = hwe["p_value"]
    result["in_hwe"] = hwe["in_hwe"]
    result["hwe_expected"] = hwe["expected"]
    result["hwe_observed"] = hwe["observed"]
    result["n_samples"] = n_hom_ref + n_het + n_hom_alt
    return result


# ---------------------------------------------------------------------------
# VCF handling
# ---------------------------------------------------------------------------

_GT_SPLIT = re.compile(r"[/|]")


def _classify_genotype(gt):
    """
    Convert a VCF genotype string into a category.
    Returns 'hom_ref', 'het', 'hom_alt', or None (missing / multi-allelic / haploid).
    """
    allele_calls = _GT_SPLIT.split(gt)
    if len(allele_calls) != 2 or "." in allele_calls:
        return None  # missing or haploid call
    a, b = allele_calls
    if a not in ("0", "1") or b not in ("0", "1"):
        return None  # multi-allelic alt (2, 3, ...) is skipped in this simple tool
    alt_copies = int(a) + int(b)
    return {0: "hom_ref", 1: "het", 2: "hom_alt"}[alt_copies]


def parse_vcf(file_like):
    """
    Read a VCF (text mode file-like object, or list of lines) and return a list
    of dicts, one per bi-allelic variant, with genotype counts.
    """
    records = []
    sample_names = []
    for raw in file_like:
        line = raw.decode("utf-8") if isinstance(raw, bytes) else raw
        line = line.rstrip("\r\n")
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            sample_names = line.split("\t")[9:]
            continue

        fields = line.split("\t")
        if len(fields) < 10:
            continue
        chrom, pos, vid, ref, alt = fields[:5]
        if "," in alt:  # multi-allelic site -> skip
            continue

        format_keys = fields[8].split(":")
        if "GT" not in format_keys:
            continue
        gt_index = format_keys.index("GT")

        counts = {"hom_ref": 0, "het": 0, "hom_alt": 0, "missing": 0}
        for sample_field in fields[9:]:
            parts = sample_field.split(":")
            gt = parts[gt_index] if gt_index < len(parts) else "."
            category = _classify_genotype(gt)
            counts[category if category else "missing"] += 1

        records.append({
            "chrom": chrom,
            "pos": int(pos),
            "id": vid,
            "ref": ref,
            "alt": alt,
            "n_hom_ref": counts["hom_ref"],
            "n_het": counts["het"],
            "n_hom_alt": counts["hom_alt"],
            "n_missing": counts["missing"],
        })
    return records


def classify_maf(maf):
    """Common (>=5%), low-frequency (1-5%), or rare (<1%)."""
    if maf >= 0.05:
        return "common"
    if maf >= 0.01:
        return "low-frequency"
    return "rare"


def analyze_vcf(file_like, maf_threshold=0.05):
    """Parse a VCF and return a DataFrame with frequencies, MAF and HWE per variant."""
    rows = []
    for rec in parse_vcf(file_like):
        n_called = rec["n_hom_ref"] + rec["n_het"] + rec["n_hom_alt"]
        if n_called == 0:
            continue  # nothing called for this variant
        res = analyze_counts(rec["n_hom_ref"], rec["n_het"], rec["n_hom_alt"])
        rows.append({
            "CHROM": rec["chrom"],
            "POS": rec["pos"],
            "ID": rec["id"],
            "REF": rec["ref"],
            "ALT": rec["alt"],
            "N_called": n_called,
            "N_missing": rec["n_missing"],
            "HomRef": rec["n_hom_ref"],
            "Het": rec["n_het"],
            "HomAlt": rec["n_hom_alt"],
            "Ref_freq (p)": round(res["p_ref"], 4),
            "Alt_freq (q)": round(res["q_alt"], 4),
            "MAF": round(res["maf"], 4),
            "Category": classify_maf(res["maf"]),
            "HWE_chi2": round(res["hwe_chi_square"], 4),
            "HWE_p": round(res["hwe_p_value"], 4),
            "HWE_pass": res["in_hwe"],
            "Passes_MAF_filter": res["maf"] >= maf_threshold,
        })
    return pd.DataFrame(rows)
