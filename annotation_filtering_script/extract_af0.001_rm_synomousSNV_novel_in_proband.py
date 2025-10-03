#!/usr/bin/env python3
import argparse
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd


def normalize_gt_anyalt(gt_token: str) -> Optional[str]:
    if gt_token in (".", "./.", ".|."):
        return None
    token = gt_token.replace("|", "/").strip()
    if token == "":
        return None
    alleles = token.split("/")
    if len(alleles) == 1:
        alleles = [alleles[0], alleles[0]]
    mapped: List[str] = []
    for a in alleles:
        a = a.strip()
        if a == "" or a == ".":
            return None
        if a.isdigit():
            mapped.append("0" if a == "0" else "1")
        else:
            return None
    has_ref = "0" in mapped
    has_alt = "1" in mapped
    if has_alt and not has_ref:
        return "1/1"
    if has_alt and has_ref:
        return "0/1"
    if has_ref and not has_alt:
        return "0/0"
    return None


def gt_class(gt_norm: str) -> int:
    return {"0/0": 0, "0/1": 1, "1/1": 2}[gt_norm]


def mendelian_consistent_anyalt(child_cls: int, mom_cls: int, dad_cls: int) -> bool:
    if mom_cls == 0 and dad_cls == 0:
        return child_cls == 0
    if mom_cls == 2 and dad_cls == 2:
        return child_cls == 2
    if {mom_cls, dad_cls} == {0, 2}:
        return child_cls == 1
    if mom_cls == 1 and dad_cls == 1:
        return child_cls in (0, 1, 2)
    if {mom_cls, dad_cls} == {0, 1}:
        return child_cls in (0, 1)
    if {mom_cls, dad_cls} == {1, 2}:
        return child_cls in (1, 2)
    return False


def get_gt_index(fmt: str, cache: Dict[str, int]) -> int:
    if not isinstance(fmt, str) or fmt == "":
        return -1
    if fmt in cache:
        return cache[fmt]
    parts = fmt.split(":")
    try:
        idx = parts.index("GT")
    except ValueError:
        idx = -1
    cache[fmt] = idx
    return idx


def extract_gt_from_sample(sample_field: str, gt_idx: int) -> str:
    if not isinstance(sample_field, str) or sample_field == "":
        return ""
    parts = sample_field.split(":")
    if gt_idx < 0 or gt_idx >= len(parts):
        return parts[0] if parts else ""
    return parts[gt_idx]


def autodetect_columns(df_columns: List[str]) -> Tuple[str, str, str, str]:
    cols = list(df_columns)
    name_to_idx = {c: i for i, c in enumerate(cols)}
    if all(k in name_to_idx for k in ("Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15")):
        return "Otherinfo12", "Otherinfo13", "Otherinfo14", "Otherinfo15"
    if all(k in name_to_idx for k in ("Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12")):
        return "Otherinfo9", "Otherinfo10", "Otherinfo11", "Otherinfo12"
    fmt_idx = None
    for cand in cols:
        if "Otherinfo" in cand or cand.upper() == "FORMAT":
            fmt_idx = name_to_idx[cand]
            break
    if fmt_idx is None:
        fmt_idx = max(0, len(cols) - 6)
    fmt_col = cols[fmt_idx]
    if fmt_idx + 3 < len(cols):
        p_col = cols[fmt_idx + 1]
        m_col = cols[fmt_idx + 2]
        f_col = cols[fmt_idx + 3]
    else:
        p_col, m_col, f_col = cols[-3], cols[-2], cols[-1]
    return fmt_col, p_col, m_col, f_col


def exonic_filter_ok(exonic_func: str) -> bool:
    if not isinstance(exonic_func, str):
        return False
    s = exonic_func.strip().lower()
    return s not in ("", ".", "synonymous snv", "unknown")


def process_trio_data(
    input_file: str,
    output_prefix: str,
    format_col: Optional[str] = None,
    proband_col: Optional[str] = None,
    mother_col: Optional[str] = None,
    father_col: Optional[str] = None,
    af_col: str = "gnomad41_genome_AF",
    af_threshold: float = 0.001,
    chunksize: int = 100000,
):
    denovo_file = f"{output_prefix}_denovo.tsv"
    merr_file = f"{output_prefix}_mendelian_errors.tsv"
    altcarriers_file = f"{output_prefix}_altcarriers.tsv"

    counters = dict(
        total=0,
        af_pass=0,
        exonic_pass=0,
        rows_with_gt=0,
        all_norm_ok=0,
        child_has_alt=0,
        denovo=0,
        merr=0,
        invalid_gt=0,
    )

    gt_index_cache: Dict[str, int] = {}
    denovo_header_written = False
    merr_header_written = False
    altcarriers_header_written = False

    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)

    try:
        first_chunk = next(reader)
    except StopIteration:
        print("[INFO] Empty input.", file=sys.stderr)
        return

    if all(col is not None for col in (format_col, proband_col, mother_col, father_col)):
        fmt_c, p_c, m_c, f_c = format_col, proband_col, mother_col, father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(list(first_chunk.columns))

    def handle_chunk(chunk: pd.DataFrame):
        nonlocal denovo_header_written, merr_header_written, altcarriers_header_written
        counters["total"] += len(chunk)
        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        af_mask = af_vals.isna() | (af_vals < af_threshold)
        chunk = chunk[af_mask]
        counters["af_pass"] += int(af_mask.sum())
        if chunk.empty:
            return
        if "ExonicFunc.refGene" in chunk.columns:
            exonic_mask = chunk["ExonicFunc.refGene"].apply(exonic_filter_ok)
            chunk = chunk[exonic_mask]
            counters["exonic_pass"] += int(exonic_mask.sum())
            if chunk.empty:
                return
        for col in (fmt_c, p_c, m_c, f_c):
            if col not in chunk.columns:
                return
        fmt_series = chunk[fmt_c].astype(str)
        p_series = chunk[p_c].astype(str)
        m_series = chunk[m_c].astype(str)
        f_series = chunk[f_c].astype(str)
        counters["rows_with_gt"] += len(chunk)

        results = []
        for fmt, ps, ms, fs in zip(fmt_series, p_series, m_series, f_series):
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                counters["invalid_gt"] += 1
                results.append((None, None, None))
                continue
            p_raw = extract_gt_from_sample(ps, gt_idx)
            m_raw = extract_gt_from_sample(ms, gt_idx)
            f_raw = extract_gt_from_sample(fs, gt_idx)
            p_norm = normalize_gt_anyalt(p_raw)
            m_norm = normalize_gt_anyalt(m_raw)
            f_norm = normalize_gt_anyalt(f_raw)
            if p_norm is None or m_norm is None or f_norm is None:
                counters["invalid_gt"] += 1
                results.append((None, None, None))
                continue
            results.append((p_norm, m_norm, f_norm))

        chunk = chunk.assign(_norm_trio=results)
        ok_mask = chunk["_norm_trio"].apply(lambda t: all(val in ("0/0", "0/1", "1/1") for val in t))
        chunk = chunk[ok_mask]
        counters["all_norm_ok"] += int(ok_mask.sum())
        if chunk.empty:
            return
        has_alt_mask = chunk["_norm_trio"].apply(lambda t: t[0] != "0/0")
        chunk = chunk[has_alt_mask]
        counters["child_has_alt"] += int(has_alt_mask.sum())
        if chunk.empty:
            return

        norm_cols = chunk["_norm_trio"].apply(pd.Series)
        norm_cols.columns = ["NormGT_Proband", "NormGT_Mother", "NormGT_Father"]
        altcarriers_out = pd.concat([chunk.drop(columns=["_norm_trio"]), norm_cols], axis=1)
        altcarriers_out.to_csv(
            altcarriers_file, sep="\t", mode="a", header=(not altcarriers_header_written), index=False
        )
        altcarriers_header_written = True

        def is_denovo(trio) -> bool:
            p, m, d = trio
            return m == "0/0" and d == "0/0" and p in ("0/1", "1/1")

        def is_merr(trio) -> bool:
            p, m, d = trio
            pc, mc, dc = gt_class(p), gt_class(m), gt_class(d)
            return not mendelian_consistent_anyalt(pc, mc, dc)

        denovo_mask = chunk["_norm_trio"].apply(is_denovo)
        merr_mask = chunk["_norm_trio"].apply(lambda t: is_merr(t) and not is_denovo(t))

        denovo_n = int(denovo_mask.sum())
        merr_n = int(merr_mask.sum())
        counters["denovo"] += denovo_n
        counters["merr"] += merr_n

        if denovo_n > 0:
            chunk.loc[denovo_mask].drop(columns=["_norm_trio"]).to_csv(
                denovo_file, sep="\t", mode="a", header=not denovo_header_written, index=False
            )
            denovo_header_written = True

        if merr_n > 0:
            chunk.loc[merr_mask].drop(columns=["_norm_trio"]).to_csv(
                merr_file, sep="\t", mode="a", header=not merr_header_written, index=False
            )
            merr_header_written = True

    handle_chunk(first_chunk)
    for chunk in reader:
        handle_chunk(chunk)

    print("Processing complete. Summary:")
    print(f"  Total variants read:                               {counters['total']}")
    print(f"  Passed AF (<{af_threshold}) filter:                {counters['af_pass']}")
    print(f"  Passed ExonicFunc filter:                          {counters['exonic_pass']}")
    print(f"  Rows with candidate trio columns:                  {counters['rows_with_gt']}")
    print(f"  Fully genotyped and normalized trio rows:          {counters['all_norm_ok']}")
    print(f"  Proband carries ALT (kept):                        {counters['child_has_alt']} -> {altcarriers_file}")
    print(f"  Strict de novo variants:                           {counters['denovo']} -> {denovo_file}")
    print(f"  Mendelian error variants:                          {counters['merr']} -> {merr_file}")
    print(f"  Rows skipped for invalid/missing GT:               {counters['invalid_gt']}")


def main():
    ap = argparse.ArgumentParser(
        description="Filter trio data (Any-ALT collapse) for rare, protein-impacting de novos and Mendelian errors."
    )
    ap.add_argument("-i", "--input", required=True)
    ap.add_argument("-o", "--output-prefix", required=True)
    ap.add_argument("--format-col", default=None)
    ap.add_argument("--proband-col", default=None)
    ap.add_argument("--mother-col", default=None)
    ap.add_argument("--father-col", default=None)
    ap.add_argument("--af-col", default="gnomad41_genome_AF")
    ap.add_argument("--af-threshold", type=float, default=0.001)
    ap.add_argument("--chunksize", type=int, default=100000)
    args = ap.parse_args()

    process_trio_data(
        input_file=args.input,
        output_prefix=args.output_prefix,
        format_col=args.format_col,
        proband_col=args.proband_col,
        mother_col=args.mother_col,
        father_col=args.father_col,
        af_col=args.af_col,
        af_threshold=args.af_threshold,
        chunksize=args.chunksize,
    )


if __name__ == "__main__":
    main()

