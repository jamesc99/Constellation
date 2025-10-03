#!/usr/bin/env python3
import argparse
import sys
from typing import Dict, List, Optional, Tuple, Set, Iterable
import pandas as pd


def normalize_gt_anyalt(gt_token: str) -> Optional[str]:
    if gt_token in (".", "./.", ".|."):
        return None
    tok = gt_token.replace("|", "/").strip()
    if tok == "":
        return None
    alleles = tok.split("/")
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


def is_sex_chrom(chrom: str) -> bool:
    if not isinstance(chrom, str):
        return False
    c = chrom.strip()
    c = c[3:] if c.lower().startswith("chr") else c
    return c.upper() in ("X", "Y")


def exonic_filter_ok(exonic_func: str) -> bool:
    if not isinstance(exonic_func, str):
        return False
    s = exonic_func.strip().lower()
    return s not in ("", ".", "synonymous snv", "unknown")


def func_refgene_has_exonic(func_refgene: str) -> bool:
    if not isinstance(func_refgene, str):
        return False
    return "exonic" in func_refgene.strip().lower()


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
        p_col, m_col, f_col = cols[fmt_idx + 1], cols[fmt_idx + 2], cols[fmt_idx + 3]
    else:
        p_col, m_col, f_col = cols[-3], cols[-2], cols[-1]
    return fmt_col, p_col, m_col, f_col


def load_gene_models(path: str) -> Tuple[List[Tuple[str, str]], Set[str], Set[str], Set[str]]:
    df = pd.read_csv(path, sep="\t", header=None, names=["Gene", "Model"], dtype=str)
    df["Gene"] = df["Gene"].str.strip()
    df["Model"] = df["Model"].str.strip().str.upper()
    pairs = list(df.itertuples(index=False, name=None))
    ad = set(df.loc[df["Model"] == "AD", "Gene"])
    ar = set(df.loc[df["Model"] == "AR", "Gene"])
    adar = set(df.loc[df["Model"] == "AD/AR", "Gene"])
    return pairs, ad, ar, adar


def gene_in_cell(target_gene: str, cell: str) -> bool:
    if not isinstance(cell, str) or cell.strip() in ("", "."):
        return False
    return any(tok.strip() == target_gene for tok in cell.split(";" ))


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


def per_gene_pass1_count(
    input_file: str,
    target_gene: str,
    fmt_c: str,
    p_c: str,
    m_c: str,
    f_c: str,
    af_col: str,
    af_threshold: float,
    chunksize: int,
    require_exonic: bool,
) -> int:
    gt_index_cache: Dict[str, int] = {}
    total = 0
    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)
    for chunk in reader:
        if "Gene.refGene" not in chunk.columns:
            continue
        chunk = chunk[chunk["Gene.refGene"].apply(lambda s: gene_in_cell(target_gene, s))]
        if chunk.empty:
            continue
        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        chunk = chunk[af_vals.isna() | (af_vals < af_threshold)]
        if chunk.empty:
            continue
        if require_exonic:
            if "Func.refGene" in chunk.columns:
                chunk = chunk[chunk["Func.refGene"].apply(func_refgene_has_exonic)]
                if chunk.empty:
                    continue
            if "ExonicFunc.refGene" in chunk.columns:
                chunk = chunk[chunk["ExonicFunc.refGene"].apply(exonic_filter_ok)]
                if chunk.empty:
                    continue
        if fmt_c not in chunk.columns or p_c not in chunk.columns:
            continue
        fmt_s = chunk[fmt_c].astype(str)
        p_s = chunk[p_c].astype(str)
        keep = 0
        for fmt, ps in zip(fmt_s, p_s):
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                continue
            p_norm = normalize_gt_anyalt(extract_gt_from_sample(ps, gt_idx))
            if p_norm is None:
                continue
            keep += 1
        total += keep
    return total


def per_gene_pass2_emit(
    input_file: str,
    target_gene: str,
    model: str,
    fmt_c: str,
    p_c: str,
    m_c: str,
    f_c: str,
    af_col: str,
    af_threshold: float,
    chunksize: int,
    require_exonic: bool,
    comphet_path: str,
    denovo_merr_path: str,
    comphet_header_written: List[bool],
    denovo_header_written: List[bool],
    skip_sex_chrom_for_comphet: bool = True,
):
    gt_index_cache: Dict[str, int] = {}
    gene_maternal: Set[Tuple] = set()
    gene_paternal: Set[Tuple] = set()
    row_store: Dict[Tuple, Dict[str, str]] = {}
    row_order: Dict[Tuple, int] = {}
    order_counter = 0

    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)
    for chunk in reader:
        if "Gene.refGene" not in chunk.columns:
            continue
        chunk = chunk[chunk["Gene.refGene"].apply(lambda s: gene_in_cell(target_gene, s))]
        if chunk.empty:
            continue
        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        chunk = chunk[af_vals.isna() | (af_vals < af_threshold)]
        if chunk.empty:
            continue
        if require_exonic:
            if "Func.refGene" in chunk.columns:
                chunk = chunk[chunk["Func.refGene"].apply(func_refgene_has_exonic)]
                if chunk.empty:
                    continue
            if "ExonicFunc.refGene" in chunk.columns:
                chunk = chunk[chunk["ExonicFunc.refGene"].apply(exonic_filter_ok)]
                if chunk.empty:
                    continue
        for col in (fmt_c, p_c, m_c, f_c, "Chr", "Start", "End", "Ref", "Alt"):
            if col not in chunk.columns:
                chunk = chunk.iloc[0:0]
                break
        if chunk.empty:
            continue

        fmt_s = chunk[fmt_c].astype(str)
        p_s = chunk[p_c].astype(str)
        m_s = chunk[m_c].astype(str)
        f_s = chunk[f_c].astype(str)

        trios: List[Tuple[Optional[str], Optional[str], Optional[str]]] = []
        for fmt, ps, ms, fs in zip(fmt_s, p_s, m_s, f_s):
            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                trios.append((None, None, None))
                continue
            p_norm = normalize_gt_anyalt(extract_gt_from_sample(ps, gt_idx))
            m_norm = normalize_gt_anyalt(extract_gt_from_sample(ms, gt_idx))
            f_norm = normalize_gt_anyalt(extract_gt_from_sample(fs, gt_idx))
            trios.append((p_norm, m_norm, f_norm))

        chunk = chunk.assign(_norm_trio=trios)
        full = chunk[chunk["_norm_trio"].apply(lambda t: all(gt in ("0/0", "0/1", "1/1") for gt in t))].copy()

        if not full.empty:
            def is_denovo(trio: Tuple[str, str, str]) -> bool:
                p, m, d = trio
                return (p in ("0/1", "1/1")) and m == "0/0" and d == "0/0"

            def is_merr(trio: Tuple[str, str, str]) -> bool:
                p, m, d = trio
                pc, mc, dc = gt_class(p), gt_class(m), gt_class(d)
                return not mendelian_consistent_anyalt(pc, mc, dc)

            denovo_mask = full["_norm_trio"].apply(is_denovo)
            merr_mask = full["_norm_trio"].apply(lambda t: is_merr(t) and not is_denovo(t))
            out_dm = full[denovo_mask | merr_mask].drop(columns=["_norm_trio"]).copy()
            if not out_dm.empty:
                out_dm.insert(0, "GenePanelHit", target_gene)
                out_dm.to_csv(
                    denovo_merr_path,
                    sep="\t",
                    mode="a",
                    header=(not denovo_header_written[0]),
                    index=False,
                )
                denovo_header_written[0] = True

        comp_src = full
        if not comp_src.empty:
            for _, row in comp_src.iterrows():
                if skip_sex_chrom_for_comphet and is_sex_chrom(row["Chr"]):
                    continue
                p_norm, m_norm, f_norm = row["_norm_trio"]
                if p_norm != "0/1":
                    continue
                m_alt = (m_norm in ("0/1", "1/1"))
                f_alt = (f_norm in ("0/1", "1/1"))
                m_non = (m_norm in ("0/0", None))
                f_non = (f_norm in ("0/0", None))
                tag_m = m_alt and f_non
                tag_p = f_alt and m_non
                if not (tag_m or tag_p):
                    continue
                rk = (row["Chr"], row["Start"], row["End"], row["Ref"], row["Alt"])
                if rk not in row_store:
                    base = row.drop(labels=["_norm_trio"]).to_dict()
                    base["GenePanelHits"] = target_gene
                    row_store[rk] = base
                    row_order[rk] = order_counter
                    order_counter += 1
                if tag_m:
                    gene_maternal.add(rk)
                if tag_p:
                    gene_paternal.add(rk)

    if model in ("AR", "AD/AR"):
        if gene_maternal and gene_paternal:
            rows = sorted(gene_maternal | gene_paternal, key=lambda k: row_order.get(k, 1_000_000_000_000))
            out = pd.DataFrame([row_store[rk] for rk in rows])
            if not out.empty:
                out.to_csv(
                    comphet_path,
                    sep="\t",
                    mode="a",
                    header=(not comphet_header_written[0]),
                    index=False,
                )
                comphet_header_written[0] = True


def main():
    ap = argparse.ArgumentParser(
        description="Memory-lean, gene-by-gene compound-het and de novo/MERR extraction."
    )
    ap.add_argument("-i", "--input", required=True, help="ANNOVAR TSV (trio multianno)")
    ap.add_argument("-g", "--genelist", required=True, help="Two-column TSV: Gene<TAB>Model (AD|AR|AD/AR)")
    ap.add_argument("-o", "--output-prefix", required=True, help="Creates <prefix>_comphet.tsv and <prefix>_denovo_merr.tsv")
    ap.add_argument("--format-col", default=None)
    ap.add_argument("--proband-col", default=None)
    ap.add_argument("--mother-col", default=None)
    ap.add_argument("--father-col", default=None)
    ap.add_argument("--af-col", default="gnomad41_genome_AF")
    ap.add_argument("--af-threshold", type=float, default=0.001)
    ap.add_argument("--chunksize", type=int, default=100000)
    ap.add_argument("--require-exonic", action="store_true")
    ap.add_argument("--no-skip-sex-for-comphet", dest="no_skip_sex_for_comphet", action="store_true")
    args = ap.parse_args()

    try:
        head_df = pd.read_csv(args.input, sep="\t", dtype=str, nrows=1, low_memory=False)
    except Exception as e:
        print(f"[ERROR] Could not read input header: {e}", file=sys.stderr)
        sys.exit(1)

    if all(v is not None for v in (args.format_col, args.proband_col, args.mother_col, args.father_col)):
        fmt_c, p_c, m_c, f_c = args.format_col, args.proband_col, args.mother_col, args.father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(list(head_df.columns))

    pairs, AD, AR, ADAR = load_gene_models(args.genelist)
    comphet_path = f"{args.output_prefix}_comphet.tsv"
    denovo_merr_path = f"{args.output_prefix}_denovo_merr.tsv"
    comphet_header_written = [False]
    denovo_header_written = [False]

    for gene, model in pairs:
        count = per_gene_pass1_count(
            input_file=args.input,
            target_gene=gene,
            fmt_c=fmt_c,
            p_c=p_c,
            m_c=m_c,
            f_c=f_c,
            af_col=args.af_col,
            af_threshold=args.af_threshold,
            chunksize=args.chunksize,
            require_exonic=args.require_exonic,
        )
        m = model.upper()
        if m == "AR":
            allow = (count >= 2)
        elif m in ("AD", "AD/AR"):
            allow = (count >= 1)
        else:
            allow = False
        if not allow:
            continue

        per_gene_pass2_emit(
            input_file=args.input,
            target_gene=gene,
            model=m,
            fmt_c=fmt_c,
            p_c=p_c,
            m_c=m_c,
            f_c=f_c,
            af_col=args.af_col,
            af_threshold=args.af_threshold,
            chunksize=args.chunksize,
            require_exonic=args.require_exonic,
            comphet_path=comphet_path,
            denovo_merr_path=denovo_merr_path,
            comphet_header_written=comphet_header_written,
            denovo_header_written=denovo_header_written,
            skip_sex_chrom_for_comphet=(not args.no_skip_sex_for_comphet),
        )

    print("Done.")
    print(f"  Comp-het: {comphet_path}")
    print(f"  De novo / MERR: {denovo_merr_path}")


if __name__ == "__main__":
    main()

