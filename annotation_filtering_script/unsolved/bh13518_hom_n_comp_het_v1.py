#!/usr/bin/env python3
import argparse
import sys
from typing import Dict, List, Optional, Tuple, Set
import pandas as pd

# =========================
# Genotype normalization
# =========================
def normalize_gt_anyalt(gt_token: str) -> Optional[str]:
    """
    Collapse GT to '0/0', '0/1', or '1/1' using the any-ALT rule:
      - '|' treated as '/'
      - any non-zero allele index -> '1'
      - haploid 'a' -> 'a/a'
      - missing '.', './.', '.|.' -> None
    """
    if gt_token in (".", "./.", ".|."):
        return None
    tok = gt_token.replace("|", "/").strip()
    if tok == "":
        return None
    alleles = tok.split("/")
    if len(alleles) == 1:  # haploid -> duplicate
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

# =========================
# FORMAT / GT helpers
# =========================
def get_gt_index(fmt: str, cache: Dict[str, int]) -> int:
    """Return index of 'GT' within FORMAT, caching by FORMAT value. -1 if absent/invalid."""
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
    """Extract GT token from a sample column using GT index derived from FORMAT."""
    if not isinstance(sample_field, str) or sample_field == "":
        return ""
    parts = sample_field.split(":")
    if gt_idx < 0 or gt_idx >= len(parts):
        # Fallback: some callers keep GT first
        return parts[0] if parts else ""
    return parts[gt_idx]

# =========================
# Column autodetection
# =========================
def autodetect_columns(df_columns: List[str]) -> Tuple[str, str, str, str]:
    """
    Find FORMAT + trio sample columns in ANNOVAR TSVs.
    Prefers Otherinfo12..15, then Otherinfo9..12, else guesses FORMAT-like column and next 3.
    """
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
        fmt_idx = max(0, len(cols) - 6)  # heuristic

    fmt_col = cols[fmt_idx]
    if fmt_idx + 3 < len(cols):
        p_col, m_col, f_col = cols[fmt_idx + 1], cols[fmt_idx + 2], cols[fmt_idx + 3]
    else:
        p_col, m_col, f_col = cols[-3], cols[-2], cols[-1]
    return fmt_col, p_col, m_col, f_col

# =========================
# Utilities
# =========================
def parse_gene_cell(cell: str) -> List[str]:
    if not isinstance(cell, str) or cell.strip() in ("", "."):
        return []
    # ANNOVAR uses ';' between gene symbols when multiple overlap
    parts = [g.strip() for g in cell.split(";")]
    return [g for g in parts if g and g != "."]

def make_row_key_from_row(row: pd.Series) -> Tuple:
    return (row.get("Chr"), row.get("Start"), row.get("End"), row.get("Ref"), row.get("Alt"))

# =========================
# Main processing
# =========================
def process_biallelic_madd_panel(
    input_file: str,
    output_prefix: str,
    genes_of_interest: List[str],
    format_col: Optional[str] = None,
    proband_col: Optional[str] = None,
    mother_col: Optional[str] = None,
    father_col: Optional[str] = None,
    af_col: str = "gnomad41_genome_AF",
    af_threshold: float = 0.001,
    chunksize: int = 100000,
):
    """
    BH13518-1: Search for biallelic variants in ETFA, ETFB, ETFDH, FLAD1, SLC52A1, SLC52A2, SLC52A3
      1) Population AF < threshold (or missing)
      2) Variant overlaps interested gene(s) AND all trio samples have genotypes (no missing)
      3.1) Homozygous: proband 1/1, BOTH parents 0/1  -> emit single-variant rows
      3.2) Compound-het (lenient): proband 0/1, and per-variant inheritance tags:
              - maternal-tag: mother in {0/1,1/1} AND father == 0/0
              - paternal-tag: father in {0/1,1/1} AND mother == 0/0
           Gene qualifies if it has at least one maternal-tagged row AND one paternal-tagged row.
           (No requirement on Func.refGene being 'exonic'.)
    Outputs:
      <prefix>_hom.tsv
      <prefix>_comphet.tsv
    """
    hom_out = f"{output_prefix}_hom.tsv"
    comphet_out = f"{output_prefix}_comphet.tsv"

    panel = {g.upper() for g in genes_of_interest}

    # Global header and storage
    header_cols: List[str] = []
    gt_index_cache: Dict[str, int] = {}

    # Homozygous candidates: set of row keys
    hom_rows: Set[Tuple] = set()
    # For comphet: gene -> sets of row_keys
    gene_maternal: Dict[str, Set[Tuple]] = {}
    gene_paternal: Dict[str, Set[Tuple]] = {}

    # Row key -> original row dict; and input order
    row_store: Dict[Tuple, Dict[str, str]] = {}
    row_order: Dict[Tuple, int] = {}
    order_counter = 0

    reader = pd.read_csv(input_file, sep="\t", dtype=str, chunksize=chunksize, low_memory=False)

    # consume first chunk to set columns and optionally detect trio cols
    try:
        first_chunk = next(reader)
    except StopIteration:
        print("[INFO] Empty input.", file=sys.stderr)
        # still emit empty headers
        pd.DataFrame(columns=[]).to_csv(hom_out, sep="\t", index=False)
        pd.DataFrame(columns=[]).to_csv(comphet_out, sep="\t", index=False)
        return

    header_cols = list(first_chunk.columns)

    if all(col is not None for col in (format_col, proband_col, mother_col, father_col)):
        fmt_c, p_c, m_c, f_c = format_col, proband_col, mother_col, father_col
    else:
        fmt_c, p_c, m_c, f_c = autodetect_columns(header_cols)

    # Helper to push a row into store once
    def ensure_store(row: pd.Series) -> Tuple:
        nonlocal order_counter
        rk = make_row_key_from_row(row)
        if rk not in row_store:
            row_store[rk] = {col: (row[col] if col in row else "") for col in header_cols}
            row_order[rk] = order_counter
            order_counter += 1
        return rk

    # iterate chunks
    for chunk in [first_chunk] + list(reader):
        if chunk.empty:
            continue

        # AF filter (missing passes)
        af_vals = pd.to_numeric(chunk.get(af_col, pd.Series([None] * len(chunk))), errors="coerce")
        chunk = chunk[af_vals.isna() | (af_vals < af_threshold)]
        if chunk.empty:
            continue

        # Must have needed columns
        needed = [fmt_c, p_c, m_c, f_c, "Gene.refGene", "Chr", "Start", "End", "Ref", "Alt"]
        if any(col not in chunk.columns for col in needed):
            # If missing, skip this chunk
            continue

        # Restrict to interested genes (Gene.refGene overlaps)
        genes_series = chunk["Gene.refGene"].astype(str)
        keep_mask = genes_series.apply(
            lambda s: any((g.strip().upper() in panel) for g in parse_gene_cell(s))
        )
        chunk = chunk[keep_mask]
        if chunk.empty:
            continue

        # Trio columns
        fmt_s = chunk[fmt_c].astype(str)
        p_s   = chunk[p_c].astype(str)
        m_s   = chunk[m_c].astype(str)
        f_s   = chunk[f_c].astype(str)

        # Row-wise pass
        for i, (fmt, ps, ms, fs) in enumerate(zip(fmt_s, p_s, m_s, f_s)):
            row = chunk.iloc[i]

            gt_idx = get_gt_index(fmt, gt_index_cache)
            if gt_idx == -1:
                continue

            p_raw = extract_gt_from_sample(ps, gt_idx)
            m_raw = extract_gt_from_sample(ms, gt_idx)
            f_raw = extract_gt_from_sample(fs, gt_idx)

            p_norm = normalize_gt_anyalt(p_raw)
            m_norm = normalize_gt_anyalt(m_raw)
            f_norm = normalize_gt_anyalt(f_raw)

            # All samples must be genotyped
            if p_norm is None or m_norm is None or f_norm is None:
                continue

            # parse genes
            genes = parse_gene_cell(row.get("Gene.refGene"))
            if not genes:
                continue

            # 3.1 Homozygous: proband 1/1; BOTH parents 0/1
            if p_norm == "1/1" and m_norm == "0/1" and f_norm == "0/1":
                rk = ensure_store(row)
                hom_rows.add(rk)

            # 3.2 Compound-het (lenient): proband 0/1 with parental tagging
            if p_norm == "0/1":
                m_alt = (m_norm in ("0/1", "1/1"))
                f_alt = (f_norm in ("0/1", "1/1"))
                m_noncarrier = (m_norm == "0/0")   # NOTE: require true 0/0 (no missing allowed)
                f_noncarrier = (f_norm == "0/0")

                tag_maternal = m_alt and f_noncarrier
                tag_paternal = f_alt and m_noncarrier

                if not (tag_maternal or tag_paternal):
                    continue

                rk = ensure_store(row)
                for g in genes:
                    if tag_maternal:
                        gene_maternal.setdefault(g, set()).add(rk)
                    if tag_paternal:
                        gene_paternal.setdefault(g, set()).add(rk)

    # Collect outputs
    # Homozygous rows
    hom_out_rows = [row_store[rk] for rk in sorted(hom_rows, key=lambda k: row_order.get(k, 10**18))]

    # Compound-het: genes with both sides
    qualifying_genes = set(
        g for g in gene_maternal
        if gene_maternal[g] and g in gene_paternal and gene_paternal[g]
    )
    comphet_rows: Set[Tuple] = set()
    for g in qualifying_genes:
        comphet_rows.update(gene_maternal[g])
        comphet_rows.update(gene_paternal[g])
    comphet_out_rows = [row_store[rk] for rk in sorted(comphet_rows, key=lambda k: row_order.get(k, 10**18))]

    # Emit with original header order
    if header_cols:
        pd.DataFrame(hom_out_rows, columns=header_cols).to_csv(hom_out, sep="\t", index=False)
        pd.DataFrame(comphet_out_rows, columns=header_cols).to_csv(comphet_out, sep="\t", index=False)
    else:
        pd.DataFrame(hom_out_rows).to_csv(hom_out, sep="\t", index=False)
        pd.DataFrame(comphet_out_rows).to_csv(comphet_out, sep="\t", index=False)

    # Summary
    print("Biallelic MADD-panel scan complete.")
    print(f"  Panel genes: {','.join(sorted(panel))}")
    print(f"  Homozygous candidate variants: {len(hom_rows)} -> {hom_out}")
    print(f"  Genes meeting lenient comp-het (both sides): {len(qualifying_genes)}")
    if qualifying_genes:
        print("   ", ", ".join(sorted(qualifying_genes)))
    print(f"  Compound-het contributing variants: {len(comphet_rows)} -> {comphet_out}")

def main():
    ap = argparse.ArgumentParser(
        description="Find biallelic candidates (homozygous + compound-het) in MADD-related genes from ANNOVAR trio TSV."
    )
    ap.add_argument("-i", "--input", required=True,
                    help="Input ANNOVAR TSV (from merged trio VCF; includes FORMAT and trio sample columns)")
    ap.add_argument("-o", "--output-prefix", required=True,
                    help="Output prefix (creates <prefix>_hom.tsv and <prefix>_comphet.tsv)")
    ap.add_argument("--format-col", default=None,
                    help="FORMAT column name (auto-detected if omitted; e.g., Otherinfo12)")
    ap.add_argument("--proband-col", default=None,
                    help="Proband sample column (auto-detected if omitted; e.g., Otherinfo13)")
    ap.add_argument("--mother-col", default=None,
                    help="Mother sample column (auto-detected if omitted; e.g., Otherinfo14)")
    ap.add_argument("--father-col", default=None,
                    help="Father sample column (auto-detected if omitted; e.g., Otherinfo15)")
    ap.add_argument("--af-col", default="gnomad41_genome_AF",
                    help="Allele frequency column (default: gnomad41_genome_AF)")
    ap.add_argument("--af-threshold", type=float, default=0.001,
                    help="Keep if AF is missing OR < threshold (default: 0.001)")
    ap.add_argument("--genes", nargs="+",
                    default=["ETFA", "ETFB", "ETFDH", "FLAD1", "SLC52A1", "SLC52A2", "SLC52A3"],
                    help="Gene symbols to include (default: MADD-related panel)")
    ap.add_argument("--chunksize", type=int, default=100000,
                    help="Chunk size for large files (default: 100000)")
    args = ap.parse_args()

    process_biallelic_madd_panel(
        input_file=args.input,
        output_prefix=args.output_prefix,
        genes_of_interest=args.genes,
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

