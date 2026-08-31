#!/usr/bin/env python3
"""
Telomere Profiler
-----------------
Publication-ready pipeline for telomeric motif detection at chromosome ends.
Author: (your name)
"""

import argparse
import os
import re
import sys
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq

# ----------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------
def find_overlapping(seq: str, motif: str):
    """Return 0-based start coordinates of all overlapping hits."""
    return [m.start() for m in re.finditer(f"(?={re.escape(motif)})", seq)]


def longest_tract(positions, motif_len):
    """Length (bp) of the longest perfect tandem run."""
    if not positions:
        return 0
    best = cur = 1
    for i in range(1, len(positions)):
        if positions[i] - positions[i - 1] == motif_len:
            cur += 1
        else:
            best = max(best, cur)
            cur = 1
    return max(best, cur) * motif_len


def rc_positions(fwd_len, rel_pos, motif_len):
    """Convert RC-relative coordinates back to forward-strand coordinates."""
    return [fwd_len - p - motif_len for p in rel_pos]


# ----------------------------------------------------------------------
# main
# ----------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Telomere motif analysis – publication pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Genome FASTA")
    parser.add_argument(
        "-m", "--motifs", default="TTTAGGG,CCCTAAA",
        help="Comma-separated motifs",
    )
    parser.add_argument("-o", "--outdir", default="Telomere_Output")
    parser.add_argument(
        "-w", "--window", type=int, default=10000,
        help="bp taken from each chromosome end",
    )
    parser.add_argument(
        "--rc", action="store_true",
        help="Also search reverse-complement of each end",
    )
    parser.add_argument("--min-len", type=int, default=1000,
                        help="Skip chromosomes shorter than this")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "fasta").mkdir(exist_ok=True)

    logging.basicConfig(
        filename=outdir / "pipeline.log",
        level=logging.INFO,
        format="%(asctime)s  %(levelname)s  %(message)s",
    )
    log = logging.getLogger()

    motifs = [m.strip().upper() for m in args.motifs.split(",") if m.strip()]
    log.info(f"Motifs: {motifs}  |  window: {args.window} bp  |  RC: {args.rc}")

    counts_rows = []
    pos_rows = []
    chrom_len = {}

    start_fh = open(outdir / "fasta" / f"starts_{args.window}bp.fa", "w")
    end_fh   = open(outdir / "fasta" / f"ends_{args.window}bp.fa", "w")

    for rec in SeqIO.parse(args.input, "fasta"):
        chrom = rec.id
        seq = str(rec.seq).upper()
        L = len(seq)
        chrom_len[chrom] = L

        if L < args.min_len:
            log.warning(f"{chrom} too short ({L} bp) – skipped")
            continue

        w = min(args.window, L)
        start_seq = seq[:w]
        end_seq   = seq[-w:]
        end_off   = L - w          # 0-based offset of the end window

        start_fh.write(f">{chrom}_start_{w}bp\n{start_seq}\n")
        end_fh.write(f">{chrom}_end_{w}bp\n{end_seq}\n")

        regions = [
            ("Start", start_seq, 0, "+"),
            ("End",   end_seq,   end_off, "+"),
        ]
        if args.rc:
            regions += [
                ("Start", str(Seq(start_seq).reverse_complement()), 0, "-"),
                ("End",   str(Seq(end_seq).reverse_complement()),   end_off, "-"),
            ]

        for region_name, region_seq, offset, strand in regions:
            for motif in motifs:
                rel = find_overlapping(region_seq, motif)

                # convert RC hits back to forward coordinates
                if strand == "-":
                    rel = rc_positions(len(region_seq), rel, len(motif))
                    rel.sort()

                abs_pos = [p + offset + 1 for p in rel]          # 1-based
                tract   = longest_tract(rel, len(motif))

                counts_rows.append({
                    "Chromosome": chrom,
                    "Length": L,
                    "Region": region_name,
                    "Strand": strand,
                    "Motif": motif,
                    "Count": len(rel),
                    "Longest_Tract_bp": tract,
                })

                for r, a in zip(rel, abs_pos):
                    pos_rows.append({
                        "Chromosome": chrom,
                        "Region": region_name,
                        "Strand": strand,
                        "Motif": motif,
                        "Relative_Position": r + 1,
                        "Absolute_Position": a,
                    })

    start_fh.close()
    end_fh.close()

    # ------------------------------------------------------------------
    # tables
    # ------------------------------------------------------------------
    df_counts = pd.DataFrame(counts_rows)
    df_pos    = pd.DataFrame(pos_rows)

    df_counts.to_csv(outdir / "telomere_summary.tsv", sep="\t", index=False)
    df_pos.to_csv(outdir / "motif_positions.tsv", sep="\t", index=False)

    with pd.ExcelWriter(outdir / "Telomere_Supplementary.xlsx") as xls:
        df_counts.to_excel(xls, sheet_name="Summary", index=False)
        df_pos.to_excel(xls, sheet_name="Positions", index=False)

    # nice console preview
    print("\n=== Count summary (first 20 rows) ===")
    print(df_counts.head(20).to_string(index=False))
    print(f"\nFull tables written to {outdir}/")

    # ------------------------------------------------------------------
    # plots (300 dpi, paper-ready)
    # ------------------------------------------------------------------
    sns.set_theme(style="ticks", font_scale=1.1)

    # 1. grouped barplot
    if not df_counts.empty:
        plt.figure(figsize=(12, 6))
        sns.barplot(
            data=df_counts[df_counts.Strand == "+"],
            x="Chromosome", y="Count", hue="Region",
            palette="Set2", errorbar=None,
        )
        plt.xticks(rotation=45, ha="right")
        plt.title(f"Telomeric motif counts (terminal {args.window} bp)")
        plt.tight_layout()
        plt.savefig(outdir / "plots" / "01_counts_barplot.png", dpi=300)
        plt.close()

        # 2. heatmap of longest tract
        pivot = df_counts[df_counts.Strand == "+"].pivot_table(
            index="Chromosome",
            columns=["Region", "Motif"],
            values="Longest_Tract_bp",
            fill_value=0,
        )
        plt.figure(figsize=(10, max(6, len(pivot) * 0.4)))
        sns.heatmap(pivot, annot=True, fmt=".0f", cmap="YlOrRd",
                    cbar_kws={"label": "Longest perfect tract (bp)"})
        plt.title("Longest tandem telomeric tract")
        plt.tight_layout()
        plt.savefig(outdir / "plots" / "02_tract_heatmap.png", dpi=300)
        plt.close()

        # 3. chromosomal landscape
        if not df_pos.empty:
            plt.figure(figsize=(12, max(4, len(chrom_len) * 0.35)))
            y_labels = list(chrom_len.keys())
            for i, c in enumerate(y_labels):
                plt.plot([0, chrom_len[c]], [i, i], color="0.8", lw=6, zorder=1)

            sns.scatterplot(
                data=df_pos[df_pos.Strand == "+"],
                x="Absolute_Position", y="Chromosome",
                hue="Motif", style="Region", s=55, alpha=0.85, zorder=2,
            )
            plt.xlabel("Genomic position (bp)")
            plt.title("Location of telomeric motifs")
            plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
            plt.tight_layout()
            plt.savefig(outdir / "plots" / "03_chromosome_map.png", dpi=300)
            plt.close()

    log.info("Finished successfully")
    print(f"\n[✓] Done. Results in  {outdir.resolve()}")


if __name__ == "__main__":
    main()
