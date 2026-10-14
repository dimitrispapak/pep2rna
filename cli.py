"""Command-line interface."""
import argparse
import sys

from .patterns import (DEFAULT_SHORT_IDS, LONG_IDS,
                       read_peptides_tsv, split_short_long)
from .pipeline import Pipeline

def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Search RNA-seq BAMs for reads whose translation "
                    "contains any of a set of target peptides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--peptides", required=True)
    p.add_argument("--nuorf", required=True)
    p.add_argument("--metamorf", required=True)
    p.add_argument("--bam-dir", required=True)
    p.add_argument("--pattern-ids", default=None,
                   help="Comma-separated subset of ids to keep.")
    p.add_argument("--output", default="hits.tsv")
    p.add_argument("--workers", type=int, default=28)
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--quiet", action="store_true")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    keep_ids = (set(args.pattern_ids.split(","))
                if args.pattern_ids else DEFAULT_SHORT_IDS)

    pattern_map = read_peptides_tsv(args.peptides, keep_ids)
    if not pattern_map:
        print("[error] no peptides matched the requested ids", file=sys.stderr)
        sys.exit(1)

    patterns = split_short_long(pattern_map)
    print(f"[info] {len(patterns.exact)} short (exact) patterns, "
          f"{len(patterns.flexible)} long (flexible) patterns",
          file=sys.stderr)

    pipe = Pipeline(args)
    peptides = pipe.scan_references(patterns)
    n_with_loci = sum(1 for info in peptides.values() if info.loci)
    print(f"[info] {n_with_loci} / {len(peptides)} peptides have loci",
          file=sys.stderr)

    flexible_set = set(patterns.flexible)
    bam_files, total_hits, rows = pipe.scan_bams(peptides, flexible_set)
    print(f"[info] scanned {len(bam_files)} BAM files", file=sys.stderr)

    with open(args.output, "w") as f:
        f.write("bam_file\tpeptide\tchrom\tstart\tend\tstrand\treads\thits\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")

    if not args.quiet:
        for row in rows:
            print("\t".join(map(str, row)))

    print(f"Total Hits: {total_hits}", file=sys.stderr)
    print(f"[info] wrote {len(rows)} rows to {args.output}", file=sys.stderr)
