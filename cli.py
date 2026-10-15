"""Command-line interface."""
import argparse
import sys

from .patterns import FLEX_MIN_LEN, read_id_list, read_peptides_tsv, split_patterns
from .pipeline import Pipeline


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        description="Search RNA-seq BAMs for reads whose translation "
                    "contains any of a set of target peptides.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--peptides", required=True,
                   help="TSV file: peptide <tab> identifier (with header).")
    p.add_argument("--nuorf", required=True,
                   help="nuORFdb protein FASTA.")
    p.add_argument("--metamorf", required=True,
                   help="MetamORF protein FASTA.")
    p.add_argument("--bam-dir", required=True,
                   help="Directory containing .bam files.")

    p.add_argument("--flexible-ids", default=None,
                   help="Optional file listing ids that should use flexible "
                        "(substring, min-length) matching. Ids not listed "
                        "use strict exact matching. If omitted, all ids use "
                        "exact matching.")
    p.add_argument("--flex-min-len", type=int, default=FLEX_MIN_LEN,
                   help="Minimum substring length for flexible matching.")

    p.add_argument("--output", default="hits.tsv",
                   help="Output TSV path.")
    p.add_argument("--workers", type=int, default=28,
                   help="Number of parallel worker processes.")
    p.add_argument("--min-mapq", type=int, default=0,
                   help="Minimum mapping quality for a read to be kept.")
    p.add_argument("--quiet", action="store_true",
                   help="Do not echo the report to stdout.")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    flexible_ids = (read_id_list(args.flexible_ids)
                    if args.flexible_ids else None)

    pattern_map = read_peptides_tsv(args.peptides)
    if not pattern_map:
        print("[error] no peptides found in --peptides", file=sys.stderr)
        sys.exit(1)

    patterns = split_patterns(pattern_map, flexible_ids=flexible_ids)

    print(f"[info] {len(patterns.exact)} exact, "
          f"{len(patterns.flexible)} flexible peptide patterns",
          file=sys.stderr)

    pipe = Pipeline(args, flex_min_len=args.flex_min_len)
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
