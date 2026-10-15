"""Orchestration: build patterns, scan references, dispatch BAM workers."""
import os
from multiprocessing import Pool
from typing import Set

from .automaton import build_exact, build_flexible
from .bam_scan import process_bam
from .patterns import FLEX_MIN_LEN, PatternSets
from .references import (MetamorfScanner, NuorfScanner,
                         collapse_overlapping_loci, merge_reference_results)


class Pipeline:
    def __init__(self, args, flex_min_len: int = FLEX_MIN_LEN):
        self.args = args
        self.flex_min_len = flex_min_len

    def scan_references(self, patterns: PatternSets):
        scans = []

        if patterns.exact:
            A = build_exact(patterns.exact)
            scans += [
                NuorfScanner(A, patterns.exact).scan(self.args.nuorf),
                MetamorfScanner(A, patterns.exact).scan(self.args.metamorf),
            ]

        if patterns.flexible:
            A = build_flexible(patterns.flexible, self.flex_min_len)
            scans += [
                NuorfScanner(A, patterns.flexible).scan(self.args.nuorf),
                MetamorfScanner(A, patterns.flexible).scan(self.args.metamorf),
            ]

        peptides = merge_reference_results(*scans)
        collapse_overlapping_loci(peptides)
        return peptides

    def scan_bams(self, peptides, flexible_set: Set[str]):
        bam_files = [x for x in os.listdir(self.args.bam_dir)
                     if x.endswith("bam")]
        tasks = [(bf, self.args.bam_dir, peptides, self.args.min_mapq,
                  flexible_set, self.flex_min_len) for bf in bam_files]
        total_hits, rows = 0, []
        with Pool(processes=self.args.workers) as pool:
            for _, hits, worker_rows in pool.starmap(process_bam, tasks,
                                                     chunksize=1):
                total_hits += hits
                rows.extend(worker_rows)
        return bam_files, total_hits, rows
