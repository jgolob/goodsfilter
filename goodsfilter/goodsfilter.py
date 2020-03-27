#!/usr/bin/env python
import csv
from collections import defaultdict
import random
import os
import logging
import argparse

logging.basicConfig(format="Good's Filter:%(levelname)s:%(asctime)s:%(message)s", level=logging.INFO)

def main():
    args_parser = argparse.ArgumentParser(description="""Filter DADA2-style sequence tables using
    Good's coverage saturation to remove likely contaminant and PCR error sequence variants
    """)

    args_parser.add_argument('--seqtable', '-s',
                             help="Sequence table from dada2, in CSV format",
                             required=True, type=argparse.FileType('rt'))
    args_parser.add_argument('--seqtable_filtered', '-O',
                             help="Write filtered dada2 style seqtable to this file (CSV)",
                             type=argparse.FileType('wt'),
                             required=True,
    )
    args_parser.add_argument('--converged_file', '-C',
                             help="Write if Good's coverged for each specimen, CSV format. (indicates of there was adequate read depth)",
                             type=argparse.FileType('wt'),
    )
    args_parser.add_argument('--keep_nonconverged', '-KNC',
                            help="Keep specimens whose Good's coverage failed to converge (ie. of insufficient depth)",
                            action='store_true'
    )
    args_parser.add_argument('--iteration_cutoff', '-IC',
                            help="Good's iteration cutoff (default=0.0001)",
                            default=0.0001
    )
    args_parser.add_argument('--min_prev', '-MP',
                            help="Minimum Prevalence of a SV to be retained (default=1, no filter)",
                            default=1
    )
    args_parser.add_argument('--min_reads', '-MR',
                            help="Minimum Reads before triggering Good's filter (default=10)",
                            default=10
    )
    args_parser.add_argument('--curves_path', '-cp',
                             help="Path to write the collector's curves in path/(specimen)_collector.csv",
    )
    

    args = args_parser.parse_args()

    # Dict with specimen as a key, and count list as a value
    filtered_counts = {}
    rareifaction_curves = {}
    goods_converged = {}

    dada2_r = csv.reader(args.seqtable)
    sv = next(dada2_r)[1:]
    # Each specimen is in a row
    # For each specimen, extract the counts, make a virual read file (based on the counts)
    # shuffle those reads, and make a collector's curve
    for row in dada2_r:
        specimen = row[0]
        sp_counts = [int(v) for v in row[1:]]
        # Make virual reads
        sp_reads = [
            item for sublist 
            in [ [sv[c_i]]*c for c_i, c in enumerate(sp_counts) ]
            for item in sublist
        ]
        # Shuffle
        random.shuffle(sp_reads)
        
        # Here is where we make our curve
        sp_esv_count = defaultdict(int)
        sp_rareifaction_curve = []
        sp_goods_threshold = None

        for sr_i, sr in enumerate(sp_reads):
            sp_esv_count[sr] += 1
            sp_total_reads = sr_i + 1
            sp_esv_total = len(sp_esv_count)
            sp_esv_singleton = len([
                v for v in 
                sp_esv_count.values()
                if v == 1
            ])
            sp_rareifaction_curve.append([
                sp_total_reads,
                sp_esv_total,
                sp_esv_singleton,
                1.0 - float(sp_esv_singleton) / sp_total_reads
            ])
            if (sp_goods_threshold is None) and sp_total_reads > int(args.min_reads) and abs(sp_rareifaction_curve[-1][3] - sp_rareifaction_curve[-2][3]) <= float(args.iteration_cutoff):
                sp_goods_threshold = sp_total_reads        
        # Let's use our threshold to filter
        
        if sp_goods_threshold is not None:
            goods_converged[specimen] = True
            filtered_counts[specimen] = [c if c >= sp_goods_threshold else 0 for c in sp_counts]
        else:
            goods_converged[specimen] = False
            filtered_counts[specimen] = sp_counts
        
        rareifaction_curves[specimen] = sp_rareifaction_curve

    # Filter out globally unique (SV only found in less than MIN_PREVALENCE specimen)
    sv_pass_prev_filter = [
        len([
            s_counts[sv_i]
            for s_i, s_counts in enumerate(filtered_counts.values())
            if s_counts[sv_i] > 0 and list(goods_converged.values())[s_i]
        ]) >= int(args.min_prev)
        for sv_i in range((len(sv)))
    ]

     # Integrate all this together and output the filtered seqtab
    filtered_header = [""]+[sv for sv_i, sv in enumerate(sv) if sv_pass_prev_filter[sv_i]]
    out_w = csv.writer(args.seqtable_filtered)
    out_w.writerow(filtered_header)
    if args.keep_nonconverged:
        out_w.writerows([   
            [sp]+[c for c_i, c in enumerate(sp_filt_count) if sv_pass_prev_filter[c_i]]
            for sp, sp_filt_count in filtered_counts.items()
        ])
    else:
        out_w.writerows([   
            [sp]+[c for c_i, c in enumerate(sp_filt_count) if sv_pass_prev_filter[c_i]]
            for sp, sp_filt_count in filtered_counts.items()
            if goods_converged[sp] is True
        ])
        

    if args.converged_file is not None:
        converged_w = csv.writer(args.converged_file)
        converged_w.writerow([
            'specimen',
            'goods_converged'
        ])
        converged_w.writerows([
            r for r in goods_converged.items()
        ])
    
    # Output collector's curves
    if args.curves_path is not None:
        try:
            os.makedirs(args.curves_path)
        except:
            pass
        for specimen, curves in rareifaction_curves.items():
            with open(os.path.join(args.curves_path, "{}_collector.csv".format(specimen)), 'wt') as sp_collect_h:
                sp_collect_w = csv.writer(sp_collect_h)
                sp_collect_w.writerow([
                    'read_num',
                    'esv_total',
                    'esv_singleton',
                    'goods_coverage'
                ])
                sp_collect_w.writerows(
                    curves
                )

# Boilerplate method to run this as a script
if __name__ == '__main__':
    main()