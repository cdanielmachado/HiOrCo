#!/usr/bin/env python

import argparse
import textwrap
import pandas as pd

from hiorco.method import compute

def main():

    parser = argparse.ArgumentParser(description="Compute higher-order sets of co-occurring species.",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('samples', help=textwrap.dedent(
        """
        Species abundance table.
        
        Should be a tab (default) or comma separated file with species as rows and samples as columns. 
        """
        ))

    parser.add_argument('--sep', default='\t', help="Separator character (default: '\t')." )
    parser.add_argument('-k', type=int, default=10, help="Compute species sets up to size k (default: 10).")
    parser.add_argument('-n', type=int, default=100, help="Store up to n species sets per size (default: 100).")
    parser.add_argument('-p', type=int, default=1000, help="Solution pool size to carry over to next iteration (default: 1000).")
    parser.add_argument('--cutoff', type=float, default=0.001, help="Relative abundance cutoff value (default: 0.001).")
    parser.add_argument('-o', '--output', dest="output", default=".", help="Output folder (optional).")
    parser.add_argument('--part-size', type=int, default=100, help="Split output files with a maximum number of species sets per file (default: 100).")
    parser.add_argument('--single', action='store_true', help="Disable multi-threaded mode.")

    args = parser.parse_args()

    try:
        data = pd.read_csv(args.samples, sep=args.sep, index_col=0)
    except Exception as e:
        print(e)
        exit()

    compute(data, k=args.k, n=args.n, p=args.p, abundance_cutoff=args.cutoff,
        output_folder=args.output, part_size=args.part_size, parallel=(not args.single))


if __name__ == '__main__':
    main()