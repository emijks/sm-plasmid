#! /usr/bin/env python3

import argparse
import snakemake


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-x', '--table', required=True)
    argparser.add_argument('-s', '--seq-dir', required=True)
    argparser.add_argument('-r', '--ref-dir', required=True)
    argparser.add_argument('-o', '--out-dir', required=True)
    argparser.add_argument('-t', '--tmp-dir', default='./tmp')
    argparser.add_argument('-@', '--threads', type=int, default=4)
    argparser.add_argument('-f', '--seq-format', type=int, default=1, help='<1|2|3|4>: xlsx table to seq-fpaths format')
    argparser.add_argument('-d', '--dry-run', action='count', default=0)
    args = argparser.parse_args()

    config = vars(args)

    snakemake.snakemake(
        snakefile='Snakefile',
        config=config,
        cores=config['threads'],
        use_singularity=True,
        singularity_args="-B $HOME -B /primary",
        dryrun=bool(args.dry_run),
        quiet=True,
        keepgoing=True,
        keep_incomplete=True,
        force_incomplete=True,
        # shadow_prefix='/primary/home/emijks/dev/sm-plasmid-stuff'
        shadow_prefix='.snakemake'
    )
