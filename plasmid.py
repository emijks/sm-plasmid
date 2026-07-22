#! /usr/bin/env python3

import argparse
import snakemake
import subprocess
from os import path
import glob


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-x', '--table', required=True)
    argparser.add_argument('-s', '--seq-dir', required=True)
    argparser.add_argument('-r', '--ref-dir', required=True)
    argparser.add_argument('-o', '--out-dir', required=True)
    argparser.add_argument('-t', '--tmp-dir', default='./tmp')
    argparser.add_argument('-@', '--threads', type=int, default=4)
    argparser.add_argument('-f', '--seq-format', type=int, default=1, help='heeelp')
    argparser.add_argument('-d', '--dry-run', action='count', default=0)
    args = argparser.parse_args()

    config = vars(args)


    try:
        sucess = snakemake.snakemake(
        snakefile='Snakefile',
        config=config,
        cores=config['threads'],
        use_singularity=True,
        singularity_args=f"-B {args.seq_dir}:{args.seq_dir}:ro\
                           -B {args.ref_dir}:{args.ref_dir}:ro\
                           -B {args.tmp_dir}:{args.tmp_dir}:rw\
                           -B {args.out_dir}:{args.out_dir}:rw\
                           -B {args.table}:{args.table}:ro",
        dryrun=bool(args.dry_run),
        keepgoing=True,
        keep_incomplete=True,
        force_incomplete=True,
        #printshellcmds=True,
        quiet=True
    )
    finally:
        if not args.dry_run:
            subprocess.run(
                [
                    "python3",
                    "scripts/make_status_table.py",
                    "--table",
                    args.table,
                    "--seq-dir",
                    args.seq_dir,
                    "--ref-dir",
                    args.ref_dir,
                    "--seq-format",
                    str(args.seq_format),
                    "--out-dir",
                    args.out_dir,
                    "--output",
                    f"{args.out_dir}/pipeline_status.tsv"
                ],
                check=True    
            )

