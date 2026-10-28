#!/usr/bin/env python3

import argparse
import subprocess
import numpy as np


parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True)
parser.add_argument("--sample", required=True)
parser.add_argument("--output", required=True)
args = parser.parse_args()


e_score = "NA"
mean_depth = "NA"
status = "FAILED"


try:
    result = subprocess.run(
        [
            "samtools",
            "depth",
            "-aa",
            "-d",
            "100000",
            args.bam
        ],
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        status = "SAMTOOLS_FAILED"

    else:
        cov = np.array(
            [
                int(line.split()[2])
                for line in result.stdout.splitlines()
            ]
        )

        if cov.size == 0:
            status = "NO_DEPTH"

        else:
            mean_cov = cov.mean()
            mean_depth = f"{mean_cov:.2f}"

            if mean_cov == 0:
                status = "ZERO_COVERAGE"

            else:
                e_score = f"{np.minimum(cov / mean_cov, 1).mean():.4f}"
                status = "OK"


except Exception:
    status = "FAILED"


with open(args.output, "w") as out:
    out.write(
        "sample\te_score\tmean_depth\te_score_status\n"
    )

    out.write(
        f"{args.sample}\t{e_score}\t{mean_depth}\t{status}\n"
    )