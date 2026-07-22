#!/usr/bin/env python3

import argparse
import glob
from pathlib import Path

import pandas as pd


def get_samples(table, seq_dir, ref_dir, seq_format):

    samples = []

    raw_data = pd.read_excel(table)

    for _, row in raw_data.iterrows():

        ref_fname = str(row["Map_name"]).strip()

        if seq_format == 1:
            seq_prefix = (
                f"{row['ID']}/"
                f"{row['Folder']}/"
                f"{row['Sample']}"
            )
        else:
            seq_prefix = str(row["Sample"])


        ref_ok = (
            Path(ref_dir) / ref_fname
        ).exists()


        fastq_ok = (
            len(
                glob.glob(
                    f"{seq_dir}/{seq_prefix}*"
                )
            )
            == 2
        )


        samples.append(
            {
                "sample": seq_prefix,
                "reference":
                    "OK" if ref_ok else "MISSING",
                "fastq":
                    "OK" if fastq_ok else "MISSING",
            }
        )

    return samples



def check_file(path):

    path = Path(path)

    if not path.exists():
        return "MISSING"

    if path.stat().st_size == 0:
        return "EMPTY"

    return "OK"



def check_bam(path):

    path = Path(path)

    if not path.exists():
        return "MISSING"

    if path.stat().st_size == 0:
        return "EMPTY"

    return "OK"



def check_fasta(path):

    path = Path(path)

    if not path.exists():
        return "MISSING"

    if path.stat().st_size == 0:
        return "EMPTY"

    with open(path) as f:

        for line in f:

            if line.startswith(">"):
                return "OK"

    return "INVALID"



def check_vcf(path):

    """
    Check brute-force caller VCF.

    Expected INFO field:
        TYPE=SNP
        TYPE=INS
        TYPE=DEL
    """

    path = Path(path)

    if not path.exists():

        return {
            "status": "MISSING",
            "snps": 0,
            "indels": 0
        }


    if path.stat().st_size == 0:

        return {
            "status": "EMPTY",
            "snps": 0,
            "indels": 0
        }


    snps = 0
    indels = 0


    with open(path) as f:

        for line in f:

            if line.startswith("#"):
                continue


            if "TYPE=SNP" in line:

                snps += 1


            elif (
                "TYPE=INS" in line
                or "TYPE=DEL" in line
            ):

                indels += 1



    if snps + indels == 0:

        status = "NO_VARIANTS"

    else:

        status = "OK"


    return {
        "status": status,
        "snps": snps,
        "indels": indels
    }



def skipped():

    return {

        "assembly": "SKIPPED",

        "bam": "SKIPPED",

        "bam_index": "SKIPPED",

        "variant_calling": "SKIPPED",

        "snp_count": 0,

        "indel_count": 0,

        "consensus": "SKIPPED",

        "report": "SKIPPED",
    }



def main():

    parser = argparse.ArgumentParser()


    parser.add_argument(
        "--table",
        required=True
    )

    parser.add_argument(
        "--seq-dir",
        required=True
    )

    parser.add_argument(
        "--ref-dir",
        required=True
    )

    parser.add_argument(
        "--seq-format",
        type=int,
        required=True
    )

    parser.add_argument(
        "--out-dir",
        required=True
    )

    parser.add_argument(
        "--output",
        required=True
    )


    args = parser.parse_args()


    samples = get_samples(
        args.table,
        args.seq_dir,
        args.ref_dir,
        args.seq_format
    )


    out_dir = Path(args.out_dir)

    rows = []


    for sample in samples:

        name = sample["sample"]


        row = {

            "sample": name,

            "reference":
                sample["reference"],

            "fastq":
                sample["fastq"],

        }


        if (
            sample["reference"] != "OK"
            or sample["fastq"] != "OK"
        ):

            row.update(
                skipped()
            )


        else:

            variant_result = check_vcf(
                out_dir /
                f"{name}.vcf"
            )


            row.update(
                {

                    "assembly":
                        check_fasta(
                            out_dir /
                            f"{name}.assembly.fasta"
                        ),


                    "bam":
                        check_bam(
                            out_dir /
                            f"{name}.bam"
                        ),


                    "bam_index":
                        check_file(
                            out_dir /
                            f"{name}.bam.bai"
                        ),


                    "variant_calling":
                        variant_result["status"],


                    "snp_count":
                        variant_result["snps"],


                    "indel_count":
                        variant_result["indels"],


                    "consensus":
                        check_fasta(
                            out_dir /
                            f"{name}.consensus.fasta"
                        ),


                    "report":
                        check_file(
                            out_dir /
                            f"{name}.report.html"
                        ),
                }
            )


        rows.append(row)


    df = pd.DataFrame(rows)


    df.to_csv(
        args.output,
        sep="\t",
        index=False
    )


if __name__ == "__main__":
    main()