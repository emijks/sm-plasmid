#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from samples import map_samples  # noqa: E402


def get_samples(table, seq_dir, ref_dir, seq_format):
    # те же правила сопоставления, что и в Snakefile (scripts/samples.py)
    return [
        {
            "sample": rec["name"],
            "order": rec["order"] or "-",
            "out_key": rec["out_key"],
            "reference": "OK" if rec["ref_ok"] else "MISSING",
            "fastq": "OK" if rec["fastq_ok"] else "MISSING",
        }
        for rec in map_samples(table, seq_dir, ref_dir, seq_format)
    ]



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



def make_status_table(table, seq_dir, ref_dir, seq_format, out_dir, output):
    """Строит статус-таблицу и записывает её в output (TSV). Возвращает DataFrame."""

    samples = get_samples(
        table,
        seq_dir,
        ref_dir,
        seq_format
    )


    out_dir = Path(out_dir)

    rows = []


    for sample in samples:

        name = sample["sample"]
        key = sample["out_key"]  # <OrderID>/<имя образца>, как в Snakefile


        row = {

            "sample": name,

            "order":
                sample["order"],

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
                f"{key}.vcf"
            )


            row.update(
                {

                    "assembly":
                        check_fasta(
                            out_dir /
                            f"{key}.assembly.fasta"
                        ),


                    "bam":
                        check_bam(
                            out_dir /
                            f"{key}.bam"
                        ),


                    "bam_index":
                        check_file(
                            out_dir /
                            f"{key}.bam.bai"
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
                            f"{key}.consensus.fasta"
                        ),


                    "report":
                        check_file(
                            out_dir /
                            f"{key}.report.html"
                        ),
                }
            )


        rows.append(row)


    df = pd.DataFrame(rows)


    df.to_csv(
        output,
        sep="\t",
        index=False
    )

    return df


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


    make_status_table(
        args.table,
        args.seq_dir,
        args.ref_dir,
        args.seq_format,
        args.out_dir,
        args.output
    )


if __name__ == "__main__":
    main()