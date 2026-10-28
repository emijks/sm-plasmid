#!/usr/bin/env python3

import sys
import argparse
import re
from collections import Counter


def parse_bases(ref, bases):
    """
    Parse samtools mpileup base column.

    Returns:
        SNP counts
        insertions
        deletions
    """

    counts = Counter()
    insertions = Counter()
    deletions = Counter()

    i = 0

    while i < len(bases):

        c = bases[i]

        # read start marker: ^ + mapping quality
        if c == "^":
            i += 2
            continue

        # read end marker
        if c == "$":
            i += 1
            continue

        # insertion/deletion
        if c in "+-":

            sign = c
            i += 1

            number = ""

            while i < len(bases) and bases[i].isdigit():
                number += bases[i]
                i += 1

            if number == "":
                continue

            length = int(number)

            seq = bases[i:i+length].upper()

            if sign == "+":
                insertions[seq] += 1
            else:
                deletions[seq] += 1

            i += length
            continue


        # matches reference
        if c == "." or c == ",":
            counts[ref.upper()] += 1

        # mismatch
        elif c.upper() in "ACGTN":
            counts[c.upper()] += 1


        i += 1


    return counts, insertions, deletions



def write_variant(out,
                  chrom,
                  pos,
                  ref,
                  alt,
                  dp,
                  ad,
                  af,
                  vtype):

    info = (
        f"DP={dp};"
        f"AD={ad};"
        f"AF={af:.4f};"
        f"TYPE={vtype}"
    )

    out.write(
        f"{chrom}\t"
        f"{pos}\t"
        f".\t"
        f"{ref}\t"
        f"{alt}\t"
        f".\t"
        f"PASS\t"
        f"{info}\n"
    )



def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--min-depth",
        type=int,
        default=10
    )

    parser.add_argument(
        "--min-af",
        type=float,
        default=0.5
    )

    args = parser.parse_args()


    out = sys.stdout


    out.write(
        "##fileformat=VCFv4.2\n"
    )

    out.write(
        "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
    )

    out.write(
        "##INFO=<ID=AD,Number=1,Type=Integer,Description=\"Alt reads\">\n"
    )

    out.write(
        "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele frequency\">\n"
    )

    out.write(
        "##INFO=<ID=TYPE,Number=1,Type=String,Description=\"Variant type\">\n"
    )

    out.write(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    )


    for line in sys.stdin:

        if line.startswith("#"):
            continue

        fields = line.rstrip().split("\t")

        if len(fields) < 6:
            continue


        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[2].upper()
        depth = int(fields[3])
        bases = fields[4]


        if depth < args.min_depth:
            continue


        counts, insertions, deletions = parse_bases(
            ref,
            bases
        )


        # SNPs

        for base, count in counts.items():

            if base == ref:
                continue

            af = count / depth

            if af >= args.min_af:

                write_variant(
                    out,
                    chrom,
                    pos,
                    ref,
                    base,
                    depth,
                    count,
                    af,
                    "SNP"
                )


        # Insertions

        for seq, count in insertions.items():

            af = count / depth

            if af >= args.min_af:

                write_variant(
                    out,
                    chrom,
                    pos,
                    ref,
                    ref + seq,
                    depth,
                    count,
                    af,
                    "INS"
                )


        # Deletions

        for seq, count in deletions.items():

            af = count / depth

            if af >= args.min_af:

                write_variant(
                    out,
                    chrom,
                    pos,
                    ref + seq,
                    ref,
                    depth,
                    count,
                    af,
                    "DEL"
                )


if __name__ == "__main__":
    main()