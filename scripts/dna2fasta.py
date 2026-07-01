from Bio import SeqIO
from pathlib import Path
import sys


def dna2fasta(dna_fpath, fasta_fpath) -> None:
    for i, record in enumerate(SeqIO.parse(dna_fpath, "snapgene")):
        record.id = f'{Path(dna_fpath).name}_{i}'
        record.name = f'{Path(dna_fpath).name}_{i}'
        record.description = ''
        SeqIO.write(record, fasta_fpath, "fasta")
    print(f'Converted: {dna_fpath} -> {fasta_fpath}')


if __name__ == '__main__':
    dna_fpath = sys.argv[1] 
    fasta_fpath = sys.argv[2] if len(sys.argv) > 2 else f"{dna_fpath.rstrip('.dna')}.fa"
    dna2fasta(dna_fpath, fasta_fpath)
