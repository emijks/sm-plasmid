import pandas as pd
import glob
import sys
from pathlib import Path


# # ============== utils ============================


seq_dir = str(Path(config['seq_dir']).resolve().absolute())
ref_dir = str(Path(config['ref_dir']).resolve().absolute())
tmp_dir = str(Path(config['tmp_dir']).resolve().absolute())
out_dir = str(Path(config['out_dir']).resolve().absolute())


def map_input() -> pd.DataFrame:
    input_map_data = []
    raw_data = pd.read_excel(config['table'])
    for i, row in raw_data.iterrows():
        ref_fname = str(row['Map_name']).strip()
        ref_prefix = str(Path(ref_fname).stem)
        # То одна структура данных
        if config['seq_format'] == 1:
            seq_prefix = f"{row['ID']}/{row['Folder']}/{row['Sample']}"
        # То другая
        else:
            seq_id = str(Path(seq_dir).name).rstrip('-fastq')
            seq_prefix = f"{row['Sample']}"
            #seq_prefix = f"{row['ID']}/{seq_id}_{row['Folder']}_{row['Sample']}"
            #seq_prefix = f"{row['ID']}/{seq_id}_{row['Sample']}"
        # Лишние и недостающие символы, кириллица в латинице - база эксель данных
        if not Path(f'{ref_dir}/{ref_fname}').exists():
            print(f'Reference "{ref_fname}" not found in "{ref_dir}"')
            continue
        if not len(glob.glob(f'{seq_dir}/{seq_prefix}*')) == 2:
            print(f'Sequence "{seq_prefix}" not found as pair in "{seq_dir}"')
            continue
        input_map_data.append({
            'ref_prefix': ref_prefix,
            'seq_prefix': seq_prefix
        })
    input_map = pd.DataFrame(input_map_data)
    print(input_map)
    return input_map

input_map = map_input()
if input_map.empty:
    sys.exit('Nothing to process')


def get_seq(suffix_pattern: str, wildcards: snakemake.io.Wildcards) -> str:
    seq = glob.glob(f'{seq_dir}/{wildcards.seq_prefix}{suffix_pattern}')
    return seq[0]


def get_ref(suffix: str, wildcards: snakemake.io.Wildcards) -> str:
    ref_prefix = next(iter(input_map[input_map['seq_prefix'] == wildcards.seq_prefix]['ref_prefix']))
    ref = f'{tmp_dir}/{ref_prefix}{suffix}'
    return ref


# # ============== rules ============================


rule all:
    input:
        expand(out_dir + '/{seq_prefix}.report.html', seq_prefix=input_map['seq_prefix']),
        expand(out_dir + '/{seq_prefix}.consensus.fasta', seq_prefix=input_map['seq_prefix']),
        expand(out_dir + '/{seq_prefix}.assembly.fasta', seq_prefix=input_map['seq_prefix']),
        expand(out_dir + '/{seq_prefix}.vcf', seq_prefix=input_map['seq_prefix']),
        out_dir + '/summary.tsv'


rule dna2fa:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_dna = ref_dir + '/{ref_prefix}.dna',
    output:
        ref_fa = tmp_dir + '/{ref_prefix}.fa'
    shell:
        '''
        python3 scripts/dna2fasta.py "{input.ref_dna}" "{output.ref_fa}"
        '''

rule fa2fa:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = ref_dir + '/{ref_prefix}.fa',
    output:
        ref_fa = tmp_dir + '/{ref_prefix}.fa'
    shell:
        '''
        cp "{input.ref_fa}" "{output.ref_fa}"
        '''


rule bwa_index:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = tmp_dir + '/{ref_prefix}.fa',
    output:
        ref_0123 = tmp_dir + '/{ref_prefix}.fa.0123',
        ref_amb = tmp_dir + '/{ref_prefix}.fa.amb',
        ref_ann = tmp_dir + '/{ref_prefix}.fa.ann',
        ref_bwt = tmp_dir + '/{ref_prefix}.fa.bwt.2bit.64',
        ref_pac = tmp_dir + '/{ref_prefix}.fa.pac',
    shell:
        '''
        bwa-mem2 index "{input.ref_fa}"
        '''


rule fastqc:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = lambda wildcards: get_seq('*R1*', wildcards),
        r2 = lambda wildcards: get_seq('*R2*', wildcards),
    output:
        zip1  = temp(out_dir + '/{seq_prefix}.R1.fastqc.zip'),
        zip2  = temp(out_dir + '/{seq_prefix}.R2.fastqc.zip'),
    shell:
        '''
        tmp_seqdir=$(dirname "{tmp_dir}/{wildcards.seq_prefix}")
        mkdir -p "$tmp_seqdir"
        fastqc "{input.r1}" "{input.r2}" -o "$tmp_seqdir" -q -t {threads}
        mv "{tmp_dir}/{wildcards.seq_prefix}"*_R1_*fastqc.zip "{output.zip1}"
        mv "{tmp_dir}/{wildcards.seq_prefix}"*_R2_*fastqc.zip "{output.zip2}"
        rm -f "{tmp_dir}/{wildcards.seq_prefix}"*_fastqc.html
        '''


rule trim:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = lambda wildcards: get_seq('*R1*', wildcards),
        r2 = lambda wildcards: get_seq('*R2*', wildcards),
    output:
        r1p = temp(out_dir + '/{seq_prefix}_R1_paired.fq.gz'),
        r1u = temp(out_dir + '/{seq_prefix}_R1_unpaired.fq.gz'),
        r2p = temp(out_dir + '/{seq_prefix}_R2_paired.fq.gz'),
        r2u = temp(out_dir + '/{seq_prefix}_R2_unpaired.fq.gz'),
    params:
        adapters = 'references/TruSeq-illumina-PE.fa:2:30:10'
    shell:
        '''
        trimmomatic PE -threads {threads} \
            "{input.r1}" \
            "{input.r2}" \
            "{output.r1p}" \
            "{output.r1u}" \
            "{output.r2p}" \
            "{output.r2u}" \
            ILLUMINACLIP:{params.adapters} \
            LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:50
        '''


rule align:
    threads: 4
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = rules.trim.output.r1p,
        r2 = rules.trim.output.r2p,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
        ref_bwt = lambda wildcards: get_ref('.fa.bwt.2bit.64', wildcards),
    output:
        bam = out_dir + '/{seq_prefix}.bam',
        bai = out_dir + '/{seq_prefix}.bam.bai',

    shell:
        '''
        bwa-mem2 mem -t {threads} "{input.ref_fa}" -R "@RG\tID:sample\tPL:Illumina\tSM:pool" "{input.r1}" "{input.r2}" \
        | samtools view -bh -@ {threads} - \
        | samtools sort -n -@ {threads} -T {tmp_dir} - \
        | samtools fixmate -m -@ {threads} - - \
        | samtools sort -@ {threads} -T {tmp_dir} - \
        | samtools markdup -@ {threads} -T {tmp_dir} -s - "{output.bam}"
        samtools index -@ {threads} "{output.bam}"
        '''


rule call:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        bam = rules.align.output.bam,
        ref_fa = lambda wildcards: get_ref(".fa", wildcards),
    output:
        vcf = out_dir + "/{seq_prefix}.vcf",
    shell:
        '''
        samtools mpileup \
            -aa \
            -A \
            -B \
            -d 100000 \
            -f "{input.ref_fa}" \
            "{input.bam}" \
        | python3 scripts/brute-force_caller.py \
            --min-depth 10 \
            --min-af 0.5 \
            > "{output.vcf}"
       '''

rule consensus:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        vcf = rules.call.output.vcf,
        bam = rules.align.output.bam,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
    output:
        fa = out_dir + '/{seq_prefix}.consensus.fasta'
    shell:
        '''
        scripts/consensus.sh "{tmp_dir}/{wildcards.seq_prefix}" \
            "{input.bam}" "{input.vcf}" "{input.ref_fa}" "{output.fa}"
        '''


rule assembly:
    threads: 4
    singularity: 'docker://kazakovaem/spades:1.0.0'
    input:
        r1 = rules.trim.output.r1p,
        r2 = rules.trim.output.r2p,
    output:
        fasta = out_dir + '/{seq_prefix}.assembly.fasta',
    shell:
        """
        /usr/bin/python3.8 /SPAdes-4.2.0/bin/spades.py \
            -1 "{input.r1}" \
            -2 "{input.r2}" \
            -t {threads} \
            -o "{tmp_dir}/{wildcards.seq_prefix}.spades"

        cp "{tmp_dir}/{wildcards.seq_prefix}.spades/contigs.fasta" "{output.fasta}"

        rm -rf "{tmp_dir}/{wildcards.seq_prefix}.spades"
        """


rule igv:
    threads: 1
    singularity: 'docker://albertea/igv:2.16.2'
    input:
        bam = rules.align.output.bam,
        vcf = rules.call.output.vcf,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
    output:
        txt = temp(out_dir + '/{seq_prefix}.igv.png.txt'),
        png = temp(out_dir + '/{seq_prefix}.igv.png')
    shell:
        '''
        scripts/igv.sh "{tmp_dir}" "{wildcards.seq_prefix}" \
            "{input.bam}" "{input.vcf}" "{input.ref_fa}" "{output.png}"
        '''


rule report:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        fastqc_r1 = rules.fastqc.output.zip1,
        fastqc_r2 = rules.fastqc.output.zip2,
        bam = rules.align.output.bam,
        bai = rules.align.output.bai,
        png = rules.igv.output.png,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
    output:
        html = out_dir + '/{seq_prefix}.report.html'
    shell:
        '''
        tmp_seqdir=$(dirname "{tmp_dir}/{wildcards.seq_prefix}")
        mkdir -p "$tmp_seqdir"
        samtools flagstat "{input.bam}" > "{tmp_dir}/{wildcards.seq_prefix}_flagstat.txt"
        samtools depth -a "{input.bam}" > "{tmp_dir}/{wildcards.seq_prefix}_depth.txt"
        Rscript scripts/report_generation.r \
            "$(basename "{wildcards.seq_prefix}")" \
            "{input.ref_fa}" \
            "$tmp_seqdir" \
            "{input.fastqc_r1}" \
            "{input.fastqc_r2}"
        cp "{tmp_dir}/{wildcards.seq_prefix}_report.html" "{output.html}"
        
        rm -f \
            "{tmp_dir}/{wildcards.seq_prefix}_flagstat.txt" \
            "{tmp_dir}/{wildcards.seq_prefix}_depth.txt"
        '''
#rm "{tmp_dir}/{wildcards.seq_prefix}"*

rule calculate_evenness:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        bam = out_dir + "/{seq_prefix}.bam",
        bai = out_dir + "/{seq_prefix}.bam.bai"
    output:
        tsv = tmp_dir + "/{seq_prefix}.e_score.tsv"
    shell:
        '''
        python3 scripts/calculate_evenness.py \
            --bam {input.bam} \
            --sample {wildcards.seq_prefix} \
            --output {output.tsv}
        '''


rule summary_table:
    input:
        expand(
            tmp_dir + "/{seq_prefix}.e_score.tsv",
            seq_prefix=input_map["seq_prefix"]
        )
    output:
        out_dir + "/summary.tsv"
    shell:
        '''
        head -n1 {input[0]} > {output}
        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        '''
