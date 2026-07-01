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
        # # То одна структура данных
        # if config['seq_format'] == 1:
        #     seq_prefix = f"{row['ID']}/{row['Folder']}/{row['Sample']}"
        # # То другая
        # else:
        #     seq_id = str(Path(seq_dir).name).rstrip('-fastq')
        #     seq_prefix = f"{row['ID']}/{seq_id}_{row['Folder']}_{row['Sample']}"
        #     # seq_prefix = f"{row['ID']}/{seq_id}_{row['Sample']}"
        # То одна структура данных, то другая ><
        if config['seq_format'] == 1:
            seq_id = str(Path(seq_dir).name).rstrip('-fastq')
            seq_prefix = f"{row['ID']}/{row['Sample']}"
        elif config['seq_format'] == 2:
            seq_id = str(Path(seq_dir).name).rstrip('-fastq')
            seq_prefix = f"{row['ID']}/{seq_id}_{row['Folder']}_{row['Sample']}"
        elif config['seq_format'] == 3:
            seq_id = str(Path(seq_dir).name).rstrip('-fastq')
            seq_prefix = f"{row['ID']}/{seq_id}_{row['Sample']}"
        elif config['seq_format'] == 4:
            seq_prefix = f"{row['ID']}/{row['Folder']}/{row['Sample']}"
        elif config['seq_format'] == 5:
            seq_prefix = f"{row['Sample']}"
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
    ref = f'{ref_dir}/{ref_prefix}{suffix}'
    return ref


# # ============== rules ============================


rule all:
    input:
        expand(out_dir + '/{seq_prefix}.report.html', seq_prefix=input_map['seq_prefix']),
        expand(out_dir + '/{seq_prefix}.consensus.fasta', seq_prefix=input_map['seq_prefix'])


rule dna2fa:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_dna = ref_dir + '/{ref_prefix}.dna',
    output:
        ref_fa = ref_dir + '/{ref_prefix}.fa'
    shell:
        '''
        python3 scripts/dna2fasta.py "{input.ref_dna}" "{output.ref_fa}"
        '''


rule bwa_index:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = ref_dir + '/{ref_prefix}.fa',
    output:
        ref_0123 = ref_dir + '/{ref_prefix}.fa.0123',
        ref_amb = ref_dir + '/{ref_prefix}.fa.amb',
        ref_ann = ref_dir + '/{ref_prefix}.fa.ann',
        ref_bwt = ref_dir + '/{ref_prefix}.fa.bwt.2bit.64',
        ref_pac = ref_dir + '/{ref_prefix}.fa.pac',
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
        html1 = temp(out_dir + '/{seq_prefix}.R1.fastqc.html'),
        html2 = temp(out_dir + '/{seq_prefix}.R2.fastqc.html'),
    shell:
        '''
        tmp_seqdir=$(dirname "{tmp_dir}/{wildcards.seq_prefix}")
        mkdir -p "$tmp_seqdir"
        fastqc "{input.r1}" "{input.r2}" -o "$tmp_seqdir" -q -t {threads}
        mv "{tmp_dir}/{wildcards.seq_prefix}"*_R1_*fastqc.html "{output.html1}"
        mv "{tmp_dir}/{wildcards.seq_prefix}"*_R2_*fastqc.html "{output.html2}"
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
        bwa-mem2 mem -p -t {threads} "{input.ref_fa}" -R "@RG\tID:sample\tPL:Illumina\tSM:pool" "{input.r1}" "{input.r2}" \
        | samtools view -bh -@ {threads} - \
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
        bai = rules.align.output.bai,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
    output:
        vcf = out_dir + '/{seq_prefix}.vcf'
    shell:
        '''
        bcftools mpileup -f "{input.ref_fa}" "{input.bam}" \
        | bcftools view -u \
        | bcftools call --variants-only --ploidy 1 -m -o "{output.vcf}"
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
        bedtools genomecov -bga -ibam "{input.bam}" | awk '$4==0' > "{tmp_dir}/{wildcards.seq_prefix}.zerocov.bed"
        bcftools view "{input.vcf}" -Oz -o "{tmp_dir}/{wildcards.seq_prefix}.vcf.gz"
        bcftools index "{tmp_dir}/{wildcards.seq_prefix}.vcf.gz"
        bcftools consensus -f "{input.ref_fa}" -m "{tmp_dir}/{wildcards.seq_prefix}.zerocov.bed" "{tmp_dir}/{wildcards.seq_prefix}.vcf.gz" > "{output.fa}"
        '''


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
        scripts/igv.sh "{tmp_dir}" "{wildcards.seq_prefix}" "{input.bam}" "{input.vcf}" "{input.ref_fa}" "{output.png}"
        '''


rule report:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        html1 = rules.fastqc.output.html1,
        html2 = rules.fastqc.output.html2,
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
        Rscript scripts/report_generation.r "$(basename "{wildcards.seq_prefix}")" "{input.ref_fa}" "$tmp_seqdir"
        pandoc -s "{tmp_dir}/{wildcards.seq_prefix}_report.html" "{input.html1}" "{input.html2}" -o "{output.html}"
        '''
