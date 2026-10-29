import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(workflow.basedir) / 'scripts'))
from samples import map_samples  # noqa: E402  (общий модуль со scripts/make_status_table.py)


# # ============== utils ============================


seq_dir = str(Path(config['seq_dir']).resolve().absolute())   # только чтение
ref_dir = str(Path(config['ref_dir']).resolve().absolute())   # только чтение
tmp_dir = str(Path(config['tmp_dir']).resolve().absolute())   # каталог прогона (plasmid.py)
out_dir = str(Path(config['out_dir']).resolve().absolute())   # только результаты

# Раскладка tmp_dir (каталога прогона <tmp-dir>/sm-plasmid-<дата-время>-<суффикс>/):
#   refs/{ref_prefix}.fa*            — копия/конвертация референса и индексы (temp)
#   samples/{sample}/...         — промежуточные файлы образца (temp)
#   samples/{sample}/<rule>.work — рабочий каталог правила, удаляется в конце правила
# Wildcard {sample} — путь результатов относительно out_dir: <OrderID>/<имя образца>
# (см. scripts/samples.py), например SEQ-746/20117_SEQ-746-41_100003225.
ref_tmp = tmp_dir + '/refs'
smp_tmp = tmp_dir + '/samples'


def map_input() -> pd.DataFrame:
    """Образцы для обработки; печатает таблицу запуска, затем список пропущенных строк."""
    input_map_data, problems = [], []
    for rec in map_samples(config['table'], seq_dir, ref_dir, config['seq_format']):
        # Лишние и недостающие символы, кириллица в латинице - база эксель данных
        if not rec['ref_ok']:
            problems.append(f'{rec["name"]}: reference "{rec["ref_fname"]}" not found in "{ref_dir}"')
            continue
        if not rec['fastq_ok']:
            problems.append(f'{rec["name"]}: FASTQ skipped ({rec["fastq_error"]}) in "{seq_dir}"')
            continue
        if rec['order'] is None:
            problems.append(f'{rec["name"]}: warning: order not recognized in sample name '
                            f'(expected <RunID>_SEQ-<N>-<item>_<SampleID>), results go to "{out_dir}/{rec["out_key"]}.*"')
        # ruleorder dna2fa > fa2fa: при наличии обоих файлов используется .dna
        if Path(rec['ref_fname']).suffix == '.fa' and Path(f'{ref_dir}/{rec["ref_prefix"]}.dna').exists():
            problems.append(f'{rec["name"]}: warning: "{rec["ref_fname"]}" requested, but "{rec["ref_prefix"]}.dna" '
                            f'exists in "{ref_dir}" and will be used instead (ruleorder: dna2fa > fa2fa)')
        input_map_data.append({**rec, 'sample': rec['out_key']})
    input_map = pd.DataFrame(input_map_data, columns=['sample', 'name', 'order', 'ref_fname', 'ref_prefix', 'r1', 'r2'])

    launch = pd.DataFrame({
        'sample': input_map['name'],
        'order': input_map['order'].fillna('-'),
        'reference': input_map['ref_fname'],
        'files': [', '.join(Path(f).name for pair in zip(r1, r2) for f in pair)
                  for r1, r2 in zip(input_map['r1'], input_map['r2'])],
    })
    print(launch.to_string(index=False) if not launch.empty else 'No samples to process', flush=True)
    if problems:
        print(f'\nSkipped rows and warnings ({len(problems)}):', flush=True)
        print('\n'.join(f'  {p}' for p in problems), flush=True)
    return input_map

input_map = map_input()
if input_map.empty:
    sys.exit('Nothing to process')


def get_sample(wildcards) -> pd.Series:
    return input_map[input_map['sample'] == wildcards.sample].iloc[0]


def get_ref(suffix: str, wildcards) -> str:
    return f'{ref_tmp}/{get_sample(wildcards)["ref_prefix"]}{suffix}'


# Логи инструментов: вывод (stdout+stderr) каждого правила всегда пишется в его
# лог-файл; с --verbose-tools он дублируется в stdout. Сообщения Snakemake и
# обвязки (plasmid.py, таблица запуска) идут в stdout всегда.
#   <out-dir>/logs/<OrderID>/<имя образца>/<rule>.log
#   <out-dir>/logs/refs/<ref_prefix>/<rule>.log
#   <out-dir>/logs/summary_table.log
log_dir = out_dir + '/logs'
if config.get('quiet_tools', True):
    LOG = 'exec > "{log}" 2>&1\n'
else:
    LOG = 'exec > >(tee "{log}") 2>&1\n'


def sample_log(rule_name: str) -> str:
    return log_dir + '/{sample}/' + rule_name + '.log'


def ref_log(rule_name: str) -> str:
    return log_dir + '/refs/{ref_prefix}/' + rule_name + '.log'


def work_dir(rule_name: str) -> str:
    """Рабочий каталог правила для образца (в shell: {params.work})."""
    return smp_tmp + '/{sample}/' + rule_name + '.work'


# # ============== rules ============================


ruleorder: dna2fa > fa2fa


rule all:
    input:
        expand(out_dir + '/{sample}.bam', sample=input_map['sample']),
        expand(out_dir + '/{sample}.vcf', sample=input_map['sample']),
        expand(out_dir + '/{sample}.consensus.fasta', sample=input_map['sample']),
        expand(out_dir + '/{sample}.assembly.fasta', sample=input_map['sample']),
        expand(out_dir + '/{sample}.report.html', sample=input_map['sample']),
        out_dir + '/summary.tsv'


rule dna2fa:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_dna = ref_dir + '/{ref_prefix}.dna',
    output:
        ref_fa = temp(ref_tmp + '/{ref_prefix}.fa')
    log:
        ref_log('dna2fa')
    shell:
        LOG + '''
        python3 scripts/dna2fasta.py "{input.ref_dna}" "{output.ref_fa}"
        '''


rule fa2fa:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = ref_dir + '/{ref_prefix}.fa',
    output:
        ref_fa = temp(ref_tmp + '/{ref_prefix}.fa')
    log:
        ref_log('fa2fa')
    shell:
        LOG + '''
        cp "{input.ref_fa}" "{output.ref_fa}"
        '''


rule bwa_index:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = ref_tmp + '/{ref_prefix}.fa',
    output:
        ref_0123 = temp(ref_tmp + '/{ref_prefix}.fa.0123'),
        ref_amb = temp(ref_tmp + '/{ref_prefix}.fa.amb'),
        ref_ann = temp(ref_tmp + '/{ref_prefix}.fa.ann'),
        ref_bwt = temp(ref_tmp + '/{ref_prefix}.fa.bwt.2bit.64'),
        ref_pac = temp(ref_tmp + '/{ref_prefix}.fa.pac'),
    log:
        ref_log('bwa_index')
    shell:
        LOG + '''
        bwa-mem2 index "{input.ref_fa}"
        '''


rule faidx:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        ref_fa = ref_tmp + '/{ref_prefix}.fa',
    output:
        ref_fai = temp(ref_tmp + '/{ref_prefix}.fa.fai'),
    log:
        ref_log('faidx')
    shell:
        LOG + '''
        samtools faidx "{input.ref_fa}"
        '''


rule concat_lanes:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        # точные списки лейнов из scripts/samples.py, R1[i] и R2[i] — один лейн
        r1 = lambda wildcards: get_sample(wildcards)['r1'],
        r2 = lambda wildcards: get_sample(wildcards)['r2'],
    output:
        r1 = temp(smp_tmp + '/{sample}/R1.concat.fq.gz'),
        r2 = temp(smp_tmp + '/{sample}/R2.concat.fq.gz'),
    log:
        sample_log('concat_lanes')
    shell:
        LOG + '''
        cat {input.r1:q} > "{output.r1}"
        cat {input.r2:q} > "{output.r2}"
        '''


rule fastqc:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = rules.concat_lanes.output.r1,
        r2 = rules.concat_lanes.output.r2,
    output:
        zip1 = temp(smp_tmp + '/{sample}/R1.fastqc.zip'),
        zip2 = temp(smp_tmp + '/{sample}/R2.fastqc.zip'),
    params:
        work = work_dir('fastqc'),
    log:
        sample_log('fastqc')
    shell:
        LOG + '''
        rm -rf "{params.work}" && mkdir -p "{params.work}"
        fastqc "{input.r1}" "{input.r2}" -o "{params.work}" -q -t {threads}
        mv "{params.work}/R1.concat_fastqc.zip" "{output.zip1}"
        mv "{params.work}/R2.concat_fastqc.zip" "{output.zip2}"
        rm -rf "{params.work}"
        '''


rule trim:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = rules.concat_lanes.output.r1,
        r2 = rules.concat_lanes.output.r2,
    output:
        r1p = temp(smp_tmp + '/{sample}/R1.paired.fq.gz'),
        r1u = temp(smp_tmp + '/{sample}/R1.unpaired.fq.gz'),
        r2p = temp(smp_tmp + '/{sample}/R2.paired.fq.gz'),
        r2u = temp(smp_tmp + '/{sample}/R2.unpaired.fq.gz'),
    params:
        adapters = 'references/TruSeq-illumina-PE.fa:2:30:10',
        minlen = config.get('trim_min_len', 50),
    log:
        sample_log('trim')
    shell:
        LOG + '''
        trimmomatic PE -threads {threads} \
            "{input.r1}" \
            "{input.r2}" \
            "{output.r1p}" \
            "{output.r1u}" \
            "{output.r2p}" \
            "{output.r2u}" \
            ILLUMINACLIP:{params.adapters} \
            LEADING:10 TRAILING:10 SLIDINGWINDOW:4:20 MINLEN:{params.minlen}
        '''


rule align:
    threads: 4
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        r1 = rules.trim.output.r1p,
        r2 = rules.trim.output.r2p,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
        ref_idx = lambda wildcards: [get_ref(s, wildcards) for s in
                                     ('.fa.0123', '.fa.amb', '.fa.ann', '.fa.bwt.2bit.64', '.fa.pac')],
    output:
        bam = out_dir + '/{sample}.bam',
        bai = out_dir + '/{sample}.bam.bai',
    params:
        work = work_dir('align'),
    log:
        sample_log('align')
    shell:
        LOG + '''
        rm -rf "{params.work}" && mkdir -p "{params.work}"
        bwa-mem2 mem -t {threads} "{input.ref_fa}" -R "@RG\tID:sample\tPL:Illumina\tSM:pool" "{input.r1}" "{input.r2}" \
        | samtools view -bh -@ {threads} - \
        | samtools sort -n -@ {threads} -T "{params.work}/sort_name" - \
        | samtools fixmate -m -@ {threads} - - \
        | samtools sort -@ {threads} -T "{params.work}/sort_coord" - \
        | samtools markdup -@ {threads} -T "{params.work}/markdup" -s - "{output.bam}"
        samtools index -@ {threads} "{output.bam}"
        rm -rf "{params.work}"
        '''


rule call:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        bam = rules.align.output.bam,
        ref_fa = lambda wildcards: get_ref(".fa", wildcards),
        ref_fai = lambda wildcards: get_ref('.fa.fai', wildcards),
    output:
        vcf = out_dir + "/{sample}.vcf",
    log:
        sample_log('call')
    shell:
        LOG + '''
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
        ref_fai = lambda wildcards: get_ref('.fa.fai', wildcards)
    output:
        fa = out_dir + '/{sample}.consensus.fasta'
    params:
        work = work_dir('consensus'),
    log:
        sample_log('consensus')
    shell:
        LOG + '''
        scripts/consensus.sh "{params.work}" \
            "{input.bam}" "{input.vcf}" "{input.ref_fa}" "{output.fa}"
        '''


rule assembly:
    threads: 4
    singularity: 'docker://kazakovaem/spades:1.0.0'
    input:
        r1 = rules.trim.output.r1p,
        r2 = rules.trim.output.r2p,
    output:
        fasta = out_dir + '/{sample}.assembly.fasta',
    params:
        work = work_dir('assembly'),
    log:
        sample_log('assembly')
    shell:
        LOG + """
        rm -rf "{params.work}"
        /usr/bin/python3.8 /SPAdes-4.2.0/bin/spades.py \
            -1 "{input.r1}" \
            -2 "{input.r2}" \
            -t {threads} \
            -o "{params.work}"

        cp "{params.work}/contigs.fasta" "{output.fasta}"

        rm -rf "{params.work}"
        """


rule igv:
    threads: 1
    singularity: 'docker://albertea/igv:2.16.2'
    input:
        bam = rules.align.output.bam,
        bai = rules.align.output.bai,
        vcf = rules.call.output.vcf,
        ref_fa = lambda wildcards: get_ref('.fa', wildcards),
        ref_fai = lambda wildcards: get_ref('.fa.fai', wildcards)
    output:
        png = temp(smp_tmp + '/{sample}/igv.png')
    params:
        work = work_dir('igv'),
    log:
        sample_log('igv')
    shell:
        LOG + '''
        scripts/igv.sh "{params.work}" \
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
        ref_fai = lambda wildcards: get_ref('.fa.fai', wildcards)
    output:
        html = out_dir + '/{sample}.report.html'
    params:
        work = work_dir('report'),
    log:
        sample_log('report')
    shell:
        LOG + '''
        rm -rf "{params.work}" && mkdir -p "{params.work}"
        name="$(basename "{wildcards.sample}")"
        samtools flagstat "{input.bam}" > "{params.work}/${{name}}_flagstat.txt"
        samtools depth -a "{input.bam}" > "{params.work}/${{name}}_depth.txt"
        cp "{input.png}" "{params.work}/${{name}}.igv.png"
        Rscript scripts/report_generation.r \
            "$name" \
            "{input.ref_fa}" \
            "{params.work}" \
            "{input.fastqc_r1}" \
            "{input.fastqc_r2}"
        cp "{params.work}/${{name}}_report.html" "{output.html}"
        rm -rf "{params.work}"
        '''


rule calculate_evenness:
    threads: 1
    singularity: 'docker://kazakovaem/plasmid:1.0.0'
    input:
        bam = rules.align.output.bam,
        bai = rules.align.output.bai,
    output:
        tsv = temp(smp_tmp + "/{sample}/e_score.tsv")
    log:
        sample_log('calculate_evenness')
    shell:
        LOG + '''
        python3 scripts/calculate_evenness.py \
            --bam "{input.bam}" \
            --sample "$(basename "{wildcards.sample}")" \
            --output "{output.tsv}"
        '''


rule summary_table:
    input:
        expand(
            smp_tmp + "/{sample}/e_score.tsv",
            sample=input_map["sample"]
        )
    output:
        out_dir + "/summary.tsv"
    log:
        log_dir + '/summary_table.log'
    shell:
        LOG + '''
        files=({input:q})
        head -n1 "${{files[0]}}" > "{output}"
        for f in "${{files[@]}}"; do
            tail -n +2 "$f" >> "{output}"
        done
        '''
