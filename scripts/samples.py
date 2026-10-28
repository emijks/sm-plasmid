"""Сопоставление строк xlsx-таблицы с референсами и FASTQ-файлами.

Общий модуль для Snakefile и scripts/make_status_table.py, чтобы правила
поиска файлов были одни и те же.

FASTQ ищутся по точному имени образца, а не по префиксу (раньше маска
"S1*R1*" захватывала и файлы S10). После имени образца допускаются только
стандартные для Illumina части имени:

    {Sample}[_S<n>][_L<lane>](_|.)R<1|2>[_<NNN>].(fastq|fq).gz

Примеры: S1_S3_L001_R1_001.fastq.gz, S1_L002_R2.fq.gz, S1_R1.fastq.gz, S1.R2.fq.gz,
20113_SEQ-821-6_100003544_L03_R1.fq.gz

Имя образца имеет вид <RunID>_<OrderID>-<item>_<SampleID>, например
20117_SEQ-746-41_100003225: RunID=20117, OrderID=SEQ-746, item=41,
SampleID=100003225. Результаты образца кладутся в <out-dir>/<OrderID>/.
"""

import re
from pathlib import Path

import pandas as pd

FASTQ_SUFFIX_RE = r'(?:_S\d+)?(?:_L\d+)?[._]R(?P<read>[12])(?:_\d+)?\.f(?:ast)?q\.gz$'
SAMPLE_NAME_RE = re.compile(r'^(?P<run>\d+)_(?P<order>SEQ-\d+)-(?P<item>\d+)_(?P<sample_id>.+)$')
NO_ORDER = 'NO_ORDER'  # подкаталог <out-dir> для имён, из которых не извлекается заказ


def parse_sample_name(name: str):
    """{'run', 'order', 'item', 'sample_id'} или None, если имя не по формату."""
    m = SAMPLE_NAME_RE.match(name)
    return m.groupdict() if m else None


def out_key(name: str, order) -> str:
    """Путь результатов образца относительно <out-dir> без расширения: <OrderID>/<имя образца>."""
    return f'{order or NO_ORDER}/{name}'


def seq_prefix_for(row, seq_format: int) -> str:
    if seq_format == 1:
        return f"{row['ID']}/{row['Folder']}/{str(row['Sample']).strip()}"
    return str(row['Sample']).strip()


def find_fastq_lanes(seq_dir: str, seq_prefix: str):
    """Возвращает (r1_files, r2_files, error).

    Списки отсортированы так, что r1_files[i] и r2_files[i] — один лейн.
    error — текст ошибки или None.
    """
    target = Path(seq_dir) / seq_prefix
    parent, name = target.parent, target.name
    if not parent.is_dir():
        return [], [], f'directory "{parent}" not found'

    rx = re.compile('^' + re.escape(name) + FASTQ_SUFFIX_RE)
    lanes = {'1': {}, '2': {}}
    for f in sorted(parent.iterdir()):
        m = rx.match(f.name)
        if not m or not f.is_file():
            continue
        # ключ лейна: имя файла без номера рида
        key = f.name[:m.start('read')] + f.name[m.end('read'):]
        lanes[m['read']][key] = str(f)

    r1, r2 = lanes['1'], lanes['2']
    if not r1 and not r2:
        return [], [], 'no FASTQ files found'
    if set(r1) != set(r2):
        unpaired = sorted(set(r1) ^ set(r2))
        return [], [], f'R1/R2 lanes do not match (unpaired: {unpaired})'
    keys = sorted(r1)
    return [r1[k] for k in keys], [r2[k] for k in keys], None


def map_samples(table: str, seq_dir: str, ref_dir: str, seq_format: int) -> list:
    """Одна запись на строку таблицы: пути к референсу и FASTQ плюс статусы."""
    records = []
    for _, row in pd.read_excel(table).iterrows():
        # пустые строки (в т.ч. ячейки с "" от формул) пропускаем
        if any(pd.isna(row[c]) or not str(row[c]).strip() for c in ('Sample', 'Map_name')):
            continue
        ref_fname = str(row['Map_name']).strip()
        seq_prefix = seq_prefix_for(row, seq_format)
        name = Path(seq_prefix).name
        parsed = parse_sample_name(name)
        order = parsed['order'] if parsed else None
        ref_ok = (Path(ref_dir) / ref_fname).exists()
        r1, r2, fastq_error = find_fastq_lanes(seq_dir, seq_prefix)
        records.append({
            'seq_prefix': seq_prefix,
            'name': name,
            'order': order,
            'out_key': out_key(name, order),
            'ref_fname': ref_fname,
            'ref_prefix': Path(ref_fname).stem,
            'ref_ok': ref_ok,
            'r1': r1,
            'r2': r2,
            'fastq_ok': fastq_error is None,
            'fastq_error': fastq_error,
        })
    return records
