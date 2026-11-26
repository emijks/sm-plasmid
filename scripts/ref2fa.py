#!/usr/bin/env python3
"""Конвертация референса (SnapGene .dna / GenBank / FASTA) в FASTA для bwa-mem2 и samtools.

Использование:
    ref2fasta.py <in> <out.fa> [--format auto|snapgene|genbank|fasta]

Формат по умолчанию определяется по расширению (см. FORMATS).

Имена последовательностей в выходном FASTA:
  snapgene, genbank — <имя файла>_<номер записи>, например pX.dna_0, pX.gbk_0
                      (так же, как в прежнем dna2fasta.py);
  fasta             — исходные ID сохраняются (как при прежнем cp в fa2fa).

Проверки (ошибка и ненулевой код возврата):
  * нет ни одной записи или запись с пустой последовательностью;
  * символы вне алфавита IUPAC для ДНК;
  * повторяющиеся ID (samtools faidx на них падает).

Модуль импортируется из Snakefile ради FORMATS / ref_format(): Biopython
импортируется лениво, внутри convert(), чтобы окружению Snakemake он был не нужен.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Расширение (в нижнем регистре) -> формат Biopython. Единый список для Snakefile.
FORMATS = {
    '.dna': 'snapgene',
    '.gbk': 'genbank',
    '.gb': 'genbank',
    '.gbff': 'genbank',
    '.genbank': 'genbank',
    '.fa': 'fasta',
    '.fasta': 'fasta',
    '.fna': 'fasta',
    '.fas': 'fasta',
}

# Однобуквенные коды IUPAC для нуклеотидов (U допускаем для РНК-записей).
IUPAC_DNA = set('ACGTURYSWKMBDHVN')

FASTA_LINE_WIDTH = 60


def ref_format(path: str) -> str | None:
    """Формат по расширению файла или None, если расширение не поддерживается."""
    return FORMATS.get(Path(path).suffix.lower())


def convert(src: str, dst: str, fmt: str) -> list[tuple[str, int]]:
    """Читает src в формате fmt, проверяет и пишет FASTA в dst. Возвращает [(id, длина)]."""
    from Bio import SeqIO
    from Bio.SeqRecord import SeqRecord

    fname = Path(src).name
    records, errors, seen = [], [], set()

    try:
        parsed = list(SeqIO.parse(src, fmt))
    except Exception as e:  # битый файл: Biopython бросает ValueError/AssertionError без имени файла
        raise ValueError(f'{src} [{fmt}]: cannot parse: {type(e).__name__}: {e}') from e

    for i, rec in enumerate(parsed):
        rec_id = rec.id if fmt == 'fasta' else f'{fname}_{i}'
        description = rec.description if fmt == 'fasta' else ''
        try:
            seq = str(rec.seq).upper()
        except Exception as e:  # GenBank без ORIGIN: последовательность не определена
            errors.append(f'record #{i} ({rec.id}): no sequence ({e})')
            continue

        if not seq:
            errors.append(f'record #{i} ({rec_id}): empty sequence')
        bad = sorted(set(seq) - IUPAC_DNA)
        if bad:
            errors.append(f'record #{i} ({rec_id}): non-IUPAC characters: {"".join(bad)!r}')
        if rec_id in seen:
            errors.append(f'record #{i}: duplicate ID {rec_id!r}')
        seen.add(rec_id)

        # Новый SeqRecord без аннотаций: в FASTA они не нужны, а описание
        # SnapGene/GenBank иначе попало бы в заголовок.
        records.append(SeqRecord(rec.seq.upper(), id=rec_id, name=rec_id, description=description))

    if not records and not errors:
        errors.append('no records found (wrong format or empty file?)')
    if errors:
        raise ValueError(f'{src} [{fmt}]:\n' + '\n'.join(f'  {e}' for e in errors))

    with open(dst, 'w') as out:
        for rec in records:
            header = rec.id if not rec.description or rec.description == rec.id else rec.description
            if not header.startswith(rec.id):
                header = f'{rec.id} {header}'
            out.write(f'>{header}\n')
            seq = str(rec.seq)
            for j in range(0, len(seq), FASTA_LINE_WIDTH):
                out.write(seq[j:j + FASTA_LINE_WIDTH] + '\n')

    return [(rec.id, len(rec.seq)) for rec in records]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('src', help='входной файл референса')
    parser.add_argument('dst', help='выходной FASTA')
    parser.add_argument('--format', default='auto', choices=['auto', *sorted(set(FORMATS.values()))],
                        help='формат входного файла (по умолчанию — по расширению)')
    args = parser.parse_args()

    fmt = ref_format(args.src) if args.format == 'auto' else args.format
    if fmt is None:
        sys.exit(f'{args.src}: unsupported extension, expected one of: {", ".join(FORMATS)} '
                 f'(or pass --format)')

    try:
        stats = convert(args.src, args.dst, fmt)
    except ValueError as e:
        Path(args.dst).unlink(missing_ok=True)
        sys.exit(str(e))

    print(f'Converted [{fmt}]: {args.src} -> {args.dst}')
    for rec_id, length in stats:
        print(f'  {rec_id}: {length} bp')


if __name__ == '__main__':
    main()
