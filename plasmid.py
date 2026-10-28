#! /usr/bin/env python3

import argparse
import shlex
import shutil
import sys
from pathlib import Path

import snakemake
from snakemake.resources import DefaultResources

sys.path.insert(0, str(Path(__file__).resolve().parent / 'scripts'))
from make_status_table import make_status_table  # noqa: E402


# Подкаталоги tmp-dir, которыми владеет пайплайн (удаляются после успешного прогона)
TMP_SUBDIRS = ('refs', 'samples', 'system')


def check_dirs(seq_dir: Path, ref_dir: Path, out_dir: Path, tmp_dir: Path) -> None:
    """-s/-r только читаются; -o и -t не должны пересекаться ни с ними, ни друг с другом."""
    def nested(a: Path, b: Path) -> bool:
        return a == b or a.is_relative_to(b) or b.is_relative_to(a)

    pairs = [('--out-dir', out_dir, '--seq-dir', seq_dir), ('--out-dir', out_dir, '--ref-dir', ref_dir),
             ('--tmp-dir', tmp_dir, '--seq-dir', seq_dir), ('--tmp-dir', tmp_dir, '--ref-dir', ref_dir),
             ('--tmp-dir', tmp_dir, '--out-dir', out_dir)]
    for name_a, a, name_b, b in pairs:
        if nested(a, b):
            sys.exit(f'{name_a} "{a}" and {name_b} "{b}" must not be the same or nested')


def bind_args(ro: list, rw: list) -> str:
    """Аргументы -B для singularity без повторов вложенных путей."""
    binds, seen = [], []
    for path, mode in [(p, 'ro') for p in ro] + [(p, 'rw') for p in rw]:
        if any(path == s or path.is_relative_to(s) for s in seen):
            continue
        seen.append(path)
        binds.append(f'-B {shlex.quote(f"{path}:{path}:{mode}")}')
    return ' '.join(binds)


def cleanup_tmp(tmp_dir: Path, created: bool) -> None:
    for sub in TMP_SUBDIRS:
        shutil.rmtree(tmp_dir / sub, ignore_errors=True)
    if created:
        try:
            tmp_dir.rmdir()  # только если пуст
        except OSError:
            pass


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-x', '--table', required=True)
    argparser.add_argument('-s', '--seq-dir', required=True, help='FASTQ directory (read-only)')
    argparser.add_argument('-r', '--ref-dir', required=True, help='references directory (read-only)')
    argparser.add_argument('-o', '--out-dir', required=True, help='result files')
    argparser.add_argument('-t', '--tmp-dir', default='./tmp',
                           help='intermediate files; removed after a successful run, kept after a failed one')
    argparser.add_argument('-@', '--threads', type=int, default=8)
    argparser.add_argument('-f', '--seq-format', type=int, default=0,
                           help='1: FASTQ in <seq-dir>/<ID>/<Folder>/<Sample>_*; otherwise: <seq-dir>/<Sample>_*')
    argparser.add_argument('-l', '--trim-min-len', type=int, default=25, help='Trimmomatic MINLEN')
    argparser.add_argument('-q', '--quiet-tools', action='store_true',
                           help='write tool output only to log files (<out-dir>/logs/...), not to stdout; '
                                'Snakemake and pipeline messages are always printed')
    argparser.add_argument('-d', '--dry-run', action='count', default=0)
    args = argparser.parse_args()

    seq_dir, ref_dir, out_dir, tmp_dir, table = (
        Path(p).resolve() for p in (args.seq_dir, args.ref_dir, args.out_dir, args.tmp_dir, args.table))
    check_dirs(seq_dir, ref_dir, out_dir, tmp_dir)

    tmp_created = not tmp_dir.exists()
    sys_tmp = tmp_dir / 'system'  # $TMPDIR для задач (sort, R, SPAdes и т.п.)
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        sys_tmp.mkdir(parents=True, exist_ok=True)

    config = vars(args)
    config.update(seq_dir=str(seq_dir), ref_dir=str(ref_dir), out_dir=str(out_dir),
                  tmp_dir=str(tmp_dir), table=str(table))

    success = False
    try:
        success = snakemake.snakemake(
            snakefile='Snakefile',
            config=config,
            cores=config['threads'],
            use_singularity=True,
            singularity_args=bind_args(ro=[seq_dir, ref_dir, table], rw=[tmp_dir, out_dir]),
            default_resources=DefaultResources([f'tmpdir={str(sys_tmp)!r}']),
            dryrun=bool(args.dry_run),
            keepgoing=True,
            keep_incomplete=True,
            force_incomplete=True,
            printshellcmds=False,
            quiet=True if bool(args.dry_run) else False,
            # rerun_triggers=['input'],
        )
    finally:
        if not args.dry_run:
            make_status_table(
                table=str(table),
                seq_dir=str(seq_dir),
                ref_dir=str(ref_dir),
                seq_format=args.seq_format,
                out_dir=str(out_dir),
                output=str(out_dir / 'pipeline_status.tsv'),
            )

    if not args.dry_run:
        if success:
            cleanup_tmp(tmp_dir, tmp_created)
        else:
            print(f'Pipeline finished with errors; intermediate files are kept in "{tmp_dir}"', file=sys.stderr)
    sys.exit(0 if success else 1)
