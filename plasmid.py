#! /usr/bin/env python3

import argparse
import fcntl
import hashlib
import shlex
import shutil
import sys
from pathlib import Path

import snakemake
from snakemake.resources import DefaultResources

sys.path.insert(0, str(Path(__file__).resolve().parent / 'scripts'))
from make_status_table import make_status_table  # noqa: E402


# Каталог прогона внутри -t: <tmp-dir>/sm-plasmid-<хеш пути -o>/.
# Имя постоянно для данного -o: Snakemake видит те же пути промежуточных файлов,
# поэтому повторный запуск ничего не пересчитывает, а после сбоя продолжает работу.
# Пайплайн пишет только в этот каталог и удаляет только его (после успешного прогона).
RUN_PREFIX = 'sm-plasmid-'
RUN_MARKER = '.sm-plasmid-run'  # отличает каталог прогона от постороннего каталога
RUN_LOCK = '.lock'              # один прогон на -o одновременно


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


def run_dir_for(base: Path, out_dir: Path) -> Path:
    return base / f'{RUN_PREFIX}{hashlib.sha256(str(out_dir).encode()).hexdigest()[:12]}'


def open_run_dir(run: Path):
    """Создаёт (или берёт существующий) каталог прогона и блокирует его. Возвращает файл блокировки."""
    if run.exists() and not (run / RUN_MARKER).is_file():
        sys.exit(f'"{run}" exists but is not a sm-plasmid run directory (no {RUN_MARKER})')
    run.mkdir(parents=True, exist_ok=True)
    (run / RUN_MARKER).touch()
    lock = open(run / RUN_LOCK, 'w')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        sys.exit(f'Another run with the same --out-dir is in progress (run directory "{run}")')
    return lock


if __name__ == '__main__':
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-x', '--table', required=True)
    argparser.add_argument('-s', '--seq-dir', required=True, help='FASTQ directory (read-only)')
    argparser.add_argument('-r', '--ref-dir', required=True, help='references directory (read-only)')
    argparser.add_argument('-o', '--out-dir', required=True, help='result files')
    argparser.add_argument('-t', '--tmp-dir', default='/tmp',
                           help=f'base directory for intermediate files: the run works in {RUN_PREFIX}<hash of -o>/ '
                                'inside it, removed after a successful run and kept after a failed one '
                                '(rerunning the same command resumes)')
    argparser.add_argument('-@', '--threads', type=int, default=8)
    argparser.add_argument('-f', '--seq-format', type=int, default=0,
                           help='1: FASTQ in <seq-dir>/<ID>/<Folder>/<Sample>_*; otherwise: <seq-dir>/<Sample>_*')
    argparser.add_argument('-l', '--trim-min-len', type=int, default=25, help='Trimmomatic MINLEN')
    # Вывод инструментов всегда пишется в <out-dir>/logs/...; в stdout — только с -v.
    # Сообщения Snakemake и plasmid.py печатаются всегда.
    argparser.add_argument('-q', '--quiet-tools', dest='quiet_tools', action='store_true', default=True,
                           help='tool output only to log files, not to stdout (default)')
    argparser.add_argument('-v', '--verbose-tools', dest='quiet_tools', action='store_false',
                           help='also print tool output to stdout')
    argparser.add_argument('-d', '--dry-run', action='count', default=0)
    args = argparser.parse_args()

    seq_dir, ref_dir, out_dir, table = (
        Path(p).resolve() for p in (args.seq_dir, args.ref_dir, args.out_dir, args.table))

    base_tmp = Path(args.tmp_dir).resolve()
    check_dirs(seq_dir, ref_dir, out_dir, base_tmp)
    base_created = not base_tmp.exists()
    run_tmp = run_dir_for(base_tmp, out_dir)
    if not args.dry_run:
        run_lock = open_run_dir(run_tmp)  # держится до выхода из процесса

    sys_tmp = run_tmp / 'system'  # $TMPDIR для задач (sort, R, SPAdes и т.п.)
    if not args.dry_run:
        out_dir.mkdir(parents=True, exist_ok=True)
        sys_tmp.mkdir(parents=True, exist_ok=True)
    print(f'Run tmp: {run_tmp}', flush=True)

    config = vars(args)
    config.update(seq_dir=str(seq_dir), ref_dir=str(ref_dir), out_dir=str(out_dir),
                  tmp_dir=str(run_tmp), table=str(table))

    success = False
    try:
        success = snakemake.snakemake(
            snakefile='Snakefile',
            config=config,
            cores=config['threads'],
            use_singularity=True,
            singularity_args=bind_args(ro=[seq_dir, ref_dir, table], rw=[run_tmp, out_dir]),
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
            shutil.rmtree(run_tmp)
            if base_created:
                try:
                    base_tmp.rmdir()  # только если пуст
                except OSError:
                    pass
        else:
            print(f'Pipeline finished with errors; intermediate files are kept in "{run_tmp}".\n'
                  'To resume, run the same command again', file=sys.stderr, flush=True)
    sys.exit(0 if success else 1)