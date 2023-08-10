import os
import pathlib

from ftplib import FTP
from typing import Callable, Iterable

from rich.progress import track
from unsync import unsync

from ensembl_cli.util import (
    atomic_write,
    dont_checksum,
    get_sig_calc_func,
    get_signature_data,
    is_signature,
)


def configured_ftp(host: str = "ftp.ensembl.org") -> FTP:
    ftp = FTP(host)
    ftp.login()
    return ftp


def listdir(host: str, path: str, pattern: Callable = None):
    """returns directory listing"""
    pattern = pattern or (lambda x: True)
    ftp = configured_ftp(host=host)
    ftp.cwd(path)
    for fn in ftp.nlst():
        if pattern(fn):
            yield f"{path}/{fn}"
    ftp.close()


@unsync
def unsynced_copy_to_local(
    host: str, src: os.PathLike, dest: os.PathLike
) -> os.PathLike:
    if dest.exists():
        return dest
    ftp = configured_ftp(host=host)
    # pass in checksum and keep going until it's correct?
    with atomic_write(dest, mode="wb") as outfile:
        ftp.retrbinary(f"RETR {src}", outfile.write)

    ftp.close()
    return dest


def download_data(
    *,
    host: str,
    local_dest: os.PathLike,
    remote_paths: Iterable[os.PathLike],
    description,
    do_checksum: bool,
) -> bool:
    tasks = [
        unsynced_copy_to_local(host, path, local_dest / pathlib.Path(path).name)
        for path in remote_paths
    ]
    saved_paths = [
        task.result() for task in track(tasks, description=description, transient=True)
    ]

    # load the signature data and sig calc keyed by parent dir
    all_checksums = {}
    all_check_funcs = {}
    for path in saved_paths:
        if is_signature(path):
            all_checksums[str(path.parent)] = get_signature_data(path)
            all_check_funcs[str(path.parent)] = get_sig_calc_func(path.name)

    if do_checksum:
        for path in track(
            saved_paths, description="Validating checksums...", transient=True
        ):
            if dont_checksum(path):
                continue
            key = str(path.parent)
            expect_sig = all_checksums[key][path.name]
            calc_sig = all_check_funcs[key]
            signature = calc_sig(path.read_bytes(), path.stat().st_size)
            assert signature == expect_sig, path

    return True
