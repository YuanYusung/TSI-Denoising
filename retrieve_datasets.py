#!/Users/yusong/miniforge3/bin/python
"""Download and unpack the public tutorial inputs from the USTC Pan share.

Existing dataset directories are left untouched, so the script is safe to run
again after an interrupted or completed download.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
import tarfile
from pathlib import Path
from urllib.request import Request, urlopen

SHARE_ID = "b63c60ec1c2d46f5bc69"
SHARE_DETAIL_URL = "https://pan.ustc.edu.cn/api/v1/share/get_share_detail"
DATASETS = (
    Path("tutorial/RR_Array/input_public"),
    Path("tutorial/MARS_DAS/input_public"),
)
CHUNK_SIZE = 1024 * 1024


def fetch_download_url() -> str:
    """Return the current download URL exposed by the public share."""
    request = Request(f"{SHARE_DETAIL_URL}?share_id={SHARE_ID}")
    with urlopen(request, timeout=60) as response:
        payload = json.load(response)

    if payload.get("code") != 1 or not payload.get("data", {}).get("download_link"):
        raise RuntimeError(f"Unable to get a download link: {payload.get('message', payload)}")
    return payload["data"]["download_link"]


def download(url: str, destination: Path) -> None:
    """Download an archive, resuming an existing partial download when possible."""
    downloaded = destination.stat().st_size if destination.exists() else 0
    request = Request(url)

    if downloaded:
        request.add_header("Range", f"bytes={downloaded}-")
        print(f"Resuming from {downloaded / 1024 / 1024:.1f} MiB ...")
    else:
        print("Downloading public tutorial inputs ...")

    with urlopen(request, timeout=120) as response:
        if downloaded and response.status != 206:
            print("Server does not support resume; restarting download.")
            downloaded = 0

        content_range = response.headers.get("Content-Range", "")
        total: int | None = None
        if "/" in content_range:
            total_text = content_range.rpartition("/")[2]
            if total_text.isdigit():
                total = int(total_text)
        elif response.headers.get("Content-Length", "").isdigit():
            total = downloaded + int(response.headers["Content-Length"])

        with destination.open("ab" if downloaded else "wb") as output:
            while chunk := response.read(CHUNK_SIZE):
                output.write(chunk)
                downloaded += len(chunk)

                if total:
                    print(
                        f"\rDownloaded {downloaded / 1024 / 1024:.1f} / "
                        f"{total / 1024 / 1024:.1f} MiB "
                        f"({min(downloaded / total * 100, 100):.1f}%)",
                        end="",
                        flush=True,
                    )
    print()


def extract_dataset(archive: tarfile.TarFile, root: Path, dataset: Path) -> None:
    """Extract exactly one dataset directory, rejecting unsafe archive members."""
    prefix = f"{dataset.as_posix()}/"
    members = [member for member in archive.getmembers() if member.name == prefix[:-1] or member.name.startswith(prefix)]
    if not members:
        raise RuntimeError(f"Archive does not contain {dataset}")

    root_resolved = root.resolve()
    for member in members:
        if member.issym() or member.islnk() or member.isdev():
            raise RuntimeError(f"Unsupported archive member: {member.name}")
        target = (root / member.name).resolve()
        if not target.is_relative_to(root_resolved):
            raise RuntimeError(f"Unsafe path in archive: {member.name}")

    print(f"Extracting {dataset} ...")
    archive.extractall(root, members=members)

def remove_hidden_entries(dataset_dir: Path, root: Path) -> None:
    """Remove dot-prefixed files and directories under a public input directory."""
    if not dataset_dir.is_dir():
        return

    removed = 0
    for current_path, directory_names, file_names in os.walk(
        dataset_dir, topdown=True, followlinks=False
    ):
        current = Path(current_path)

        for name in file_names:
            if name.startswith("."):
                (current / name).unlink()
                removed += 1

        for name in directory_names[:]:
            if not name.startswith("."):
                continue

            entry = current / name
            if entry.is_symlink():
                entry.unlink()
            else:
                shutil.rmtree(entry)
            directory_names.remove(name)
            removed += 1

    if removed:
        print(f"Removed {removed} hidden item(s) from {dataset_dir.relative_to(root)}.")


def main() -> int:
    root = Path(__file__).resolve().parent
    missing: list[Path] = []

    for dataset in DATASETS:
        target = root / dataset
        if target.is_dir():
            print(f"Skipping {dataset}: already exists.")
        elif target.exists():
            print(f"Refusing to overwrite non-directory path: {target}", file=sys.stderr)
            return 1
        else:
            missing.append(dataset)

    if not missing:
        print("Both public tutorial inputs are already present; nothing to download.")
        return 0

    archive_path = root / "public_tutorial_inputs.tar.gz.part"
    completed = False

    try:
        download(fetch_download_url(), archive_path)
        with tarfile.open(archive_path, "r:*") as archive:
            for dataset in missing:
                extract_dataset(archive, root, dataset)
        completed = True
    except Exception as error:
        print(f"Download failed: {error}", file=sys.stderr)
        return 1
    finally:
        # 仅在下载和解压成功后删除；失败时保留 .part 文件供下次续传。
        if completed:
            archive_path.unlink(missing_ok=True)

    for dataset in DATASETS:
        remove_hidden_entries(root / dataset, root)

    print("Public tutorial input download completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
