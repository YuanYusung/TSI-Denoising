#!/Users/yusong/miniforge3/bin/python
"""Download and unpack the tutorial raw datasets from the USTC Pan share.

Existing dataset directories are left untouched, so the script is safe to run
again after an interrupted or completed download.
"""

from __future__ import annotations

import json
import sys
import tarfile
import tempfile
from pathlib import Path
from urllib.request import Request, urlopen


SHARE_ID = "26cffdd7fe3947babbfb"
SHARE_DETAIL_URL = "https://pan.ustc.edu.cn/api/v1/share/get_share_detail"
DATASETS = (
    Path("tutorial/RR_Array/input"),
    Path("tutorial/MARS_DAS/input"),
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
    """Download *url* to *destination* without loading the archive into memory."""
    print("Downloading raw_data_inputs.tar ...")
    with urlopen(Request(url), timeout=120) as response, destination.open("wb") as output:
        while chunk := response.read(CHUNK_SIZE):
            output.write(chunk)


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
        print("Both raw datasets are already present; nothing to download.")
        return 0

    archive_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix="raw_data_inputs_", suffix=".tar", dir=root, delete=False
        ) as temporary_file:
            archive_path = Path(temporary_file.name)

        download(fetch_download_url(), archive_path)
        with tarfile.open(archive_path, "r") as archive:
            for dataset in missing:
                extract_dataset(archive, root, dataset)
    except Exception as error:
        print(f"Download failed: {error}", file=sys.stderr)
        return 1
    finally:
        if archive_path is not None:
            archive_path.unlink(missing_ok=True)

    print("Raw data download completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
