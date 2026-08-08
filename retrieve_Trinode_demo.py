#!/Users/yusong/miniforge3/bin/python
"""Download the Trinode demo from USTC Pan and unpack it into ``..``.

The destination is resolved from this script rather than the current working
directory. For a script stored in ``TSI-Denoising/``, the installed demo is
therefore ``../Trinode-demo``. An existing demo directory is left untouched.
"""

from __future__ import annotations

import json
import shutil
import stat
import sys
import tempfile
import zipfile
from pathlib import Path, PurePosixPath
from urllib.request import Request, urlopen

SHARE_ID = "6c21baec311649fc8f47"
SHARE_PAGE_URL = f"https://pan.ustc.edu.cn/share/index/{SHARE_ID}"
SHARE_DETAIL_URL = "https://pan.ustc.edu.cn/api/v1/share/get_share_detail"
TARGET_DIRECTORY = "Trinode-demo"
ARCHIVE_NAME = f"{TARGET_DIRECTORY}.zip.part"
CHUNK_SIZE = 1024 * 1024


def fetch_share_detail() -> tuple[str, int]:
    """Return the current download URL and expected archive size."""
    request = Request(
        f"{SHARE_DETAIL_URL}?share_id={SHARE_ID}",
        headers={"Accept": "application/json"},
    )
    with urlopen(request, timeout=60) as response:
        payload = json.load(response)

    data = payload.get("data", {})
    download_link = data.get("download_link")
    archive_name = data.get("obj_name")
    archive_size = data.get("size")
    can_download = data.get("permissions", {}).get("can_download")

    if payload.get("code") != 1 or not download_link:
        message = payload.get("message", payload)
        raise RuntimeError(f"Unable to get a download link: {message}")
    if data.get("is_dir") or archive_name != f"{TARGET_DIRECTORY}.zip":
        raise RuntimeError(
            f"Unexpected shared object {archive_name!r}; expected "
            f"{TARGET_DIRECTORY}.zip"
        )
    if data.get("is_expired") or can_download is False:
        raise RuntimeError(f"The USTC Pan share is not downloadable: {SHARE_PAGE_URL}")
    if not isinstance(archive_size, int) or archive_size <= 0:
        raise RuntimeError(f"Invalid archive size reported by the share: {archive_size!r}")

    return str(download_link), archive_size


def download(url: str, destination: Path, expected_size: int) -> None:
    """Download the ZIP archive, resuming a valid partial file when possible."""
    downloaded = destination.stat().st_size if destination.exists() else 0

    if downloaded == expected_size:
        print(f"Using completed archive: {destination.name}")
        return
    if downloaded > expected_size:
        print("Existing partial archive is too large; restarting download.")
        destination.unlink()
        downloaded = 0

    request = Request(url)
    if downloaded:
        request.add_header("Range", f"bytes={downloaded}-")
        print(f"Resuming from {downloaded / 1024 / 1024:.1f} MiB ...")
    else:
        print("Downloading Trinode-demo ...")

    with urlopen(request, timeout=120) as response:
        if downloaded and response.status != 206:
            print("Server does not support resume; restarting download.")
            downloaded = 0

        with destination.open("ab" if downloaded else "wb") as output:
            while chunk := response.read(CHUNK_SIZE):
                output.write(chunk)
                downloaded += len(chunk)
                print(
                    f"\rDownloaded {downloaded / 1024 / 1024:.1f} / "
                    f"{expected_size / 1024 / 1024:.1f} MiB "
                    f"({min(downloaded / expected_size * 100, 100):.1f}%)",
                    end="",
                    flush=True,
                )
    print()

    actual_size = destination.stat().st_size
    if actual_size != expected_size:
        raise RuntimeError(
            f"Incomplete archive: expected {expected_size} bytes, got {actual_size}"
        )


def _validated_demo_members(archive: zipfile.ZipFile, staging: Path):
    """Return safe Trinode-demo members while excluding macOS metadata."""
    members: list[zipfile.ZipInfo] = []
    staging_resolved = staging.resolve()

    for member in archive.infolist():
        normalized = member.filename.replace("\\", "/")
        member_path = PurePosixPath(normalized)
        if not member_path.parts or member_path.parts[0] != TARGET_DIRECTORY:
            continue
        if member_path.name == ".DS_Store":
            continue
        if member_path.is_absolute() or ".." in member_path.parts:
            raise RuntimeError(f"Unsafe path in archive: {member.filename}")
        if member.flag_bits & 0x1:
            raise RuntimeError(f"Encrypted archive member is unsupported: {member.filename}")

        file_type = stat.S_IFMT(member.external_attr >> 16)
        if file_type not in {0, stat.S_IFREG, stat.S_IFDIR}:
            raise RuntimeError(f"Unsupported archive member: {member.filename}")

        target = (staging / member_path).resolve()
        if not target.is_relative_to(staging_resolved):
            raise RuntimeError(f"Unsafe path in archive: {member.filename}")
        members.append(member)

    if not members:
        raise RuntimeError(f"Archive does not contain {TARGET_DIRECTORY}/")
    return members


def extract_demo(archive_path: Path, destination_root: Path) -> Path:
    """Safely extract the demo through a temporary staging directory."""
    target = destination_root / TARGET_DIRECTORY
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{TARGET_DIRECTORY}-extract-",
            dir=destination_root,
        )
    )

    try:
        with zipfile.ZipFile(archive_path) as archive:
            members = _validated_demo_members(archive, staging)
            print(f"Extracting {len(members)} archive entries ...")
            archive.extractall(staging, members=members)

        staged_demo = staging / TARGET_DIRECTORY
        if not staged_demo.is_dir():
            raise RuntimeError(f"Archive did not create {TARGET_DIRECTORY}/")
        if target.exists() or target.is_symlink():
            raise FileExistsError(f"Refusing to overwrite existing path: {target}")

        staged_demo.rename(target)
        return target
    finally:
        shutil.rmtree(staging, ignore_errors=True)


def main() -> int:
    script_root = Path(__file__).resolve().parent
    destination_root = script_root.parent
    target = destination_root / TARGET_DIRECTORY

    if target.is_symlink():
        print(f"Refusing to overwrite symbolic link: {target}", file=sys.stderr)
        return 1
    if target.is_dir():
        print(f"Skipping {target}: already exists.")
        return 0
    if target.exists():
        print(f"Refusing to overwrite non-directory path: {target}", file=sys.stderr)
        return 1

    archive_path = script_root / ARCHIVE_NAME
    try:
        download_url, expected_size = fetch_share_detail()
        download(download_url, archive_path, expected_size)
        installed_path = extract_demo(archive_path, destination_root)
    except zipfile.BadZipFile as error:
        archive_path.unlink(missing_ok=True)
        print(f"Invalid ZIP archive; removed partial file: {error}", file=sys.stderr)
        return 1
    except Exception as error:
        print(f"Trinode-demo retrieval failed: {error}", file=sys.stderr)
        return 1

    archive_path.unlink(missing_ok=True)
    print(f"Trinode-demo is ready at {installed_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
