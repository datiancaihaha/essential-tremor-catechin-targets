from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import requests


CHUNK_SIZE = 8 * 1024 * 1024


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_text(url: str) -> str:
    response = requests.get(url, timeout=60)
    response.raise_for_status()
    return response.text


def expected_md5(base_url: str, filename: str) -> str:
    text = fetch_text(f"{base_url}/md5sum.txt")
    pattern = re.compile(rf"^([0-9a-f]{{32}})\s+{re.escape(filename)}$", re.MULTILINE)
    match = pattern.search(text)
    if not match:
        raise ValueError(f"No checksum found for {filename}")
    return match.group(1)


def download_with_resume(url: str, destination: Path, expected_size: int) -> None:
    part = destination.with_suffix(destination.suffix + ".part")
    completed = part.stat().st_size if part.exists() else 0
    if destination.exists() and destination.stat().st_size == expected_size:
        return
    if completed > expected_size:
        raise ValueError(f"Partial file exceeds expected size: {part}")

    headers = {"Range": f"bytes={completed}-"} if completed else {}
    mode = "ab" if completed else "wb"
    with requests.get(url, headers=headers, stream=True, timeout=(30, 300)) as response:
        response.raise_for_status()
        if completed and response.status_code != 206:
            raise RuntimeError(f"Server did not honor resume request for {url}")
        with part.open(mode) as handle:
            for chunk in response.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    handle.write(chunk)
    if part.stat().st_size != expected_size:
        raise ValueError(
            f"Downloaded size mismatch for {destination.name}: "
            f"{part.stat().st_size} != {expected_size}"
        )
    os.replace(part, destination)


def validate_gzip(path: Path) -> int:
    line_count = 0
    with gzip.open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            line_count += chunk.count(b"\n")
    return line_count


def process_row(row: dict[str, str], output_root: Path) -> dict[str, Any]:
    accession = row["accession"]
    url = row["url"]
    expected_size = int(row["content_length_bytes"])
    base_url = url.rsplit("/", 1)[0]
    filename = url.rsplit("/", 1)[1]
    accession_dir = output_root / accession
    accession_dir.mkdir(parents=True, exist_ok=True)
    destination = accession_dir / filename
    metadata_name = f"{filename}-meta.yaml"
    metadata_path = accession_dir / metadata_name
    checksum_path = accession_dir / "md5sum.txt"

    checksum_text = fetch_text(f"{base_url}/md5sum.txt")
    checksum_path.write_text(checksum_text, encoding="utf-8")
    metadata_path.write_text(fetch_text(f"{base_url}/{metadata_name}"), encoding="utf-8")
    checksum = expected_md5(base_url, filename)
    download_with_resume(url, destination, expected_size)
    observed = file_md5(destination)
    if observed != checksum:
        raise ValueError(f"MD5 mismatch for {accession}: {observed} != {checksum}")
    line_count = validate_gzip(destination)
    expected_rows = int(row["expected_rows"])
    data_rows = max(line_count - 1, 0)
    if data_rows != expected_rows:
        raise ValueError(
            f"Row-count mismatch for {accession}: {data_rows} != {expected_rows}"
        )

    result = dict(row)
    result.update(
        {
            "expected_md5": checksum,
            "observed_md5": observed,
            "observed_size_bytes": destination.stat().st_size,
            "data_rows": data_rows,
            "gzip_crc": "PASS",
            "integrity_status": "PASS",
            "local_path": str(destination),
        }
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("output_root", type=Path)
    parser.add_argument("result_manifest", type=Path)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()

    with args.manifest.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    args.output_root.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_row, row, args.output_root): row for row in rows}
        for future in as_completed(futures):
            row = futures[future]
            accession = row["accession"]
            result = future.result()
            results.append(result)
            print(f"{accession}: PASS ({result['data_rows']} rows)", flush=True)

    results.sort(key=lambda item: item["accession"])
    fieldnames = list(results[0].keys())
    args.result_manifest.parent.mkdir(parents=True, exist_ok=True)
    with args.result_manifest.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
