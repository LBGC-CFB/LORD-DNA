#!/usr/bin/env python3
"""Merge Dorado's nested or legacy flat FASTQs into mapped sample FASTQs."""

import argparse
import gzip
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


def compress_fastqs(paths, destination, threads):
    def copy_reads(target):
        size = 0
        for path in paths:
            opener = gzip.open if path.suffix == ".gz" else open
            with opener(path, "rb") as source:
                while chunk := source.read(1024 * 1024):
                    target.write(chunk)
                    size += len(chunk)
        return size

    with destination.open("wb") as output:
        pigz = shutil.which("pigz")
        if pigz:
            with subprocess.Popen([pigz, "-p", str(threads), "-c"], stdin=subprocess.PIPE, stdout=output) as process:
                try:
                    size = copy_reads(process.stdin)
                finally:
                    process.stdin.close()
                if process.wait() != 0:
                    raise RuntimeError("pigz failed while compressing " + str(destination))
                return size
        with gzip.GzipFile(fileobj=output, mode="wb", compresslevel=6) as target:
            return copy_reads(target)


def merge_fastqs(demux_dir, outdir, mapping, samples, threads=1):
    demux_dir, outdir = Path(demux_dir), Path(outdir)
    mapping = {sample: sample for sample in samples} | mapping
    outdir.mkdir(parents=True, exist_ok=True)
    files = {sample: [] for sample in samples}
    sizes = dict.fromkeys(samples, 0)

    # Prepare complete replacements before touching any existing final FASTQ.
    with tempfile.TemporaryDirectory(prefix=".dorado-", dir=outdir) as scratch:
        for path in sorted(demux_dir.rglob("*")):
            if not path.is_file() or not path.name.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                continue
            parts = path.relative_to(demux_dir).parts
            if any("unclassified" in part.lower() for part in parts):
                continue
            barcode = next(
                (match.group() for part in reversed(parts)
                 if (match := re.search(r"barcode[0-9]+", part))),
                None,
            )
            sample = mapping.get(barcode)
            if sample not in files:
                continue
            files[sample].append(path)

        for sample in samples:
            sizes[sample] = compress_fastqs(files[sample], Path(scratch) / f"{sample}.fastq.gz", threads)

        missing = [sample for sample, size in sizes.items() if not size]
        if missing:
            raise ValueError("No demultiplexed reads for configured sample(s): " + ", ".join(missing))
        for sample in samples:
            (Path(scratch) / f"{sample}.fastq.gz").replace(outdir / f"{sample}.fastq.gz")
            print(f"{sample}: merged {len(files[sample])} FASTQ file(s)")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demux_dir")
    parser.add_argument("outdir")
    parser.add_argument("mapping", help="Barcode-to-sample mapping as JSON")
    parser.add_argument("samples", help="Expected sample names as JSON")
    parser.add_argument("--threads", type=int, default=1)
    args = parser.parse_args()
    merge_fastqs(args.demux_dir, args.outdir, json.loads(args.mapping), json.loads(args.samples), args.threads)


if __name__ == "__main__":
    main()
