#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
freeze_containers.py

Reads a config.yaml, resolves symlinks in config['singularity'], and writes a JSON
mapping tool -> real .simg/.sif path (a "lockfile" of used versions).

Usage:
  python3 freeze_containers.py \
    --config /path/to/config.yaml \
    --out /path/to/containers.lock.json
"""

import argparse
import json
import os
from pathlib import Path

import yaml


def resolve_path(p: str, base_dir: Path) -> str:
    p = (p or "").strip()
    if not p:
        return p

    path = Path(p)

    # If relative, resolve relative to config file directory
    if not path.is_absolute():
        path = (base_dir / path)

    # realpath() follows symlinks; we still output absolute
    return os.path.realpath(str(path))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True, help="Path to config.yaml")
    ap.add_argument("--out", required=True, help="Path to output JSON lockfile")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    if not cfg_path.exists():
        raise SystemExit(f"[ERROR] config not found: {cfg_path}")

    with cfg_path.open("r", encoding="utf-8") as f:
        cfg = yaml.safe_load(f) or {}

    sing = cfg.get("singularity")
    if not isinstance(sing, dict) or not sing:
        raise SystemExit("[ERROR] config['singularity'] missing or empty")

    base_dir = cfg_path.parent

    locked = {}
    for name, p in sing.items():
        if p is None:
            continue
        locked[str(name)] = resolve_path(str(p), base_dir)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(locked, f, indent=4, ensure_ascii=False)

    print(f"[OK] Wrote {len(locked)} container paths to: {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
