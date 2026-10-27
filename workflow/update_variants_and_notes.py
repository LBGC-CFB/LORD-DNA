#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Update notes and variant databases from Excel exports.

Modes:
  • SNV only:
      --snv-excel SNV.xlsx --snv-notes snv_notes.tsv --snv-db snv_variants_db.tsv
  • SV only:
      --sv-excel  SV.xlsx  --sv-notes  sv_notes.tsv  --sv-db  sv_variants_db.tsv
  • Both in one run (preferred):
      --snv-excel ... --snv-notes ... --snv-db ... --sv-excel ... --sv-notes ... --sv-db ...

Rules:
  - SNV updates use key = ID (e.g., CHR:POS_REF_ALT) for notes **and** variants DB.
  - SV updates use key = ID for notes DB, and key = ID_ALT (= ID + "_" + ALT) for variants DB.
  - Notes dedup uses canonical comparison: no accents, no spaces, uppercase.
  - TSV notes may use column 'notes_variant' or legacy 'notes'. If neither exists, we create 'notes_variant'.
  - Excel may carry 'pré-notes_variant' or legacy 'anciennes_notes'; and 'notes_variant' or legacy 'notes'.
"""

import argparse
import os
import re
import unicodedata
import pandas as pd

# ============== Notes normalization helpers ==============
def strip_accents(s: str) -> str:
    s = unicodedata.normalize("NFKD", s)
    return "".join(ch for ch in s if not unicodedata.combining(ch))

def canon(s) -> str:
    """Canonical form for duplicate checks: no accents, no whitespace, UPPERCASE."""
    if s is None or (isinstance(s, float) and pd.isna(s)):
        return ""
    s = str(s)
    s = strip_accents(s)
    s = re.sub(r"\s+", "", s)
    return s.upper()

def split_notes_cell(cell, sep: str = "_"):
    """Split a notes cell into list of notes. Accept several historical separators."""
    if cell is None or (isinstance(cell, float) and pd.isna(cell)):
        return []
    s = str(cell).strip()
    if not s:
        return []
    s = re.sub(r"[|;,/]+", sep, s)  # normalize alt seps to sep
    parts = [p.strip() for p in re.split(rf"{re.escape(sep)}+", s) if p.strip()]
    return parts

def dedup_notes_list(notes):
    """Order-preserving dedup using canonical keys."""
    seen = set()
    out = []
    for n in notes:
        k = canon(n)
        if k and k not in seen:
            seen.add(k)
            out.append(n.strip())
    return out

def clean_cell(cell, sep: str = "_"):
    """Retroactively clean a notes cell: normalize separators + remove dups."""
    return sep.join(dedup_notes_list(split_notes_cell(cell, sep=sep)))

def _pick_present(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None

# Accept both legacy and new column names in Excel:
EXCEL_PRENOTE_CANDS = ["pré-notes_variant", "anciennes_notes"]
EXCEL_NOTE_CANDS    = ["notes_variant", "notes"]
# Accept both legacy and new column names in TSV:
TSV_NOTE_CANDS      = ["notes_variant", "notes"]

# ============== Variant DB helpers ==============
def _load_variant_db_as_dict(database_file: str):
    """Load TSV with columns: ID, sample_ids (comma-separated). Merge duplicate IDs."""
    if not os.path.exists(database_file):
        return {}
    db_df = pd.read_csv(database_file, sep="\t", dtype=str).fillna("")
    if "ID" not in db_df.columns:
        raise ValueError(f"{database_file} must contain an 'ID' column.")
    if "sample_ids" not in db_df.columns:
        db_df["sample_ids"] = ""
    db_dict = {}
    for _, row in db_df.iterrows():
        vid = str(row["ID"]).strip()
        sids = str(row["sample_ids"]).strip()
        sid_set = set(t.strip() for t in sids.split(",") if t.strip()) if sids else set()
        db_dict.setdefault(vid, set()).update(sid_set)
    return db_dict

def _save_variant_db_from_dict(db_dict: dict, database_file: str):
    out = pd.DataFrame(
        [{"ID": vid, "sample_ids": ",".join(sorted(sids))} for vid, sids in db_dict.items()]
    )
    out = out.sort_values("ID", kind="stable")
    out.to_csv(database_file, sep="\t", index=False)
    print(f"✅ Variant database updated: {database_file} (rows={len(out)})")

def update_variants_db_from_excel(df: pd.DataFrame, database_file: str, id_col: str):
    """Update variants DB using df[id_col] and df['Sample_ID']."""
    if id_col not in df.columns or "Sample_ID" not in df.columns:
        raise ValueError(f"DataFrame must contain '{id_col}' and 'Sample_ID' columns.")
    work = df[[id_col, "Sample_ID"]].copy()
    work[id_col] = work[id_col].astype(str).str.strip()
    work["Sample_ID"] = work["Sample_ID"].astype(str).str.strip()
    work = work[(work[id_col] != "") & (work["Sample_ID"] != "")]
    new_entries = (
        work.groupby(id_col, dropna=False)["Sample_ID"]
            .agg(lambda x: set(v for v in x if v))
            .to_dict()
    )
    db_dict = _load_variant_db_as_dict(database_file)
    for vid, samples in new_entries.items():
        db_dict.setdefault(vid, set()).update(samples)
    _save_variant_db_from_dict(db_dict, database_file)

# ============== Notes updater ==============
def update_notes_from_excel(df_excel: pd.DataFrame, notes_db_file: str, sep: str = "_"):
    """Update a notes TSV keyed by ID. Dedup retroactively; accept multiple column names."""
    if "ID" not in df_excel.columns:
        raise ValueError("Excel must contain 'ID' column.")
    prenote_col = _pick_present(df_excel.columns, EXCEL_PRENOTE_CANDS)
    note_col    = _pick_present(df_excel.columns, EXCEL_NOTE_CANDS)
    if not prenote_col and not note_col:
        raise ValueError(
            f"Excel must contain at least one notes column among: "
            f"{EXCEL_PRENOTE_CANDS + EXCEL_NOTE_CANDS}"
        )
    df1 = df_excel.copy()
    df1["ID"] = df1["ID"].astype(str).str.strip()
    if prenote_col:
        df1[prenote_col] = df1[prenote_col].fillna("").astype(str)
    if note_col:
        df1[note_col] = df1[note_col].fillna("").astype(str)

    # Load or init TSV
    if os.path.exists(notes_db_file):
        df2 = pd.read_csv(notes_db_file, sep="\t", dtype=str)
    else:
        df2 = pd.DataFrame(columns=["ID"])
    tsv_notes_col = _pick_present(df2.columns, TSV_NOTE_CANDS)
    if not tsv_notes_col:
        tsv_notes_col = "notes_variant"  # create preferred name
        df2[tsv_notes_col] = ""
    if "ID" not in df2.columns:
        df2["ID"] = ""

    # Retro-clean existing notes
    df2["ID"] = df2["ID"].astype(str).str.strip()
    df2[tsv_notes_col] = df2[tsv_notes_col].fillna("").apply(lambda c: clean_cell(c, sep=sep))

    # Build updates per ID
    cols_to_collect = [c for c in [prenote_col, note_col] if c]
    df1["__all_notes"] = (
        df1[cols_to_collect]
        .fillna("")
        .agg(lambda s: sep.join([x for x in s if str(x).strip()]), axis=1)
    )
    def collect_notes(series):
        acc = []
        for x in series:
            acc.extend(split_notes_cell(x, sep=sep))
        return dedup_notes_list(acc)
    updates = df1.groupby("ID", dropna=False)["__all_notes"].apply(collect_notes).to_dict()

    # Fast map of first occurrence per ID
    id_to_idx = {}
    for idx, val in df2["ID"].items():
        id_to_idx.setdefault(val, idx)

    total_added, updated_ids, new_ids = 0, 0, 0
    for id_value, new_notes in updates.items():
        if not id_value:
            continue
        if id_value in id_to_idx:
            idx = id_to_idx[id_value]
            existing_list = split_notes_cell(df2.at[idx, tsv_notes_col], sep=sep)
            existing_keys = {canon(x) for x in existing_list}
            added_here = 0
            for n in new_notes:
                key = canon(n)
                if key and key not in existing_keys:
                    existing_list.append(n.strip())
                    existing_keys.add(key)
                    total_added += 1
                    added_here += 1
            df2.at[idx, tsv_notes_col] = sep.join(dedup_notes_list(existing_list))
            if added_here:
                updated_ids += 1
        else:
            cleaned = dedup_notes_list(new_notes)
            df2 = pd.concat(
                [df2, pd.DataFrame([{"ID": id_value, tsv_notes_col: sep.join(cleaned)}])],
                ignore_index=True,
            )
            id_to_idx[id_value] = df2.index[-1]
            new_ids += 1
            total_added += len(cleaned)

    # Collapse duplicate IDs if any already existed
    if df2.duplicated("ID").any():
        merged = (
            df2.groupby("ID", as_index=False)[tsv_notes_col]
               .apply(lambda s: sep.join(dedup_notes_list(sum((split_notes_cell(x, sep=sep) for x in s), []))))
        )
        df2 = merged

    df2 = df2.sort_values("ID", kind="stable")
    df2.to_csv(notes_db_file, sep="\t", index=False)
    print(f"✅ Notes DB updated: {notes_db_file}  |  +IDs:{new_ids}  updIDs:{updated_ids}  +notes:{total_added}")

# ============== Workflow wrappers ==============
def run_snv(snv_excel, snv_notes, snv_db):
    """Update SNV notes (key=ID) and SNV variants DB (key=ID)."""
    df = pd.read_excel(snv_excel, engine="openpyxl", dtype=str)
    update_notes_from_excel(df, snv_notes, sep="_")
    update_variants_db_from_excel(df, snv_db, id_col="ID")

def run_sv(sv_excel, sv_notes, sv_db):
    df = pd.read_excel(sv_excel, engine="openpyxl", dtype=str)
    if "ID" not in df.columns or "ALT" not in df.columns:
        raise ValueError("SV Excel must contain 'ID' and 'ALT' columns.")

    # Notes keyed by original ID
    update_notes_from_excel(df, sv_notes, sep="_")

    # Variants keyed by ID_ALT
    df2 = df.copy()
    df2["ID_ALT"] = df2["ID"].astype(str).str.strip() + "_" + df2["ALT"].astype(str).str.strip()
    update_variants_db_from_excel(df2, sv_db, id_col="ID_ALT")


# ============== CLI ==============
def main():
    ap = argparse.ArgumentParser(description="Update SNV and/or SV notes + variant DBs.")
    ap.add_argument("--snv-excel", help="Excel file for SNV (e.g., filtered_final_snv.xlsx)")
    ap.add_argument("--snv-notes", help="SNV notes TSV (key=ID)")
    ap.add_argument("--snv-db",    help="SNV variants DB TSV (key=ID, columns: ID, sample_ids)")
    ap.add_argument("--sv-excel",  help="Excel file for SV  (e.g., filtered_final_sv.xlsx)")
    ap.add_argument("--sv-notes",  help="SV notes TSV (key=ID)")
    ap.add_argument("--sv-db",     help="SV variants DB TSV (key=ID_ALT stored in 'ID')")
    args = ap.parse_args()

    have_snv = all([args.snv_excel, args.snv_notes, args.snv_db])
    have_sv  = all([args.sv_excel,  args.sv_notes,  args.sv_db])

    if not have_snv and not have_sv:
        ap.error("Provide SNV triplet (--snv-excel/--snv-notes/--snv-db) and/or SV triplet (--sv-excel/--sv-notes/--sv-db).")

    if have_snv:
        print("→ Updating SNV notes and variants DB…")
        run_snv(args.snv_excel, args.snv_notes, args.snv_db)

    if have_sv:
        print("→ Updating SV notes and variants DB…")
        run_sv(args.sv_excel, args.sv_notes, args.sv_db)

    print("✅ Done.")

if __name__ == "__main__":
    main()
