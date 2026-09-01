#!/usr/bin/env python3
"""Prepare the 1000G EUR LD panel, clump MiBioGen instruments, and audit ET overlap."""

from __future__ import annotations

import concurrent.futures
import csv
import hashlib
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import duckdb


ROOT = Path(__file__).resolve().parents[1]
MBG = ROOT / "raw" / "mibiogen" / "MBG.allHits.p1e4.ranged.txt"
ET = ROOT / "raw" / "et_gwas" / "extracted_v1" / "G250_Essential_tremor_summary"
PANEL = ROOT / "raw" / "ld_reference" / "1kg_v3_eur" / "EUR"
PANEL_ARCHIVE = ROOT / "raw" / "ld_reference" / "downloads" / "1kg.v3.tgz"
PANEL_DOWNLOAD_META = PANEL_ARCHIVE.with_suffix(PANEL_ARCHIVE.suffix + ".download.json")
PLINK = ROOT / "tools" / "plink_1.9_20250819" / "plink.exe"
PLINK_ARCHIVE = ROOT / "raw" / "ld_reference" / "downloads" / "plink_win64_20250819.zip"
OUT = ROOT / "outputs"
QA = ROOT / "qa"
LOG = ROOT / "logs"
WORK = ROOT / "work" / "stage1_ld_clump"
SUBSET = WORK / "EUR_candidate_union"

P_THRESHOLD = 1e-5
CLUMP_R2 = 0.001
CLUMP_KB = 10000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def valid_rsid(value: str) -> bool:
    return bool(value) and value not in {".", "NA"}


def complement(allele: str) -> str | None:
    allele = (allele or "").upper()
    if len(allele) != 1 or allele not in "ACGT":
        return None
    return allele.translate(str.maketrans("ACGT", "TGCA"))


def allele_compatibility(ea: str, oa: str, a1: str, a2: str) -> str:
    ea, oa, a1, a2 = ea.upper(), oa.upper(), a1.upper(), a2.upper()
    if (ea, oa) == (a1, a2):
        return "exact"
    if (ea, oa) == (a2, a1):
        return "swapped"
    cea, coa = complement(ea), complement(oa)
    if cea is not None and (cea, coa) == (a1, a2):
        return "strand_exact"
    if cea is not None and (cea, coa) == (a2, a1):
        return "strand_swapped"
    return "incompatible"


def is_palindromic(ea: str, oa: str) -> bool:
    return len(ea) == len(oa) == 1 and {ea.upper(), oa.upper()} in ({"A", "T"}, {"C", "G"})


def read_candidates() -> list[dict]:
    rows = []
    with MBG.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for source_row_id, row in enumerate(reader, start=1):
            p_value = float(row["P.weightedSumZ"])
            if p_value >= P_THRESHOLD:
                continue
            beta = float(row["beta"])
            se = float(row["SE"])
            rows.append(
                {
                    "candidate_row_id": len(rows) + 1,
                    "source_row_id": source_row_id,
                    "bac": row["bac"],
                    "exposure_chr": row["chr"],
                    "exposure_bp_grch37_inferred": int(row["bp"]),
                    "source_rsid": row["rsID"],
                    "exposure_oa": row["ref.allele"].upper(),
                    "exposure_ea": row["eff.allele"].upper(),
                    "beta_exposure": beta,
                    "se_exposure": se,
                    "p_exposure": p_value,
                    "n_exposure": int(float(row["N"])),
                    "f_statistic": (beta / se) ** 2,
                }
            )
    if len(rows) != 14587:
        raise ValueError(f"Unexpected P<1e-5 candidate count: {len(rows)}")
    return rows


def scan_panel(candidates: list[dict]) -> tuple[dict, dict, dict]:
    wanted_ids = {row["source_rsid"] for row in candidates if valid_rsid(row["source_rsid"])}
    wanted_positions = {
        (row["exposure_chr"], row["exposure_bp_grch37_inferred"])
        for row in candidates
        if not valid_rsid(row["source_rsid"])
    }
    by_id: dict[str, list[dict]] = defaultdict(list)
    by_position: dict[tuple[str, int], list[dict]] = defaultdict(list)
    chromosome_counts = Counter()
    variant_rows = 0
    known_build_rows = []
    with PANEL.with_suffix(".bim").open("r", encoding="utf-8", newline="") as handle:
        for line in handle:
            chrom, rsid, cm, bp, a1, a2 = line.rstrip("\r\n").split()
            variant_rows += 1
            chromosome_counts[chrom] += 1
            record = {
                "panel_chr": chrom,
                "panel_rsid": rsid,
                "panel_cm": cm,
                "panel_bp_grch37": int(bp),
                "panel_a1": a1.upper(),
                "panel_a2": a2.upper(),
            }
            if rsid in wanted_ids:
                by_id[rsid].append(record)
            if (chrom, int(bp)) in wanted_positions:
                by_position[(chrom, int(bp))].append(record)
            if rsid == "rs182549":
                known_build_rows.append(record)
    sample_rows = sum(1 for _ in PANEL.with_suffix(".fam").open("r", encoding="utf-8"))
    bed_bytes = PANEL.with_suffix(".bed").stat().st_size
    expected_bed_bytes = 3 + ((sample_rows + 3) // 4) * variant_rows
    if bed_bytes != expected_bed_bytes:
        raise ValueError(f"BED size mismatch: {bed_bytes} vs expected {expected_bed_bytes}")
    if not any(
        row["panel_chr"] == "2" and row["panel_bp_grch37"] == 136616754
        for row in known_build_rows
    ):
        raise ValueError(f"GRCh37 sentinel rs182549 not found at expected position: {known_build_rows}")
    metrics = {
        "samples": sample_rows,
        "variants": variant_rows,
        "bed_bytes": bed_bytes,
        "bed_expected_bytes": expected_bed_bytes,
        "bed_size_consistent": True,
        "chromosome_variant_counts": dict(sorted(chromosome_counts.items(), key=lambda x: int(x[0]))),
        "grch37_sentinel": known_build_rows,
    }
    return by_id, by_position, metrics


def panel_gate(candidates: list[dict], by_id: dict, by_position: dict) -> list[dict]:
    output = []
    for row in candidates:
        source_rsid = row["source_rsid"]
        recovered = False
        if valid_rsid(source_rsid):
            records = by_id.get(source_rsid, [])
            if not records:
                output.append({**row, "panel_gate_status": "absent_from_panel", "clump_input": "no"})
                continue
            coordinate_matches = [
                record
                for record in records
                if record["panel_chr"] == row["exposure_chr"]
                and record["panel_bp_grch37"] == row["exposure_bp_grch37_inferred"]
            ]
            if not coordinate_matches:
                output.append(
                    {
                        **row,
                        "panel_id_match_count": len(records),
                        "panel_gate_status": "rsid_coordinate_mismatch",
                        "clump_input": "no",
                    }
                )
                continue
        else:
            coordinate_matches = by_position.get(
                (row["exposure_chr"], row["exposure_bp_grch37_inferred"]), []
            )
            recovered = True
            if not coordinate_matches:
                output.append(
                    {**row, "panel_gate_status": "placeholder_not_recovered", "clump_input": "no"}
                )
                continue
        compatible = []
        for record in coordinate_matches:
            orientation = allele_compatibility(
                row["exposure_ea"], row["exposure_oa"], record["panel_a1"], record["panel_a2"]
            )
            if orientation != "incompatible":
                compatible.append((record, orientation))
        if len(compatible) != 1:
            output.append(
                {
                    **row,
                    "panel_coordinate_match_count": len(coordinate_matches),
                    "panel_compatible_match_count": len(compatible),
                    "panel_gate_status": (
                        "panel_allele_incompatible" if not compatible else "ambiguous_panel_match"
                    ),
                    "clump_input": "no",
                }
            )
            continue
        record, orientation = compatible[0]
        output.append(
            {
                **row,
                **record,
                "analysis_rsid": record["panel_rsid"],
                "panel_allele_orientation": orientation,
                "panel_gate_status": (
                    "eligible_recovered_from_chr_bp_alleles" if recovered else "eligible"
                ),
                "clump_input": "yes",
            }
        )
    return output


def run_plink(arguments: list[str], label: str) -> subprocess.CompletedProcess:
    result = subprocess.run(
        [str(PLINK), *arguments], capture_output=True, text=True, encoding="utf-8", errors="replace"
    )
    (LOG / f"{label}.stdout.txt").write_text(
        result.stdout + "\n---STDERR---\n" + result.stderr, encoding="utf-8"
    )
    if result.returncode != 0:
        raise RuntimeError(f"PLINK failed for {label} with exit code {result.returncode}")
    return result


def prepare_subset(gated_rows: list[dict]) -> dict[str, dict]:
    eligible = [row for row in gated_rows if row["clump_input"] == "yes"]
    union_rsids = sorted({row["analysis_rsid"] for row in eligible})
    union_path = WORK / "candidate_union_rsids.txt"
    union_path.write_text("\n".join(union_rsids) + "\n", encoding="ascii")
    run_plink(
        [
            "--bfile", str(PANEL), "--extract", str(union_path), "--make-bed",
            "--out", str(SUBSET), "--allow-no-sex",
        ],
        "stage1_panel_subset",
    )
    subset_ids = {}
    with SUBSET.with_suffix(".bim").open("r", encoding="utf-8") as handle:
        for line in handle:
            chrom, rsid, cm, bp, a1, a2 = line.split()
            subset_ids[rsid] = {
                "panel_chr": chrom,
                "panel_bp_grch37": int(bp),
                "panel_a1": a1.upper(),
                "panel_a2": a2.upper(),
            }
    if set(subset_ids) != set(union_rsids):
        raise ValueError(
            f"Subset ID mismatch: requested {len(union_rsids)}, obtained {len(subset_ids)}"
        )
    run_plink(
        ["--bfile", str(SUBSET), "--freq", "--out", str(SUBSET), "--allow-no-sex"],
        "stage1_panel_subset_freq",
    )
    frequencies = {}
    with SUBSET.with_suffix(".frq").open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=" ", skipinitialspace=True)
        for row in reader:
            frequencies[row["SNP"]] = {
                "panel_freq_a1": row["A1"].upper(),
                "panel_freq_a2": row["A2"].upper(),
                "panel_maf": float(row["MAF"]),
                "panel_nchrobs": int(row["NCHROBS"]),
            }
    for rsid, record in subset_ids.items():
        record.update(frequencies[rsid])
    return subset_ids


def unique_clump_inputs(gated_rows: list[dict]) -> tuple[list[dict], list[dict]]:
    groups: dict[tuple[str, str], list[dict]] = defaultdict(list)
    for row in gated_rows:
        if row["clump_input"] == "yes":
            groups[(row["bac"], row["analysis_rsid"])].append(row)
    selected = []
    duplicates = []
    for key, rows in groups.items():
        ordered = sorted(rows, key=lambda row: (row["p_exposure"], row["candidate_row_id"]))
        selected.append(ordered[0])
        for duplicate in ordered[1:]:
            duplicates.append(
                {
                    "bac": key[0],
                    "analysis_rsid": key[1],
                    "kept_candidate_row_id": ordered[0]["candidate_row_id"],
                    "removed_candidate_row_id": duplicate["candidate_row_id"],
                    "reason": "duplicate_bac_analysis_rsid_keep_lowest_p",
                }
            )
    return selected, duplicates


def parse_clumped(path: Path) -> set[str]:
    if not path.exists() or not path.read_text(encoding="utf-8", errors="replace").strip():
        return set()
    lines = [line.split() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    header = lines[0]
    snp_index = header.index("SNP")
    return {line[snp_index] for line in lines[1:]}


def run_taxon_clumping(rows: list[dict]) -> tuple[set[tuple[str, str]], list[dict]]:
    by_taxon: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_taxon[row["bac"]].append(row)
    clump_dir = WORK / "taxa"
    clump_dir.mkdir(parents=True, exist_ok=True)
    jobs = []
    mapping = []
    for index, (bac, items) in enumerate(sorted(by_taxon.items()), start=1):
        slug = f"taxon_{index:03d}"
        assoc = clump_dir / f"{slug}.assoc.txt"
        with assoc.open("w", encoding="ascii", newline="") as handle:
            handle.write("SNP P\n")
            for row in sorted(items, key=lambda x: (x["p_exposure"], x["analysis_rsid"])):
                handle.write(f"{row['analysis_rsid']} {row['p_exposure']:.17g}\n")
        prefix = clump_dir / slug
        mapping.append({"taxon_index": index, "taxon_slug": slug, "bac": bac, "input_rows": len(items)})
        jobs.append((bac, assoc, prefix))

    def execute(job: tuple[str, Path, Path]) -> tuple[str, set[str]]:
        bac, assoc, prefix = job
        result = subprocess.run(
            [
                str(PLINK), "--bfile", str(SUBSET), "--clump", str(assoc),
                "--clump-snp-field", "SNP", "--clump-field", "P",
                "--clump-p1", str(P_THRESHOLD), "--clump-p2", str(P_THRESHOLD),
                "--clump-r2", str(CLUMP_R2), "--clump-kb", str(CLUMP_KB),
                "--out", str(prefix), "--allow-no-sex",
            ],
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        (prefix.with_suffix(".stdout.txt")).write_text(
            result.stdout + "\n---STDERR---\n" + result.stderr, encoding="utf-8"
        )
        if result.returncode != 0:
            raise RuntimeError(f"PLINK clump failed for {bac}: {result.returncode}")
        return bac, parse_clumped(prefix.with_suffix(".clumped"))

    chosen: set[tuple[str, str]] = set()
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
        futures = [executor.submit(execute, job) for job in jobs]
        for future in concurrent.futures.as_completed(futures):
            bac, rsids = future.result()
            chosen.update((bac, rsid) for rsid in rsids)
    return chosen, mapping


def resolve_decode_effect_allele(raw_ea: str, oa: str, cid: str) -> str | None:
    raw_ea, oa = raw_ea.upper(), oa.upper()
    if not raw_ea.startswith("!"):
        return raw_ea
    excluded = raw_ea[1:]
    parts = cid.split("_", 2)
    if excluded != oa or len(parts) != 3:
        return None
    alternatives = [allele.upper() for allele in parts[1:] if allele.upper() != excluded]
    return alternatives[0] if len(alternatives) == 1 else None


def outcome_classification(exp_ea: str, exp_oa: str, out_ea: str, out_oa: str) -> str:
    if is_palindromic(exp_ea, exp_oa) and {exp_ea, exp_oa} == {out_ea, out_oa}:
        return "palindromic_ambiguous"
    return allele_compatibility(exp_ea, exp_oa, out_ea, out_oa)


def et_relation() -> str:
    escaped = str(ET).replace("'", "''")
    return (
        "read_csv_auto('" + escaped
        + "', delim=' ', header=true, nullstr=['NaN','nan'], sample_size=100000)"
    )


def fetch_outcomes(selected_rows: list[dict]) -> dict[int, list[dict]]:
    con = duckdb.connect()
    con.execute("SET threads=4")
    con.execute("CREATE TEMP TABLE selected (iv_row_id INTEGER, rsID VARCHAR)")
    con.executemany(
        "INSERT INTO selected VALUES (?, ?)",
        [(row["iv_row_id"], row["analysis_rsid"]) for row in selected_rows],
    )
    freq_columns = [
        "ICE_FREQ", "DNK_FREQ", "EST_FREQ", "NOR_FREQ", "UK_FREQ",
        "USINTMT_FREQ", "USEMORY_FREQ", "LIAO_FREQ",
    ]
    select_freq = ", ".join(f"e.{column}" for column in freq_columns)
    frame = con.execute(
        f"""
        SELECT s.iv_row_id, e.CHR AS outcome_chr, e.POS AS outcome_pos_grch38,
               e.cID AS outcome_cid, e.rsID, e.OA AS outcome_oa, e.EA AS outcome_ea_raw,
               e."OR" AS outcome_or, e.SE AS se_outcome, e.P AS p_outcome,
               {select_freq}
        FROM {et_relation()} AS e
        INNER JOIN selected AS s ON e.rsID = s.rsID
        ORDER BY s.iv_row_id, e.CHR, e.POS, e.cID
        """
    ).fetchdf()
    con.close()
    matches: dict[int, list[dict]] = defaultdict(list)
    for row in frame.to_dict("records"):
        matches[int(row["iv_row_id"])].append(row)
    return matches


def harmonize_selected(rows: list[dict], matches: dict[int, list[dict]]) -> list[dict]:
    priority = {
        "exact": 1,
        "swapped": 2,
        "strand_exact": 3,
        "strand_swapped": 4,
        "palindromic_ambiguous": 5,
        "incompatible": 6,
        "decode_notation_unresolved": 7,
    }
    freq_columns = [
        "ICE_FREQ", "DNK_FREQ", "EST_FREQ", "NOR_FREQ", "UK_FREQ",
        "USINTMT_FREQ", "USEMORY_FREQ", "LIAO_FREQ",
    ]
    output = []
    for row in rows:
        options = matches.get(row["iv_row_id"], [])
        if not options:
            output.append(
                {**row, "outcome_match_count": 0, "harmonization_status": "missing_outcome", "mr_usable": "no"}
            )
            continue
        classified = []
        for option in options:
            raw_ea = str(option["outcome_ea_raw"]).upper()
            out_oa = str(option["outcome_oa"]).upper()
            out_ea = resolve_decode_effect_allele(raw_ea, out_oa, str(option["outcome_cid"]))
            status = (
                outcome_classification(row["exposure_ea"], row["exposure_oa"], out_ea, out_oa)
                if out_ea is not None
                else "decode_notation_unresolved"
            )
            classified.append((priority[status], status, option, out_ea))
        _, status, option, out_ea = min(classified, key=lambda item: item[0])
        outcome_or = float(option["outcome_or"])
        se_outcome = float(option["se_outcome"])
        p_outcome = float(option["p_outcome"])
        valid_values = outcome_or > 0 and se_outcome > 0 and 0 <= p_outcome <= 1
        usable = status in {"exact", "swapped", "strand_exact", "strand_swapped"} and valid_values
        beta_outcome = math.log(outcome_or) if usable else ""
        if usable and status in {"swapped", "strand_swapped"}:
            beta_outcome = -float(beta_outcome)
        frequencies = []
        for column in freq_columns:
            value = option[column]
            if value is not None and not (isinstance(value, float) and math.isnan(value)):
                frequencies.append(float(value) / 100.0)
        outcome_eaf_median = sorted(frequencies)[len(frequencies) // 2] if frequencies else ""
        output.append(
            {
                **row,
                "outcome_match_count": len(options),
                "outcome_chr": option["outcome_chr"],
                "outcome_pos_grch38": int(option["outcome_pos_grch38"]),
                "outcome_cid": option["outcome_cid"],
                "outcome_oa": str(option["outcome_oa"]).upper(),
                "outcome_ea_raw": str(option["outcome_ea_raw"]).upper(),
                "outcome_ea": out_ea or "",
                "outcome_eaf_median": outcome_eaf_median,
                "outcome_or": outcome_or,
                "beta_outcome_aligned": beta_outcome,
                "se_outcome": se_outcome,
                "p_outcome": p_outcome,
                "harmonization_status": status if valid_values else "invalid_outcome_values",
                "mr_usable": "yes" if usable else "no",
            }
        )
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    QA.mkdir(parents=True, exist_ok=True)
    LOG.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)
    for required in (
        MBG, ET, PANEL.with_suffix(".bed"), PANEL.with_suffix(".bim"),
        PANEL.with_suffix(".fam"), PANEL_ARCHIVE, PANEL_DOWNLOAD_META, PLINK, PLINK_ARCHIVE,
    ):
        if not required.is_file():
            raise FileNotFoundError(required)

    candidates = read_candidates()
    by_id, by_position, panel_metrics = scan_panel(candidates)
    gated_rows = panel_gate(candidates, by_id, by_position)
    eligible_unique, duplicates = unique_clump_inputs(gated_rows)
    subset_records = prepare_subset(gated_rows)
    for row in eligible_unique:
        record = subset_records[row["analysis_rsid"]]
        row.update(record)
        orientation = allele_compatibility(
            row["exposure_ea"], row["exposure_oa"], record["panel_freq_a1"], record["panel_freq_a2"]
        )
        if orientation in {"exact", "strand_exact"}:
            row["panel_exposure_eaf"] = record["panel_maf"]
        elif orientation in {"swapped", "strand_swapped"}:
            row["panel_exposure_eaf"] = 1 - record["panel_maf"]
        else:
            raise ValueError(f"Unexpected panel frequency allele mismatch: {row['analysis_rsid']}")

    chosen, taxon_mapping = run_taxon_clumping(eligible_unique)
    for row in gated_rows:
        row["clump_selected"] = (
            "yes"
            if row.get("analysis_rsid") and (row["bac"], row["analysis_rsid"]) in chosen
            else "no"
        )
    selected_rows = []
    for row in eligible_unique:
        if (row["bac"], row["analysis_rsid"]) in chosen:
            selected_rows.append({**row, "iv_row_id": len(selected_rows) + 1})
    matches = fetch_outcomes(selected_rows)
    harmonized = harmonize_selected(selected_rows, matches)

    all_by_taxon: dict[str, list[dict]] = defaultdict(list)
    gated_by_taxon: dict[str, list[dict]] = defaultdict(list)
    selected_by_taxon: dict[str, list[dict]] = defaultdict(list)
    for row in candidates:
        all_by_taxon[row["bac"]].append(row)
    for row in gated_rows:
        gated_by_taxon[row["bac"]].append(row)
    for row in harmonized:
        selected_by_taxon[row["bac"]].append(row)
    taxon_rows = []
    for bac in sorted(all_by_taxon):
        all_items = all_by_taxon[bac]
        gate_items = gated_by_taxon[bac]
        selected_items = selected_by_taxon.get(bac, [])
        usable = sum(row["mr_usable"] == "yes" for row in selected_items)
        if usable == 0:
            tier = "NO_USABLE_IV"
        elif usable == 1:
            tier = "SINGLE_IV_WALD_ONLY"
        elif usable == 2:
            tier = "TWO_IV_LIMITED"
        else:
            tier = "MULTI_IV"
        taxon_rows.append(
            {
                "bac": bac,
                "candidate_rows_p_lt_1e_5": len(all_items),
                "panel_eligible_rows": sum(row["clump_input"] == "yes" for row in gate_items),
                "clumped_independent_iv_rows": len(selected_items),
                "et_matched_independent_iv_rows": sum(row["outcome_match_count"] > 0 for row in selected_items),
                "mr_usable_independent_iv_rows": usable,
                "palindromic_ambiguous_iv_rows": sum(
                    row["harmonization_status"] == "palindromic_ambiguous" for row in selected_items
                ),
                "missing_outcome_iv_rows": sum(
                    row["harmonization_status"] == "missing_outcome" for row in selected_items
                ),
                "minimum_f_statistic": min((row["f_statistic"] for row in selected_items), default=""),
                "analysis_tier": tier,
            }
        )

    gate_counts = Counter(row["panel_gate_status"] for row in gated_rows)
    harmonization_counts = Counter(row["harmonization_status"] for row in harmonized)
    taxa_any_iv = sum(row["clumped_independent_iv_rows"] > 0 for row in taxon_rows)
    taxa_usable = sum(row["mr_usable_independent_iv_rows"] > 0 for row in taxon_rows)
    taxa_multi = sum(row["mr_usable_independent_iv_rows"] >= 3 for row in taxon_rows)
    weak_iv_rows = sum(row["f_statistic"] < 10 for row in harmonized)
    if taxa_usable == 211 and weak_iv_rows == 0:
        status_code = "PASS_LD_CLUMPING_READY_FOR_FORMAL_MR"
        status_zh = "通过LD clumping与结局协调门控；211个taxa均至少有1个可用独立工具变量，可进入正式MR"
    elif taxa_usable >= 200 and weak_iv_rows == 0:
        status_code = "PASS_PARTIAL_LD_CLUMPING_READY_FOR_FORMAL_MR"
        status_zh = "LD clumping总体通过；部分taxa无可用独立工具变量，其余可进入正式MR"
    else:
        status_code = "CAUTION_LD_CLUMPING_FEASIBILITY_REVIEW"
        status_zh = "LD clumping后工具变量保留不足，需复核阈值、代理策略或疾病路线"

    download_meta = json.loads(PANEL_DOWNLOAD_META.read_text(encoding="utf-8"))
    plink_version = subprocess.run(
        [str(PLINK), "--version"], capture_output=True, text=True, encoding="utf-8"
    ).stdout.strip()
    metrics = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "status_code": status_code,
        "status_zh": status_zh,
        "parameters": {
            "exposure_p_threshold": P_THRESHOLD,
            "clump_r2": CLUMP_R2,
            "clump_kb": CLUMP_KB,
            "ld_population": "1000 Genomes Phase 3 EUR",
            "ld_build": "GRCh37",
            "palindromic_primary_rule": "conservative exclusion because MiBioGen study EAF is unavailable",
        },
        "software": {"plink": plink_version, "duckdb": duckdb.__version__},
        "panel": {
            **panel_metrics,
            "source_url": download_meta["requested_url"],
            "archive_bytes": PANEL_ARCHIVE.stat().st_size,
            "archive_sha256": sha256(PANEL_ARCHIVE),
            "archive_etag": download_meta["etag"],
            "archive_last_modified": download_meta["last_modified"],
            "panel_bed_sha256": sha256(PANEL.with_suffix(".bed")),
            "panel_bim_sha256": sha256(PANEL.with_suffix(".bim")),
            "panel_fam_sha256": sha256(PANEL.with_suffix(".fam")),
            "source_filter_note": "IEU/OpenGWAS local LD documentation: biallelic SNPs with population MAF>0.01",
        },
        "candidate_gate": {
            "candidate_rows": len(candidates),
            "panel_gate_status_counts": dict(sorted(gate_counts.items())),
            "eligible_unique_bac_rsid_rows": len(eligible_unique),
            "deduplicated_rows": len(duplicates),
            "union_panel_rsids": len(subset_records),
        },
        "clumping": {
            "selected_bac_rsid_rows": len(harmonized),
            "selected_unique_rsids": len({row["analysis_rsid"] for row in harmonized}),
            "taxa_with_any_clumped_iv": taxa_any_iv,
            "taxa_with_any_mr_usable_iv": taxa_usable,
            "taxa_with_at_least_3_mr_usable_iv": taxa_multi,
            "weak_f_lt_10_selected_rows": weak_iv_rows,
            "harmonization_status_counts": dict(sorted(harmonization_counts.items())),
        },
        "scope_boundary_zh": "已完成工具变量独立化和结局协调，但尚未计算任何MR因果效应、FDR或敏感性检验。",
    }

    candidate_fields = [
        "candidate_row_id", "source_row_id", "bac", "exposure_chr",
        "exposure_bp_grch37_inferred", "source_rsid", "exposure_oa", "exposure_ea",
        "beta_exposure", "se_exposure", "p_exposure", "n_exposure", "f_statistic",
        "panel_chr", "panel_bp_grch37", "panel_rsid", "panel_a1", "panel_a2",
        "analysis_rsid", "panel_allele_orientation", "panel_id_match_count",
        "panel_coordinate_match_count", "panel_compatible_match_count",
        "panel_gate_status", "clump_input", "clump_selected",
    ]
    selected_fields = [
        "iv_row_id", "bac", "analysis_rsid", "source_rsid", "exposure_chr",
        "exposure_bp_grch37_inferred", "exposure_oa", "exposure_ea", "beta_exposure",
        "se_exposure", "p_exposure", "n_exposure", "f_statistic", "panel_chr",
        "panel_bp_grch37", "panel_a1", "panel_a2", "panel_maf", "panel_exposure_eaf",
        "panel_nchrobs", "outcome_match_count", "outcome_chr", "outcome_pos_grch38",
        "outcome_cid", "outcome_oa", "outcome_ea_raw", "outcome_ea",
        "outcome_eaf_median", "outcome_or", "beta_outcome_aligned", "se_outcome",
        "p_outcome", "harmonization_status", "mr_usable",
    ]
    write_csv(OUT / "stage1_ld_panel_candidate_audit.csv", gated_rows, candidate_fields)
    write_csv(OUT / "stage1_clumped_independent_ivs.csv", harmonized, selected_fields)
    write_csv(OUT / "stage1_taxon_iv_summary.csv", taxon_rows, list(taxon_rows[0].keys()))
    write_csv(
        OUT / "stage1_duplicate_resolution.csv",
        duplicates,
        ["bac", "analysis_rsid", "kept_candidate_row_id", "removed_candidate_row_id", "reason"],
    )
    write_csv(
        OUT / "stage1_taxon_file_mapping.csv",
        taxon_mapping,
        ["taxon_index", "taxon_slug", "bac", "input_rows"],
    )
    (QA / "stage1_ld_clump_metrics.json").write_text(
        json.dumps(metrics, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    tier_counts = Counter(row["analysis_tier"] for row in taxon_rows)
    report = f"""# ET MR Stage 1A：LD clumping 与结局协调审计

生成时间（UTC）：{metrics['generated_at_utc']}

## 判定报告

`{status_code}`（{status_zh}）。

## 实际执行

- 暴露候选：MiBioGen `P < 1×10⁻⁵`，共 {len(candidates):,} 行、211 taxa。
- LD 参考：1000 Genomes Phase 3 EUR、GRCh37、503 人、{panel_metrics['variants']:,} 个变异；BED/BIM/FAM 结构一致。
- 软件：{plink_version}。
- clumping：`r² < {CLUMP_R2}`、窗口 `{CLUMP_KB:,} kb`，按 taxon 独立执行。
- 面板门控：{len(eligible_unique):,} 个唯一 taxon–rsID 候选进入 clumping；门控分布为 {dict(sorted(gate_counts.items()))}。
- clumping 后：{len(harmonized):,} 个 taxon–IV 行、{len(set(row['analysis_rsid'] for row in harmonized)):,} 个唯一 rsID。
- 可用性：{taxa_any_iv}/211 taxa 至少有 1 个独立 IV；{taxa_usable}/211 taxa 至少有 1 个 ET 可用独立 IV；{taxa_multi}/211 taxa 至少有 3 个。
- 分析层级：{dict(sorted(tier_counts.items()))}。
- 工具强度：所选 IV 中 `F < 10` 为 {weak_iv_rows} 行。
- 正式协调状态：{dict(sorted(harmonization_counts.items()))}。

## 关键方法边界

MiBioGen 坐标按与 GRCh37 哨兵位点一致作推断；EUR 面板原生为 GRCh37，未执行 liftover。有效 rsID 必须同时通过坐标与等位基因门控；`.`/`NA` 仅在 GRCh37 坐标唯一且等位基因兼容时恢复为面板 rsID。

IEU/OpenGWAS 面板只保留双等位 SNP 且 EUR MAF>1%，因此“面板缺失”同时包含变异类别或频率筛选造成的缺失。回文 A/T、C/G SNP 因 MiBioGen 轻量文件没有研究 EAF，在主分析中保守排除；面板频率仅记录，不替代暴露研究 EAF。

本阶段没有计算 IVW、Wald ratio、weighted median、MR-Egger、MR-PRESSO 或 BH-FDR，也没有 risk/protective taxa。只有 `mr_usable=yes` 的独立 IV 才允许进入下一阶段正式 MR。
"""
    (ROOT / "STAGE1A_LD_CLUMP_REPORT_zh.md").write_text(report, encoding="utf-8")

    manifest_path = QA / "artifact_manifest_sha256.csv"
    manifest_rows = []
    for path in sorted(ROOT.rglob("*")):
        if not path.is_file() or path == manifest_path or path.suffix == ".pyc":
            continue
        if "work\\stage1_ld_clump\\taxa" in str(path):
            continue
        manifest_rows.append(
            {
                "relative_path": str(path.relative_to(ROOT)),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_csv(manifest_path, manifest_rows, ["relative_path", "bytes", "sha256"])
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
