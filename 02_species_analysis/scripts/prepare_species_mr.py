from __future__ import annotations

import json
import math
import re
import subprocess
import urllib.request
from pathlib import Path

import pandas as pd


PROJECT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
VERSION = PROJECT / "outputs" / "v11_document_guided_strengthening_20260826"
ROOT = VERSION / "01_species_replication"
DOWNLOADS = ROOT / "downloads"
RESULTS = ROOT / "results"
WORK = PROJECT / "work" / "species_mr_20260826"
SOURCE = VERSION / "08_source_code"
TRAIT_MAP = (
    PROJECT
    / "outputs"
    / "v5_strengthening_20260823"
    / "01_evidence_basis"
    / "dekkers_2026_gwas_catalog_trait_map.csv"
)
MANIFEST = ROOT / "species_gwas_integrity_manifest.csv"
ET_GWAS = PROJECT / "raw" / "et_gwas" / "extracted_v1" / "G250_Essential_tremor_summary"
PLINK = PROJECT / "tools" / "plink_1.9_20250819" / "plink.exe"
LD_PREFIX = PROJECT / "raw" / "ld_reference" / "1kg_v3_eur" / "EUR"
GENE_LOCATION = PROJECT / "tools" / "magma_v1.10" / "NCBI37.3.gene.loc"

THRESHOLDS = {
    "study_wide": 5.4e-11,
    "genome_wide": 5.0e-8,
    "exploratory": 5.0e-6,
}
MIBIOGEN = {
    "Faecalibacterium": {
        "taxon_id": "genus.Faecalibacterium.id.2057",
        "beta": -0.210471475388647,
    },
    "Flavonifractor": {
        "taxon_id": "genus.Flavonifractor.id.2059",
        "beta": -0.270367558373701,
    },
}
EXPLICIT_HOST_GENES = {"ABO", "FUT2", "LCT", "FOXP1"}


def ensure_inputs() -> None:
    for path in [TRAIT_MAP, MANIFEST, ET_GWAS, PLINK, GENE_LOCATION]:
        if not path.exists():
            raise FileNotFoundError(path)
    for suffix in [".bed", ".bim", ".fam"]:
        if not Path(str(LD_PREFIX) + suffix).exists():
            raise FileNotFoundError(str(LD_PREFIX) + suffix)
    for path in [RESULTS, WORK, SOURCE]:
        path.mkdir(parents=True, exist_ok=True)


def fetch_text(url: str, path: Path) -> str:
    if not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            path.write_bytes(response.read())
    return path.read_text(encoding="utf-8")


def fetch_json(url: str, path: Path) -> dict:
    if not path.exists():
        request = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(request, timeout=60) as response:
            path.write_bytes(response.read())
    return json.loads(path.read_text(encoding="utf-8"))


def trait_metadata() -> pd.DataFrame:
    integrity = pd.read_csv(MANIFEST, dtype=str)
    catalog = pd.read_csv(TRAIT_MAP, dtype=str)
    selected = catalog[catalog["accession"].isin(integrity["accession"])].copy()
    if len(selected) != len(integrity):
        raise RuntimeError("The selected GWAS accessions do not match the integrity manifest")
    selected = selected.merge(
        integrity[["accession", "local_path", "integrity_status"]],
        on="accession",
        how="left",
        validate="one_to_one",
    )
    if not selected["integrity_status"].eq("PASS").all():
        raise RuntimeError("At least one species GWAS file failed integrity checks")

    records: list[dict] = []
    for row in selected.to_dict("records"):
        accession = row["accession"]
        interval = (
            "GCST90670001-GCST90671000"
            if int(accession[-5:]) <= 71000
            else "GCST90671001-GCST90672000"
        )
        meta_url = (
            "https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
            f"{interval}/{accession}/{accession}.tsv.gz-meta.yaml"
        )
        meta_path = Path(row["local_path"]).parent / f"{accession}.tsv.gz-meta.yaml"
        meta_text = fetch_text(meta_url, meta_path)
        prevalence_match = re.search(r"Prevalence:([0-9.]+)", meta_text)
        prevalence = float(prevalence_match.group(1)) if prevalence_match else math.nan

        assembly = row.get("genome_assembly")
        ncbi_tax_id: int | None = None
        ncbi_organism = ""
        if isinstance(assembly, str) and assembly.startswith(("GCF_", "GCA_")):
            assembly_path = Path(row["local_path"]).parent / f"{assembly}.ncbi.json"
            report = fetch_json(
                f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{assembly}/dataset_report",
                assembly_path,
            )
            reports = report.get("reports", [])
            if reports:
                organism = reports[0].get("organism", {})
                ncbi_tax_id = organism.get("tax_id")
                ncbi_organism = organism.get("organism_name", "")

        taxon = row["taxon_name"]
        genus = "Faecalibacterium" if taxon.startswith("Faecalibacterium") else "Flavonifractor"
        records.append(
            {
                **row,
                "prevalence": prevalence,
                "ncbi_taxonomy_id": ncbi_tax_id,
                "ncbi_organism_name": ncbi_organism,
                "mibiogen_taxon_id": MIBIOGEN[genus]["taxon_id"],
                "mibiogen_taxon_name": genus,
                "mibiogen_rank": "genus",
                "mibiogen_primary_beta": MIBIOGEN[genus]["beta"],
                "taxonomic_mapping": "same genus; species-resolved exposure",
                "interpretation_scope": (
                    "species-level external assessment; not an exact taxon replication"
                ),
                "source_doi": "10.1038/s41588-026-02512-2",
            }
        )
    result = pd.DataFrame.from_records(records)
    result.to_csv(RESULTS / "microbiome_taxonomic_crosswalk.csv", index=False)
    return result


def extract_candidate_variants(metadata: pd.DataFrame) -> pd.DataFrame:
    output_path = RESULTS / "species_variants_p_lt_5e-6.csv"
    if output_path.exists():
        return pd.read_csv(output_path)

    usecols = [
        "chromosome",
        "base_pair_location",
        "effect_allele",
        "other_allele",
        "beta",
        "standard_error",
        "effect_allele_frequency",
        "p_value",
        "rs_id",
    ]
    extracted: list[pd.DataFrame] = []
    for row in metadata.to_dict("records"):
        accession = row["accession"]
        path = Path(row["local_path"])
        trait_parts: list[pd.DataFrame] = []
        for chunk in pd.read_csv(
            path,
            sep="\t",
            usecols=usecols,
            chunksize=500_000,
            compression="gzip",
            low_memory=False,
        ):
            chunk["p_value"] = pd.to_numeric(chunk["p_value"], errors="coerce")
            keep = chunk["p_value"].le(THRESHOLDS["exploratory"])
            keep &= chunk["rs_id"].astype(str).str.match(r"^rs\d+$")
            if keep.any():
                trait_parts.append(chunk.loc[keep].copy())
        trait = pd.concat(trait_parts, ignore_index=True) if trait_parts else pd.DataFrame(columns=usecols)
        for column in ["beta", "standard_error", "effect_allele_frequency", "p_value"]:
            trait[column] = pd.to_numeric(trait[column], errors="coerce")
        trait = trait.dropna(subset=["beta", "standard_error", "p_value", "rs_id"])
        trait = trait[trait["standard_error"].gt(0)].copy()
        trait["f_statistic"] = (trait["beta"] / trait["standard_error"]) ** 2
        trait = trait[trait["f_statistic"].gt(10)].copy()
        trait.insert(0, "taxon_name", row["taxon_name"])
        trait.insert(0, "gwas_accession", accession)
        extracted.append(trait)
        print(f"{accession}: {len(trait)} variants at P < 5e-6 and F > 10", flush=True)
    result = pd.concat(extracted, ignore_index=True)
    result.to_csv(output_path, index=False)
    return result


def read_gene_locations() -> dict[str, list[tuple[int, int, str]]]:
    table = pd.read_csv(
        GENE_LOCATION,
        sep="\t",
        header=None,
        names=["entrez", "chromosome", "start", "end", "strand", "gene_symbol"],
        dtype={"chromosome": str, "gene_symbol": str},
    )
    result: dict[str, list[tuple[int, int, str]]] = {}
    for chromosome, group in table.groupby("chromosome"):
        result[str(chromosome)] = [
            (int(r.start), int(r.end), str(r.gene_symbol))
            for r in group.itertuples(index=False)
        ]
    return result


def nearest_gene(chromosome: str, position: int, genes: dict[str, list[tuple[int, int, str]]]) -> str:
    candidates = genes.get(str(chromosome).replace("chr", ""), [])
    if not candidates:
        return ""
    best_gene = ""
    best_distance = math.inf
    for start, end, symbol in candidates:
        if start <= position <= end:
            return symbol
        distance = min(abs(position - start), abs(position - end))
        if distance < best_distance:
            best_distance = distance
            best_gene = symbol
    return best_gene


def clump_variants(variants: pd.DataFrame) -> pd.DataFrame:
    output_path = RESULTS / "species_ld_clumped_instruments.csv"
    if output_path.exists():
        return pd.read_csv(output_path)
    genes = read_gene_locations()
    output: list[pd.DataFrame] = []
    for (accession, taxon), group in variants.groupby(["gwas_accession", "taxon_name"]):
        prefix = WORK / accession
        input_path = WORK / f"{accession}.clump.tsv"
        group[["rs_id", "p_value"]].rename(
            columns={"rs_id": "SNP", "p_value": "P"}
        ).drop_duplicates("SNP").to_csv(input_path, sep="\t", index=False)
        command = [
            str(PLINK),
            "--bfile",
            str(LD_PREFIX),
            "--clump",
            str(input_path),
            "--clump-snp-field",
            "SNP",
            "--clump-field",
            "P",
            "--clump-p1",
            str(THRESHOLDS["exploratory"]),
            "--clump-p2",
            str(THRESHOLDS["exploratory"]),
            "--clump-r2",
            "0.001",
            "--clump-kb",
            "10000",
            "--out",
            str(prefix),
        ]
        completed = subprocess.run(command, capture_output=True, text=True, check=False)
        (WORK / f"{accession}.command.txt").write_text(
            subprocess.list2cmdline(command), encoding="utf-8"
        )
        (WORK / f"{accession}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (WORK / f"{accession}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if completed.returncode != 0:
            raise RuntimeError(f"PLINK clumping failed for {accession}: {completed.stderr}")
        clumped_path = Path(str(prefix) + ".clumped")
        if not clumped_path.exists() or clumped_path.stat().st_size == 0:
            print(f"{accession}: no variants present in the European LD reference", flush=True)
            continue
        clumped = pd.read_csv(clumped_path, sep=r"\s+")
        if clumped.empty:
            continue
        selected = group[group["rs_id"].isin(clumped["SNP"])].copy()
        selected = selected.sort_values("p_value").drop_duplicates("rs_id")
        selected["nearest_gene_grch37"] = [
            nearest_gene(str(chromosome), int(position), genes)
            for chromosome, position in zip(
                selected["chromosome"], selected["base_pair_location"]
            )
        ]
        selected["explicit_host_locus"] = selected["nearest_gene_grch37"].where(
            selected["nearest_gene_grch37"].isin(EXPLICIT_HOST_GENES), ""
        )
        output.append(selected)
        print(f"{accession}: {len(selected)} independent instruments at P < 5e-6", flush=True)
    result = pd.concat(output, ignore_index=True) if output else pd.DataFrame()
    if not result.empty:
        for tier, threshold in THRESHOLDS.items():
            result[f"included_{tier}"] = result["p_value"].le(threshold)
    result.to_csv(output_path, index=False)
    return result


def extract_outcome(instruments: pd.DataFrame) -> pd.DataFrame:
    output_path = RESULTS / "essential_tremor_outcome_matches.csv"
    if output_path.exists():
        return pd.read_csv(output_path)
    rsids = set(instruments["rs_id"].dropna().astype(str))
    usecols = [
        "CHR",
        "POS",
        "cID",
        "rsID",
        "OA",
        "EA",
        "ICE_FREQ",
        "DNK_FREQ",
        "EST_FREQ",
        "NOR_FREQ",
        "UK_FREQ",
        "USINTMT_FREQ",
        "USEMORY_FREQ",
        "LIAO_FREQ",
        "OR",
        "SE",
        "P",
    ]
    matches: list[pd.DataFrame] = []
    for chunk in pd.read_csv(
        ET_GWAS,
        sep=r"\s+",
        usecols=usecols,
        chunksize=500_000,
        low_memory=False,
    ):
        keep = chunk["rsID"].astype(str).isin(rsids)
        if keep.any():
            matches.append(chunk.loc[keep].copy())
    result = pd.concat(matches, ignore_index=True) if matches else pd.DataFrame(columns=usecols)
    result.to_csv(output_path, index=False)
    return result


def compatible_pair(exposure_ea: str, exposure_oa: str, outcome_ea: str, outcome_oa: str) -> bool:
    valid = set("ACGT")
    alleles = [exposure_ea, exposure_oa, outcome_ea, outcome_oa]
    if any(not isinstance(x, str) or len(x) != 1 or x not in valid for x in alleles):
        return False
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    exposure_set = {exposure_ea, exposure_oa}
    outcome_set = {outcome_ea, outcome_oa}
    complement_set = {complement[outcome_ea], complement[outcome_oa]}
    return exposure_set == outcome_set or exposure_set == complement_set


def harmonization_input(instruments: pd.DataFrame, outcome: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    frequency_columns = [
        "ICE_FREQ",
        "DNK_FREQ",
        "EST_FREQ",
        "NOR_FREQ",
        "UK_FREQ",
        "USINTMT_FREQ",
        "USEMORY_FREQ",
        "LIAO_FREQ",
    ]
    meta_by_accession = metadata.set_index("accession").to_dict("index")
    outcome_by_rsid = {key: group for key, group in outcome.groupby("rsID")}
    for instrument in instruments.to_dict("records"):
        rsid = instrument["rs_id"]
        candidates = outcome_by_rsid.get(rsid, pd.DataFrame())
        compatible: list[dict] = []
        for candidate in candidates.to_dict("records"):
            exposure_ea = str(instrument["effect_allele"]).upper()
            exposure_oa = str(instrument["other_allele"]).upper()
            outcome_ea = str(candidate["EA"]).upper()
            outcome_oa = str(candidate["OA"]).upper()
            if compatible_pair(exposure_ea, exposure_oa, outcome_ea, outcome_oa):
                compatible.append(candidate)
        compatible.sort(key=lambda x: float(x["P"]) if pd.notna(x["P"]) else math.inf)
        chosen = compatible[0] if compatible else None
        meta = meta_by_accession[instrument["gwas_accession"]]
        record = dict(instrument)
        record["outcome_match_count"] = len(candidates)
        record["allele_compatible_outcome_count"] = len(compatible)
        record["prevalence"] = meta.get("prevalence")
        record["exposure_sample_size"] = 16017
        if chosen is None:
            record["outcome_selection_status"] = (
                "missing_outcome" if len(candidates) == 0 else "no_standard_allele_match"
            )
        else:
            frequencies = pd.to_numeric(
                pd.Series([chosen.get(column) for column in frequency_columns]),
                errors="coerce",
            )
            record.update(
                {
                    "outcome_selection_status": "standard_allele_match",
                    "outcome_chr_grch38": chosen["CHR"],
                    "outcome_pos_grch38": chosen["POS"],
                    "outcome_cid": chosen["cID"],
                    "outcome_other_allele": chosen["OA"],
                    "outcome_effect_allele": chosen["EA"],
                    "outcome_eaf_median": frequencies.median(skipna=True) / 100.0,
                    "outcome_or": chosen["OR"],
                    "outcome_beta": math.log(float(chosen["OR"])),
                    "outcome_se": chosen["SE"],
                    "outcome_p": chosen["P"],
                    "outcome_case_count": 16480,
                    "outcome_control_count": 1936173,
                }
            )
        rows.append(record)
    result = pd.DataFrame.from_records(rows)
    result.to_csv(RESULTS / "species_mr_harmonization_input.csv", index=False)
    return result


def write_summary(
    metadata: pd.DataFrame,
    variants: pd.DataFrame,
    instruments: pd.DataFrame,
    outcome: pd.DataFrame,
    harmonization: pd.DataFrame,
) -> None:
    summary = {
        "species_gwas_files": int(len(metadata)),
        "species_gwas_integrity_pass": int(metadata["integrity_status"].eq("PASS").sum()),
        "variants_p_lt_5e_6_f_gt_10": int(len(variants)),
        "independent_instruments_p_lt_5e_6": int(len(instruments)),
        "independent_instruments_p_lt_5e_8": int(instruments["included_genome_wide"].sum()),
        "independent_instruments_p_lt_5_4e_11": int(instruments["included_study_wide"].sum()),
        "unique_instrument_rsids": int(instruments["rs_id"].nunique()),
        "unique_et_outcome_rsids": int(outcome["rsID"].nunique()) if not outcome.empty else 0,
        "standard_allele_matches": int(
            harmonization["outcome_selection_status"].eq("standard_allele_match").sum()
        ),
        "ld_reference": "1000 Genomes Project phase 3 European subset",
        "clumping_r2": 0.001,
        "clumping_window_kb": 10000,
        "candidate_result_threshold": "nominal P < 0.05",
        "bh_fdr_role": "supplementary evidence strength only",
    }
    (RESULTS / "species_mr_preparation_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


def main() -> None:
    ensure_inputs()
    metadata = trait_metadata()
    variants = extract_candidate_variants(metadata)
    instruments = clump_variants(variants)
    if instruments.empty:
        raise RuntimeError("No independent instruments were retained")
    outcome = extract_outcome(instruments)
    harmonization = harmonization_input(instruments, outcome, metadata)
    write_summary(metadata, variants, instruments, outcome, harmonization)
    print(f"Species MR preparation completed: {RESULTS}", flush=True)


if __name__ == "__main__":
    main()
