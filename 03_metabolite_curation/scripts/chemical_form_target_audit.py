from __future__ import annotations

import csv
import json
import re
import time
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import quote

import requests


PROJECT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
VERSION = PROJECT / "outputs" / "v11_document_guided_strengthening_20260826"
ROOT = VERSION / "03_chemical_forms"
RAW = ROOT / "database_responses"
SOURCE = VERSION / "08_source_code"
INPUT = (
    PROJECT
    / "outputs"
    / "v5_strengthening_20260823"
    / "03_chemical_identity_exposure"
    / "chemical_identity_exposure_audit_data.json"
)


EXPOSURE_CLASSIFICATION = {
    "PROJECT_ACID": {
        "exact_entity_grade": "C",
        "family_grade": "C",
        "human_matrix": "plasma derivatives only",
        "identity_basis": "parent compound was not directly quantified in human plasma",
        "reference": "Li et al., Food Chemistry 2023; doi:10.1016/j.foodchem.2022.135203",
    },
    "LACTONE_UNSPECIFIED": {
        "exact_entity_grade": "B",
        "family_grade": "B",
        "human_matrix": "urine",
        "identity_basis": "aglycone detected against a reference standard in a human intervention study",
        "reference": "Hollands et al., Molecular Nutrition & Food Research 2020; doi:10.1002/mnfr.201901135",
    },
    "LACTONE_5R": {
        "exact_entity_grade": "B",
        "family_grade": "B",
        "human_matrix": "urine",
        "identity_basis": "human urinary aglycone evidence; circulating plasma forms were predominantly conjugated",
        "reference": "Hollands et al., Molecular Nutrition & Food Research 2020; doi:10.1002/mnfr.201901135",
    },
    "LACTONE_HYDROLYSIS_ACID": {
        "exact_entity_grade": "C",
        "family_grade": "C",
        "human_matrix": "plasma derivatives only",
        "identity_basis": "related phase-II derivatives were reported; the exact parent acid was not directly quantified",
        "reference": "Li et al., Food Chemistry 2023; doi:10.1016/j.foodchem.2022.135203",
    },
    "LACTONE_3_SULFATE": {
        "exact_entity_grade": "B",
        "family_grade": "A",
        "human_matrix": "plasma and urine",
        "identity_basis": "the circulating sulfate family was quantified, but positional isomers coeluted",
        "reference": "Angelino et al., American Journal of Clinical Nutrition 2023; doi:10.1016/j.ajcnut.2023.06.006",
    },
    "LACTONE_4_SULFATE": {
        "exact_entity_grade": "B",
        "family_grade": "A",
        "human_matrix": "plasma and urine",
        "identity_basis": "the circulating sulfate family was quantified, but positional isomers coeluted",
        "reference": "Angelino et al., American Journal of Clinical Nutrition 2023; doi:10.1016/j.ajcnut.2023.06.006",
    },
    "LACTONE_3_GLUCURONIDE": {
        "exact_entity_grade": "B",
        "family_grade": "A",
        "human_matrix": "plasma and urine",
        "identity_basis": "the circulating glucuronide family was detected; one positional isomer was quantified using an isomeric standard",
        "reference": "Angelino et al., American Journal of Clinical Nutrition 2023; doi:10.1016/j.ajcnut.2023.06.006",
    },
    "LACTONE_4_GLUCURONIDE": {
        "exact_entity_grade": "B",
        "family_grade": "A",
        "human_matrix": "plasma and urine",
        "identity_basis": "the circulating glucuronide family was detected; positional assignment was not supported by a matched standard for every isomer",
        "reference": "Angelino et al., American Journal of Clinical Nutrition 2023; doi:10.1016/j.ajcnut.2023.06.006",
    },
}


class ResultTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.in_table = False
        self.in_row = False
        self.in_cell = False
        self.current_cell: list[str] = []
        self.current_row: list[str] = []
        self.rows: list[list[str]] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attributes = dict(attrs)
        if tag == "table" and attributes.get("id") == "resultTable":
            self.in_table = True
        elif self.in_table and tag == "tr":
            self.in_row = True
            self.current_row = []
        elif self.in_row and tag in {"td", "th"}:
            self.in_cell = True
            self.current_cell = []

    def handle_data(self, data: str) -> None:
        if self.in_cell:
            self.current_cell.append(data)

    def handle_endtag(self, tag: str) -> None:
        if self.in_cell and tag in {"td", "th"}:
            value = " ".join(" ".join(self.current_cell).split())
            self.current_row.append(value)
            self.in_cell = False
        elif self.in_row and tag == "tr":
            if self.current_row:
                self.rows.append(self.current_row)
            self.in_row = False
        elif self.in_table and tag == "table":
            self.in_table = False


def write_csv(path: Path, rows: list[dict], fields: list[str] | None = None) -> None:
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def request_with_retries(
    session: requests.Session,
    method: str,
    url: str,
    *,
    attempts: int = 4,
    **kwargs,
) -> requests.Response:
    last_error: Exception | None = None
    for attempt in range(1, attempts + 1):
        try:
            response = session.request(method, url, **kwargs)
            if response.status_code == 429 or response.status_code >= 500:
                response.raise_for_status()
            return response
        except requests.RequestException as error:
            last_error = error
            if attempt == attempts:
                break
            time.sleep(5 * attempt)
    raise RuntimeError(f"Request failed after {attempts} attempts: {url}") from last_error


def query_json(session: requests.Session, url: str, path: Path) -> tuple[int, dict | None]:
    if path.exists() and path.stat().st_size:
        try:
            return 200, json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            pass
    response = request_with_retries(session, "GET", url, timeout=90)
    path.write_bytes(response.content)
    if response.status_code == 404:
        return 404, None
    response.raise_for_status()
    return response.status_code, response.json()


def post_json(
    session: requests.Session,
    url: str,
    payload: dict,
    path: Path,
) -> tuple[int, dict]:
    if path.exists() and path.stat().st_size:
        try:
            return 200, json.loads(path.read_text(encoding="utf-8"))
        except (json.JSONDecodeError, UnicodeDecodeError):
            pass
    response = request_with_retries(
        session,
        "POST",
        url,
        json=payload,
        headers={"Content-Type": "application/json"},
        timeout=90,
    )
    response.raise_for_status()
    path.write_bytes(response.content)
    return response.status_code, response.json()


def submit_swiss_prediction(
    session: requests.Session, record: dict, prior_job: str | None = None
) -> tuple[str, list[dict]]:
    record_id = record["record_id"]
    submission_path = RAW / f"{record_id}_swiss_submission.html"
    result_path = RAW / f"{record_id}_swiss_result.html"
    if not prior_job and submission_path.exists():
        match = re.search(
            r"result\.php\?job=(\d+)&organism=Homo_sapiens",
            submission_path.read_text(encoding="utf-8", errors="ignore"),
        )
        if match:
            prior_job = match.group(1)
    if prior_job:
        job = prior_job
    else:
        request_with_retries(
            session, "GET", "https://www.swisstargetprediction.ch/", timeout=60
        )
        response = request_with_retries(
            session,
            "POST",
            "https://www.swisstargetprediction.ch/predict.php",
            data={
                "organism": "Homo_sapiens",
                "smiles": record["isomeric_smiles"],
                "ioi": "2",
            },
            headers={
                "Origin": "https://www.swisstargetprediction.ch",
                "Referer": "https://www.swisstargetprediction.ch/",
            },
            timeout=180,
        )
        response.raise_for_status()
        submission_path.write_bytes(response.content)
        match = re.search(r"result\.php\?job=(\d+)&organism=Homo_sapiens", response.text)
        if not match:
            raise RuntimeError(f"SwissTargetPrediction did not return a result job for {record_id}")
        job = match.group(1)
    result_url = (
        f"https://www.swisstargetprediction.ch/result.php?job={job}&organism=Homo_sapiens"
    )
    if result_path.exists() and result_path.stat().st_size:
        result_text = result_path.read_text(encoding="utf-8", errors="ignore")
    else:
        response = request_with_retries(session, "GET", result_url, timeout=120)
        response.raise_for_status()
        result_path.write_bytes(response.content)
        result_text = response.text
    parser = ResultTableParser()
    parser.feed(result_text)
    if len(parser.rows) < 2:
        raise RuntimeError(f"No SwissTargetPrediction result table was found for {record_id}")
    header = parser.rows[0]
    rows: list[dict] = []
    for rank, values in enumerate(parser.rows[1:], start=1):
        if len(values) < 7:
            continue
        item = dict(zip(header, values))
        probability_match = re.search(r"[0-9]+(?:\.[0-9]+)?(?:[eE][-+]?\d+)?", item["Probability*"])
        rows.append(
            {
                "record_id": record_id,
                "standard_name": record["standard_name"],
                "pubchem_cid": record["pubchem_cid"],
                "rank": rank,
                "target_name": item["Target"],
                "gene_symbol": item["Common name"],
                "uniprot_id": item["Uniprot ID"],
                "chembl_target_id": item["ChEMBL ID"],
                "target_class": item["Target Class"],
                "probability": float(probability_match.group(0)) if probability_match else None,
                "known_actives_3d_2d": item["Known actives (3D/2D)"],
                "swiss_job": job,
                "source_url": result_url,
            }
        )
    return job, rows


def query_experimental_databases(session: requests.Session, record: dict) -> tuple[dict, list[dict]]:
    record_id = record["record_id"]
    inchikey = record["inchikey"]
    cid = record["pubchem_cid"]
    smiles = record["isomeric_smiles"]

    unichem_url = "https://www.ebi.ac.uk/unichem/api/v1/compounds"
    _, unichem = post_json(
        session,
        unichem_url,
        {"compound": inchikey, "type": "inchikey"},
        RAW / f"{record_id}_unichem_exact_mapping.json",
    )
    chembl_ids = sorted(
        {
            source["compoundId"]
            for compound in unichem.get("compounds", [])
            for source in compound.get("sources", [])
            if str(source.get("shortName", "")).lower() == "chembl"
        }
    )
    activities: list[dict] = []
    for molecule_id in chembl_ids:
        activity_url = (
            "https://www.ebi.ac.uk/chembl/api/data/activity.json"
            f"?molecule_chembl_id={quote(str(molecule_id))}&limit=1000"
        )
        _, activity_json = query_json(
            session,
            activity_url,
            RAW / f"{record_id}_{molecule_id}_chembl_activities.json",
        )
        for activity in (activity_json or {}).get("activities", []):
            activities.append(
                {
                    "record_id": record_id,
                    "standard_name": record["standard_name"],
                    "molecule_chembl_id": molecule_id,
                    "target_chembl_id": activity.get("target_chembl_id"),
                    "target_organism": activity.get("target_organism"),
                    "target_pref_name": activity.get("target_pref_name"),
                    "assay_type": activity.get("assay_type"),
                    "standard_type": activity.get("standard_type"),
                    "standard_relation": activity.get("standard_relation"),
                    "standard_value": activity.get("standard_value"),
                    "standard_units": activity.get("standard_units"),
                    "pchembl_value": activity.get("pchembl_value"),
                    "document_chembl_id": activity.get("document_chembl_id"),
                }
            )

    pubchem_url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    )
    pubchem_status, pubchem = query_json(
        session, pubchem_url, RAW / f"{record_id}_pubchem_assaysummary.json"
    )
    pubchem_rows = []
    if pubchem:
        pubchem_rows = pubchem.get("Table", {}).get("Row", [])

    bindingdb_url = (
        "https://bindingdb.org/rest/getTargetByCompound"
        f"?smiles={quote(smiles, safe='')}&cutoff=1&response=application/json"
    )
    _, bindingdb = query_json(
        session, bindingdb_url, RAW / f"{record_id}_bindingdb_exact.json"
    )
    binding_payload = (bindingdb or {}).get("getLindsByUniprotResponse", {})
    binding_affinities = binding_payload.get("bdb.affinities", []) or []

    summary = {
        "record_id": record_id,
        "standard_name": record["standard_name"],
        "inchikey": inchikey,
        "pubchem_cid": cid,
        "chembl_exact_molecule_count": len(chembl_ids),
        "chembl_exact_molecule_ids": ";".join(chembl_ids),
        "chembl_activity_count": len(activities),
        "pubchem_bioassay_count": len(pubchem_rows),
        "pubchem_query_status": pubchem_status,
        "bindingdb_exact_hit_count": int(binding_payload.get("bdb.hit", 0) or 0),
        "bindingdb_affinity_count": len(binding_affinities),
        "unichem_query_url": unichem_url,
        "pubchem_query_url": pubchem_url,
        "bindingdb_query_url": bindingdb_url,
    }
    return summary, activities


def main() -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    RAW.mkdir(parents=True, exist_ok=True)
    SOURCE.mkdir(parents=True, exist_ok=True)
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    compounds = source["chemical_identity"]
    missing = [record["record_id"] for record in compounds if record["record_id"] not in EXPOSURE_CLASSIFICATION]
    if missing:
        raise RuntimeError(f"Exposure classification is missing for: {missing}")

    session = requests.Session()
    session.headers.update(
        {
            "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/140.0 Safari/537.36",
            "Accept": "*/*",
        }
    )
    exposure_rows: list[dict] = []
    prediction_rows: list[dict] = []
    database_rows: list[dict] = []
    activity_rows: list[dict] = []
    jobs: dict[str, str] = {}

    for index, record in enumerate(compounds):
        record_id = record["record_id"]
        exposure_rows.append({**record, **EXPOSURE_CLASSIFICATION[record_id]})
        job, predictions = submit_swiss_prediction(session, record)
        jobs[record_id] = job
        prediction_rows.extend(predictions)
        database_summary, activities = query_experimental_databases(session, record)
        database_rows.append(database_summary)
        activity_rows.extend(activities)
        print(
            f"{record_id}: {len(predictions)} predicted targets; "
            f"{database_summary['chembl_activity_count']} ChEMBL activities; "
            f"{database_summary['pubchem_bioassay_count']} PubChem assays; "
            f"{database_summary['bindingdb_affinity_count']} BindingDB affinities",
            flush=True,
        )
        if index < len(compounds) - 1:
            time.sleep(5)

        write_csv(ROOT / "chemical_exposure_classification.csv", exposure_rows)
        write_csv(ROOT / "swisstargetprediction_all_chemical_forms.csv", prediction_rows)
        write_csv(ROOT / "experimental_target_database_summary.csv", database_rows)
        write_csv(
            ROOT / "chembl_experimental_activities.csv",
            activity_rows,
            fields=[
                "record_id",
                "standard_name",
                "molecule_chembl_id",
                "target_chembl_id",
                "target_organism",
                "target_pref_name",
                "assay_type",
                "standard_type",
                "standard_relation",
                "standard_value",
                "standard_units",
                "pchembl_value",
                "document_chembl_id",
            ],
        )
        (ROOT / "swisstargetprediction_jobs.json").write_text(
            json.dumps(jobs, indent=2), encoding="utf-8"
        )

    write_csv(ROOT / "chemical_exposure_classification.csv", exposure_rows)
    write_csv(ROOT / "swisstargetprediction_all_chemical_forms.csv", prediction_rows)
    write_csv(ROOT / "experimental_target_database_summary.csv", database_rows)
    write_csv(
        ROOT / "chembl_experimental_activities.csv",
        activity_rows,
        fields=[
            "record_id",
            "standard_name",
            "molecule_chembl_id",
            "target_chembl_id",
            "target_organism",
            "target_pref_name",
            "assay_type",
            "standard_type",
            "standard_relation",
            "standard_value",
            "standard_units",
            "pchembl_value",
            "document_chembl_id",
        ],
    )
    (ROOT / "swisstargetprediction_jobs.json").write_text(
        json.dumps(jobs, indent=2), encoding="utf-8"
    )

    grade_a_family = {
        row["record_id"]
        for row in exposure_rows
        if row["family_grade"] == "A"
    }
    grade_a_predictions = [
        row for row in prediction_rows if row["record_id"] in grade_a_family
    ]
    write_csv(ROOT / "circulating_conjugate_family_predictions.csv", grade_a_predictions)

    methylated_scope = {
        "status": "not_submitted_for_structure_specific_prediction",
        "reason": (
            "Human studies report methylated phenyl-gamma-valerolactone families, "
            "but the current evidence package does not establish a single positional isomer "
            "with an authoritative structure suitable for structure-specific prediction."
        ),
        "interpretation": (
            "Methylated forms remain a family-level exposure observation and are not assigned "
            "structure-specific targets."
        ),
        "reference": (
            "Angelino et al., American Journal of Clinical Nutrition 2023; "
            "standardized phenyl-gamma-valerolactone nomenclature recommendations"
        ),
    }
    (ROOT / "methylated_form_scope.json").write_text(
        json.dumps(methylated_scope, indent=2), encoding="utf-8"
    )
    print(f"Chemical-form target audit completed: {ROOT}", flush=True)


if __name__ == "__main__":
    main()
