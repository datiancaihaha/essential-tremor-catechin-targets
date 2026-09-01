from __future__ import annotations

import csv
import html
import json
import re
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RAW_HTML = (
    PROJECT
    / "raw"
    / "stage3_literature"
    / "swisstargetprediction_job_1517142978.html"
)
OUT = PROJECT / "outputs"
OUT.mkdir(parents=True, exist_ok=True)

SMILES = "C1=CC(=C(C=C1CCCCC(=O)O)O)O"
SWISS_JOB_URL = (
    "https://www.swisstargetprediction.ch/"
    "result.php?job=1517142978&organism=Homo_sapiens"
)


def clean_cell(value: str) -> str:
    value = re.sub(r"(?i)<br\s*/?>", "; ", value)
    value = re.sub(r"<[^>]+>", "", value)
    value = html.unescape(value)
    return re.sub(r"\s+", " ", value).strip()


def main() -> None:
    source = RAW_HTML.read_text(encoding="utf-8", errors="replace")
    table_match = re.search(
        r'<table[^>]*id=["\']resultTable["\'][^>]*>(.*?)</table>',
        source,
        re.S,
    )
    if not table_match:
        raise RuntimeError("SwissTargetPrediction resultTable not found")

    records: list[dict[str, object]] = []
    for row in re.findall(r"<tr[^>]*>(.*?)</tr>", table_match.group(1), re.S):
        cells = re.findall(r"<td[^>]*>(.*?)</td>", row, re.S)
        if len(cells) < 7:
            continue
        cleaned = [clean_cell(cell) for cell in cells]
        probability_match = re.search(r"-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", cleaned[5])
        probability = float(probability_match.group(0)) if probability_match else 0.0
        actives_match = re.search(r"(\d+)\s*/\s*(\d+)", cleaned[6])
        records.append(
            {
                "rank": len(records) + 1,
                "target_name": cleaned[0],
                "gene_symbol": cleaned[1],
                "uniprot_id": cleaned[2],
                "chembl_id": cleaned[3],
                "target_class": cleaned[4],
                "probability": probability,
                "known_actives_3d": int(actives_match.group(1)) if actives_match else 0,
                "known_actives_2d": int(actives_match.group(2)) if actives_match else 0,
            }
        )

    if not records:
        raise RuntimeError("SwissTargetPrediction table contained no data rows")
    records.sort(key=lambda item: (-float(item["probability"]), int(item["rank"])))
    for rank, item in enumerate(records, start=1):
        item["rank"] = rank

    fields = [
        "rank", "target_name", "gene_symbol", "uniprot_id", "chembl_id",
        "target_class", "probability", "known_actives_3d", "known_actives_2d",
    ]
    with (OUT / "stage3_swisstargetprediction_all.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(records)

    nonzero = [row for row in records if float(row["probability"]) > 0]
    with (OUT / "stage3_swisstargetprediction_nonzero.csv").open(
        "w", newline="", encoding="utf-8-sig"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(nonzero)

    status = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "compound": "5-(3,4-Dihydroxyphenyl)pentanoic acid",
        "pubchem_cid": 49831816,
        "smiles": SMILES,
        "organism": "Homo sapiens",
        "swiss_target_prediction": {
            "status": "SUCCESS",
            "job_url": SWISS_JOB_URL,
            "all_rows": len(records),
            "probability_gt_0_rows": len(nonzero),
            "interpretation_zh": "配体相似性计算预测，未经实验验证",
        },
        "sea": {
            "status": "UNAVAILABLE_TRANSPORT_FAILURE",
            "urls_attempted": [
                "https://seadev.docking.org/",
                "https://sea.docking.org/",
                "https://sea.bkslab.org/",
            ],
            "interpretation_zh": "官方SEA服务在本次执行中发生EOF或SSL连接失败，未获得可核验结果",
        },
        "consensus_target_status": "NOT_COMPUTABLE_WITHOUT_SEA",
        "consensus_target_status_zh": "因SEA未完成，不能生成SwissTargetPrediction与SEA共同靶点",
    }
    (OUT / "stage3_target_prediction_status.json").write_text(
        json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(status, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
