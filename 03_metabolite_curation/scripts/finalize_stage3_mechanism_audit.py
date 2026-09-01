from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
OUT = PROJECT / "outputs"
RAW = PROJECT / "raw" / "stage3_literature"
QA = PROJECT / "qa"
REPORT = PROJECT / "STAGE3_TARGETED_MECHANISM_AUDIT_zh.md"
OUT.mkdir(parents=True, exist_ok=True)
QA.mkdir(parents=True, exist_ok=True)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    generated_at = datetime.now(timezone.utc).isoformat()
    selected = read_csv(OUT / "stage2_nominal_selected_taxa.csv")
    faec = next(row for row in selected if row["taxon_name"] == "Faecalibacterium")
    strict_edges = read_csv(OUT / "stage2_nominal_selected_gutmgene_strict_edges.csv")
    permissive_edges = read_csv(
        OUT / "stage2_nominal_selected_gutmgene_permissive_species_edges.csv"
    )
    swiss = read_csv(OUT / "stage3_swisstargetprediction_all.csv")
    target_status = read_json(OUT / "stage3_target_prediction_status.json")
    pubchem_properties = read_json(RAW / "pubchem_CID49831816_properties.json")
    pubchem_synonyms = read_json(RAW / "pubchem_CID49831816_synonyms.json")
    compound_properties = pubchem_properties["PropertyTable"]["Properties"][0]
    synonyms = pubchem_synonyms["InformationList"]["Information"][0]["Synonym"]

    identity = {
        "generated_at_utc": generated_at,
        "compound": "5-(3,4-Dihydroxyphenyl)pentanoic acid",
        "preferred_synonym": "3,4-Dihydroxyphenylvaleric acid",
        "pubchem_cid": 49831816,
        "chebi": "CHEBI:137503",
        "hmdb": "HMDB0029233",
        "foodb": "FDB029855",
        "molecular_formula": compound_properties.get("MolecularFormula"),
        "molecular_weight": compound_properties.get("MolecularWeight"),
        "canonical_smiles": compound_properties.get("ConnectivitySMILES")
        or compound_properties.get("SMILES"),
        "inchi": compound_properties.get("InChI"),
        "inchikey": compound_properties.get("InChIKey"),
        "xlogp": compound_properties.get("XLogP"),
        "tpsa": compound_properties.get("TPSA"),
        "pubchem_synonyms": synonyms,
        "identity_boundary_zh": (
            "本化合物是无侧链羟基的开链苯基戊酸；不得与5-(3',4'-二羟基苯基)-"
            "γ-戊内酯或其水解得到的4-羟基苯基戊酸混同。"
        ),
    }
    (OUT / "stage3_metabolite_identity.json").write_text(
        json.dumps(identity, ensure_ascii=False, indent=2), encoding="utf-8"
    )

    evidence_rows: list[dict[str, object]] = [
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "taxon_to_ET",
            "entity_1": "Faecalibacterium (genus)",
            "relationship": "genetically predicted abundance associated with lower ET odds",
            "entity_2": "Essential tremor",
            "source": "Stage 2 two-sample MR",
            "identifier": "genus.Faecalibacterium.id.2057",
            "study_context": "MiBioGen exposure; deCODE ET outcome",
            "evidence_design": "IVW MR with sensitivity analyses",
            "key_result": (
                f"OR={float(faec['odds_ratio']):.3f}; 95% CI "
                f"{float(faec['or_ci_lower']):.3f}-{float(faec['or_ci_upper']):.3f}; "
                f"P={float(faec['p_value']):.4g}; BH q={float(faec['fdr_bh']):.3f}"
            ),
            "directness_to_ET": "direct disease association, nominal only",
            "causal_interpretation": "suggestive MR signal; not FDR-significant",
            "limitations": "P<0.05 selection; q>=0.05; broad instrument threshold; exposure scale is not absolute abundance",
            "use_in_decision": "retain as nominal entry candidate per user-specified gate",
            "source_url": "",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "taxon_ET_cross_validation",
            "entity_1": "Faecalibacterium (genus)",
            "relationship": "decreased in ET and correlated with fecal SCFAs",
            "entity_2": "Essential tremor",
            "source": "Huang et al., npj Parkinson's Disease, 2023",
            "identifier": "PMID:37460569; DOI:10.1038/s41531-023-00554-5",
            "study_context": "37 ET, 37 de novo PD, 35 healthy controls",
            "evidence_design": "human cross-sectional microbiome plus fecal SCFA study",
            "key_result": "lower Faecalibacterium in ET; positive association with butyrate after FDR adjustment",
            "directness_to_ET": "direct human ET evidence for taxon, not candidate metabolite",
            "causal_interpretation": "observational association only",
            "limitations": "small cohort; cross-sectional; supports a butyrate relationship, not CID 49831816",
            "use_in_decision": "supports taxon plausibility but cannot validate the nominated metabolite",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/37460569/",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "taxon_to_metabolite",
            "entity_1": "Faecalibacterium (genus)",
            "relationship": "positive Spearman correlation in Figure 5A",
            "entity_2": "5-(3,4-Dihydroxyphenyl)pentanoic acid (M6)",
            "source": "Li et al., Food Chemistry, 2023",
            "identifier": "PMID:36565551; DOI:10.1016/j.foodchem.2022.135203",
            "study_context": "M-SHIME microbiota from two selected human donors",
            "evidence_design": "quantitative microbiome-metabolite correlation",
            "key_result": "Faecalibacterium-M6 heatmap association marked significant (P<0.05)",
            "directness_to_ET": "not an ET study",
            "causal_interpretation": "correlation; does not establish direct production",
            "limitations": "two donors; community fermentation; M6 was tentatively annotated; no isolate or enzyme attribution",
            "use_in_decision": "only strict same-rank gutMGene edge; retain but downgrade mechanistic language",
            "source_url": "https://doi.org/10.1016/j.foodchem.2022.135203",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "metabolite_human_exposure",
            "entity_1": "(+)-catechin-derived gut microbial metabolism",
            "relationship": "phase-II metabolite detected in plasma",
            "entity_2": "M6-related circulating conjugate",
            "source": "Li et al., Food Chemistry, 2023",
            "identifier": "PMID:36565551",
            "study_context": "human pilot component, two preselected converter donors",
            "evidence_design": "targeted plasma metabolite detection",
            "key_result": "phase-II metabolites of M6 were detected in both donors",
            "directness_to_ET": "not an ET study; no brain measurement",
            "causal_interpretation": "supports possible systemic exposure only",
            "limitations": "n=2; conjugate rather than parent compound; no quantified brain or CSF exposure",
            "use_in_decision": "insufficient for a cerebellar mechanism claim",
            "source_url": "https://doi.org/10.1016/j.foodchem.2022.135203",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "metabolite_identity",
            "entity_1": "CID 49831816",
            "relationship": "standardized identity",
            "entity_2": "5-(3,4-Dihydroxyphenyl)pentanoic acid",
            "source": "PubChem PUG REST; gutMGene cross-identifiers",
            "identifier": "CID:49831816; CHEBI:137503; HMDB0029233",
            "study_context": "chemical database record",
            "evidence_design": "database identity mapping",
            "key_result": "C11H14O4; MW 210.23; XLogP 1.7; TPSA 77.8; exact SMILES locked",
            "directness_to_ET": "none",
            "causal_interpretation": "identity only",
            "limitations": "predicted descriptors do not prove BBB permeability or bioactivity",
            "use_in_decision": "use this exact structure for target prediction",
            "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/49831816",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "direct_ET_or_neuro_search",
            "entity_1": "CID 49831816 exact-name variants",
            "relationship": "PubMed Title/Abstract search",
            "entity_2": "ET or neurologic/brain evidence",
            "source": "PubMed E-utilities search executed 2026-08-22",
            "identifier": "stage3_literature_search_log.csv",
            "study_context": "exact term variants with ET, tremor, brain, neurologic, cerebellar, or BBB terms",
            "evidence_design": "reproducible bibliographic search",
            "key_result": "0 records for exact-compound-plus-ET and 0 for exact-compound-plus-neuro queries",
            "directness_to_ET": "no direct indexed Title/Abstract evidence found",
            "causal_interpretation": "absence of retrieved records is not proof of biological absence",
            "limitations": "does not search paywalled full text or non-indexed literature; nomenclature variants remain possible",
            "use_in_decision": "do not claim direct ET, cerebellar, or BBB evidence for the exact compound",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/",
        },
        {
            "route_scope": "ANALOG_ONLY",
            "chain_link": "brain_exposure_analog",
            "entity_1": "5-(3',4'-dihydroxyphenyl)-gamma-valerolactone",
            "relationship": "sulfated metabolite detected in brain",
            "entity_2": "rat and pig brain tissue",
            "source": "Angelino et al., Nutrients, 2019",
            "identifier": "PMID:31694297; PMCID:PMC6893823",
            "study_context": "in silico, endothelial-cell, rat, and pig models",
            "evidence_design": "multi-model BBB study",
            "key_result": "a hydroxyphenyl-gamma-valerolactone sulfate was detected in brain tissue",
            "directness_to_ET": "not ET and not the exact open-chain compound",
            "causal_interpretation": "analog evidence only",
            "limitations": "different ring form/conjugate; exact isomer uncertain; most values near detection limits",
            "use_in_decision": "may motivate testing, but must not be used as BBB proof for CID 49831816",
            "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6893823/",
        },
        {
            "route_scope": "ANALOG_ONLY",
            "chain_link": "neuroinflammation_analog",
            "entity_1": "5-(3',4'-dihydroxyphenyl)-gamma-valerolactone",
            "relationship": "with curcumin reduced LPS responses",
            "entity_2": "primary rat cortical microglia",
            "source": "Marcolin et al., Nutrients, 2025",
            "identifier": "PMID:40284180; PMCID:PMC12030566",
            "study_context": "cell culture combination treatment",
            "evidence_design": "in vitro mechanistic study",
            "key_result": "combination affected NLRP3 and NOX2/Nrf2-related readouts",
            "directness_to_ET": "not ET; different molecule; combination effect",
            "causal_interpretation": "cannot be transferred to CID 49831816",
            "limitations": "closed-ring analog, cortical rather than cerebellar cells, curcumin combination",
            "use_in_decision": "exclude from direct target evidence; cite only as analog context",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/40284180/",
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "metabolite_to_target",
            "entity_1": "CID 49831816 exact SMILES",
            "relationship": "ligand-similarity prediction",
            "entity_2": "100 human protein target rows",
            "source": "SwissTargetPrediction",
            "identifier": "job:1517142978",
            "study_context": "Homo sapiens; exact PubChem structure",
            "evidence_design": "in silico target prediction",
            "key_result": "top targets CA9, CA2, CA1, CA12, IAPP, SNCA, ALOX5; all 100 retained",
            "directness_to_ET": "computational only",
            "causal_interpretation": "hypothesis generation; no binding or functional validation",
            "limitations": "similarity-based predictions; rankings may favor well-liganded target classes",
            "use_in_decision": "retain raw list, but do not use alone for formal PPI/enrichment",
            "source_url": target_status["swiss_target_prediction"]["job_url"],
        },
        {
            "route_scope": "STRICT_PRIMARY",
            "chain_link": "target_consensus",
            "entity_1": "CID 49831816 exact SMILES",
            "relationship": "SEA prediction attempt",
            "entity_2": "human protein targets",
            "source": "SEA Search Server",
            "identifier": "RRID:SCR_023754",
            "study_context": "official development, production, and legacy endpoints",
            "evidence_design": "external web service",
            "key_result": "no result retrieved because endpoints returned EOF or SSL transport failures",
            "directness_to_ET": "not available",
            "causal_interpretation": "no Swiss-SEA consensus can be calculated",
            "limitations": "external service failure, not a negative biological result",
            "use_in_decision": "hold formal common-target/PPI step until SEA is rerun",
            "source_url": "https://seadev.docking.org/",
        },
        {
            "route_scope": "PERMISSIVE_BACKUP",
            "chain_link": "taxon_to_metabolite",
            "entity_1": "Flavonifractor plautii",
            "relationship": "directly converted catechin ring-fission intermediate",
            "entity_2": "valerolactone and 4-hydroxy phenylvaleric acid",
            "source": "Kutschera et al., Journal of Applied Microbiology, 2011",
            "identifier": "PMID:21457417",
            "study_context": "isolated human intestinal strains aK2 and DSM 6740",
            "evidence_design": "in vitro bacterial culture",
            "key_result": "direct conversion was reproduced in two F. plautii strains",
            "directness_to_ET": "not ET; species-level evidence",
            "causal_interpretation": "stronger production evidence, weaker taxonomic match to MiBioGen genus exposure",
            "limitations": "genus-to-species extrapolation; products are not CID 49831816",
            "use_in_decision": "keep as mechanistically stronger but taxonomically permissive backup",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/21457417/",
        },
        {
            "route_scope": "PERMISSIVE_BACKUP",
            "chain_link": "taxon_to_metabolite",
            "entity_1": "Methanobrevibacter smithii",
            "relationship": "plasma metabolite correlations",
            "entity_2": "three non-GABA metabolites",
            "source": "gutMGene record from rectal neuroendocrine tumor study",
            "identifier": "PMID:35265196",
            "study_context": "rectal neuroendocrine tumors",
            "evidence_design": "cross-sectional multi-omic correlation",
            "key_result": "three database correlation edges",
            "directness_to_ET": "not ET; species-level evidence",
            "causal_interpretation": "correlation only",
            "limitations": "disease-specific tumor context and genus-to-species extrapolation",
            "use_in_decision": "low-priority backup only",
            "source_url": "https://pubmed.ncbi.nlm.nih.gov/35265196/",
        },
    ]
    evidence_fields = [
        "route_scope", "chain_link", "entity_1", "relationship", "entity_2",
        "source", "identifier", "study_context", "evidence_design", "key_result",
        "directness_to_ET", "causal_interpretation", "limitations", "use_in_decision",
        "source_url",
    ]
    write_csv(OUT / "stage3_mechanism_evidence_table.csv", evidence_rows, evidence_fields)

    search_rows: list[dict[str, object]] = []
    for path in sorted(RAW.glob("pubmed_esearch_*.json")):
        payload = read_json(path)["esearchresult"]
        search_rows.append(
            {
                "source": "PubMed E-utilities",
                "tier": "T1",
                "label": path.stem.removeprefix("pubmed_esearch_"),
                "query": payload.get("querytranslation", ""),
                "reported_count": payload.get("count", ""),
                "retrieved_count": len(payload.get("idlist", [])),
                "status": "SUCCESS",
                "note_zh": "检索计数是书目检索结果，不等于生物学阴性证据",
                "generated_at_utc": generated_at,
            }
        )
    for path in sorted(RAW.glob("crossref_*.json")):
        payload = read_json(path)["message"]
        search_rows.append(
            {
                "source": "Crossref REST API",
                "tier": "T1",
                "label": path.stem.removeprefix("crossref_"),
                "query": payload.get("query", {}).get("search-terms", "")
                if isinstance(payload.get("query"), dict)
                else "",
                "reported_count": payload.get("total-results", ""),
                "retrieved_count": len(payload.get("items", [])),
                "status": "SUCCESS",
                "note_zh": "仅用于发现与元数据核验；宽泛相关性排序不能视为精确命中",
                "generated_at_utc": generated_at,
            }
        )
    for label in ("crossref_exact_metabolite_neuro", "crossref_faecalibacterium_et"):
        search_rows.append(
            {
                "source": "Crossref REST API",
                "tier": "T1",
                "label": label,
                "query": "",
                "reported_count": "",
                "retrieved_count": 0,
                "status": "TRANSPORT_FAILURE",
                "note_zh": "远端在分块响应期间断开；不得解释为零结果",
                "generated_at_utc": generated_at,
            }
        )
    search_fields = [
        "source", "tier", "label", "query", "reported_count", "retrieved_count",
        "status", "note_zh", "generated_at_utc",
    ]
    write_csv(OUT / "stage3_literature_search_log.csv", search_rows, search_fields)

    swiss_top = swiss[:15]
    top_table = "\n".join(
        f"| {row['rank']} | {row['gene_symbol']} | {row['target_name']} | {float(row['probability']):.3f} |"
        for row in swiss_top
    )

    report = f"""# ET-microbiota Stage 3: 定向机制证据与精确代谢物靶点预测审计

生成时间（UTC）：{generated_at}

## 判定报告

`CONDITIONAL_HOLD_STAGE3_NO_SEA_CONSENSUS`（条件性暂缓：严格候选与精确化学结构已锁定，SwissTargetPrediction 已完成；但严格菌-代谢物边仅为相关性，未找到精确代谢物的直接 ET/脑证据，且 SEA 外部服务失败，尚不能形成共同靶点或进入正式 PPI/富集）。

本阶段继续遵守用户指定规则：**IVW nominal P<0.05 是下游入口，BH-FDR 只报告、不作为入选门槛**。Faecalibacterium 的 IVW P={float(faec['p_value']):.4g}、BH q={float(faec['fdr_bh']):.3f}，因此它仍是“名义提示性候选”，不是“FDR 显著因果菌”。

## 决策摘要

1. **严格主线保留，但只用于假说生成。** Faecalibacterium 是唯一同分类层级的人源、非 GABA gutMGene 命中；ET 人体队列也观察到其丰度下降。然而，ET 队列关联的是丁酸等短链脂肪酸，不是本次代谢物。
2. **原始 gutMGene 边已降级为相关证据。** PMID 36565551 的 Figure 5A 显示 Faecalibacterium 与 M6 的 Spearman 显著相关，但研究使用两个预筛选供体的完整菌群 M-SHIME；作者明确指出仍需识别负责各代谢步骤的具体物种和酶。不能写成“Faecalibacterium 产生 M6”。
3. **化合物身份已锁定。** 后续统一使用 PubChem CID 49831816 的开链 `5-(3,4-Dihydroxyphenyl)pentanoic acid`，Canonical SMILES 为 `{identity['canonical_smiles']}`。它与闭环的 phenyl-gamma-valerolactone 及带 4-羟基侧链的 gamma-valeric acid 不是同一化合物。
4. **没有精确代谢物的 ET 或脑直接证据。** PubMed Title/Abstract 的精确名称加 ET/震颤检索为 0，精确名称加 brain/neurologic/cerebellar/BBB 检索也为 0。该结论只表示本次检索未找到，不表示证明不存在。
5. **脑暴露与小胶质细胞文献只能作为结构类似物背景。** PMID 31694297 和 PMID 40284180 研究的是苯基-gamma-戊内酯或其硫酸化物，不能外推成 CID 49831816 已过 BBB 或已调控 NLRP3/NOX2/Nrf2。
6. **SwissTargetPrediction 已完成，SEA 未完成。** 精确结构的人源预测得到 100 行；SEA 三个官方/历史端点均发生传输或 SSL 失败。因此没有合法的 Swiss-SEA 共同靶点，暂不运行 PPI、GO/KEGG 或 CytoHubba。

## 严格主线证据链

### 1. Faecalibacterium -> ET

- Stage 2 MR：OR {float(faec['odds_ratio']):.3f}（95% CI {float(faec['or_ci_lower']):.3f}-{float(faec['or_ci_upper']):.3f}），IVW P={float(faec['p_value']):.4g}，BH q={float(faec['fdr_bh']):.3f}；敏感性方法方向一致，未见异质性、Egger 截距或 MR-PRESSO 全局警示。
- 人体交叉验证：PMID 37460569 在 37 ET、37 初诊 PD、35 健康对照中观察到 ET 的 Faecalibacterium 降低，并在 FDR 校正后与丁酸正相关。这支持 taxon 的疾病相关性，但不验证 CID 49831816。

### 2. Faecalibacterium -> CID 49831816

- gutMGene 对应原文 PMID 36565551。该研究使用两个按儿茶素代谢速度预选的供体，建立 M-SHIME 社区发酵模型。
- M6 为 `5-(3',4'-dihydroxyphenyl)-valeric acid`，根据 MS2 碎片和既往文献作**暂定鉴定**；Figure 5A 中 Faecalibacterium 与 M6 的相关格标有 P<0.05。
- 研究的人体小型验证检测到两名供体血浆中的 M6 相关 II 相代谢物，但未测 ET、脑或脑脊液。
- 因为证据来自群落相关热图，且作者要求进一步识别物种与酶，正式措辞只能是“与 M6 相关”，不能是“产生 M6”或“通过 M6 介导 ET”。

### 3. CID 49831816 -> 脑/ET

- PubChem：C11H14O4，MW 210.23，XLogP 1.7，TPSA 77.8；这些理化描述不能替代 BBB 实测。
- HMDB 条目当前无可用组织定位；原研究只支持外周血中相关结合物的可检测性。
- 严格化合物未找到 ET、震颤、小脑、BBB 或神经系统直接研究。
- 苯基-gamma-戊内酯硫酸化物的动物脑检出，以及闭环 gamma-戊内酯在大鼠皮层小胶质细胞中的体外结果，仅为类似物证据，均已在证据表标为 `ANALOG_ONLY`。

## SwissTargetPrediction 精确结构结果

查询物种：Homo sapiens；作业：`1517142978`；原始 100 行全部保留。概率是“假定查询分子具有生物活性时，基于已知配体相似性预测该蛋白为靶点的概率”，不是结合常数或实验证据。

| rank | gene | target | probability |
|---:|---|---|---:|
{top_table}

SNCA 排名第 6 并不等于已证明该代谢物作用于 alpha-synuclein，也不能因其神经疾病知名度而优先包装。正式候选必须等待 SEA 共识、ET 疾病基因交集及后续实验可验证性筛选。

## 探索性备选路线

- **Flavonifractor plautii：** PMID 21457417 在分离株 aK2 与 DSM 6740 中直接观察到儿茶素开环中间体向 valerolactone 与 4-hydroxy phenylvaleric acid 的转化，微生物代谢证据强于 Faecalibacterium-M6 相关边；但 MiBioGen 暴露是 genus，gutMGene 是 species，且产物不是 CID 49831816。因此它是“机制更直接、分类映射更弱”的备选，不能并入严格主线。
- **Methanobrevibacter smithii：** 三条记录来自直肠神经内分泌肿瘤队列的血浆相关性，疾病情境与 ET 不匹配，优先级较低。

## 审计结论

- `PASS_EXACT_COMPOUND_IDENTITY`（通过：CID、交叉数据库编号、精确结构和同义词已锁定）。
- `PASS_ORIGINAL_EDGE_SOURCE_AUDIT`（通过：已核查 PMID 36565551 正文与 Figure 5A，确认是相关性而非单菌产物证明）。
- `PASS_SWISS_TARGET_PREDICTION`（通过：精确 SMILES 的 Homo sapiens 作业完成，100 行已解析并保留原始 HTML）。
- `DEFERRED_SEA_TRANSPORT_FAILURE`（延期：SEA 外部服务传输失败；这不是生物学阴性结果）。
- `NO_DIRECT_ET_METABOLITE_EVIDENCE_FOUND`（本次检索未找到：精确化合物的直接 ET、脑或小脑证据）。

## 阶段性结论与下一触发条件

当前不应直接进入 PPI/富集。继续条件是：

1. 使用同一精确 SMILES 成功取得 SEA 人源结果；
2. 规范化 Swiss 与 SEA 的 UniProt/基因符号后取交集；
3. 再与带来源和版本记录的 ET 疾病基因集合求交；
4. 只有得到数量可控且有脑/小脑可验证性的共同靶点，才进入 STRING、富集与动物 qPCR 候选选择。

若 SEA 反复不可用，需由用户明确同意替代方案后，才能改用另一个独立靶点预测器；不能悄悄把 Swiss 单库结果当作“双平台共识”。

## 交付文件

- `outputs/stage3_mechanism_evidence_table.csv`：逐证据边界表。
- `outputs/stage3_metabolite_identity.json`：精确化合物身份与结构边界。
- `outputs/stage3_literature_search_log.csv`：PubMed/Crossref 检索与故障记录。
- `outputs/stage3_swisstargetprediction_all.csv`：SwissTargetPrediction 100 行原始解析结果。
- `outputs/stage3_target_prediction_status.json`：Swiss/SEA 执行状态。
- `qa/stage3_qa.json` 与 `qa/stage3_artifact_manifest_sha256.csv`：质量核查与哈希清单。
"""
    REPORT.write_text(report, encoding="utf-8")

    checks = {
        "faecalibacterium_nominal_p_lt_0_05": float(faec["p_value"]) < 0.05,
        "faecalibacterium_bh_q_ge_0_05": float(faec["fdr_bh"]) >= 0.05,
        "strict_edge_count_is_1": len(strict_edges) == 1,
        "strict_edge_is_correlative": strict_edges[0]["Associative mode"] == "correlatively",
        "strict_edge_pmid_correct": strict_edges[0]["PMID"] == "36565551",
        "permissive_rows_present": len(permissive_edges) == 22,
        "pubchem_cid_correct": int(compound_properties["CID"]) == 49831816,
        "canonical_smiles_locked": identity["canonical_smiles"] == "C1=CC(=C(C=C1CCCCC(=O)O)O)O",
        "swiss_raw_rows_100": len(swiss) == 100,
        "swiss_top_gene_ca9": swiss[0]["gene_symbol"] == "CA9",
        "swiss_rank6_gene_snca": swiss[5]["gene_symbol"] == "SNCA",
        "sea_not_misreported_as_negative": target_status["sea"]["status"] == "UNAVAILABLE_TRANSPORT_FAILURE",
        "exact_metabolite_et_pubmed_count_zero": read_json(RAW / "pubmed_esearch_exact_metabolite_et.json")["esearchresult"]["count"] == "0",
        "exact_metabolite_neuro_pubmed_count_zero": read_json(RAW / "pubmed_esearch_exact_metabolite_neuro.json")["esearchresult"]["count"] == "0",
        "corrected_source_pmid_count_one": read_json(RAW / "pubmed_esearch_strict_source_pmid_corrected.json")["esearchresult"]["count"] == "1",
        "source_pdf_nonempty": (RAW / "Li_2023_FoodChem_PMID36565551.pdf").stat().st_size > 100_000,
    }
    qa_status = "PASS_STAGE3_QA" if all(checks.values()) else "FAIL_STAGE3_QA"
    qa_payload = {
        "generated_at_utc": generated_at,
        "status": qa_status,
        "status_zh": "Stage 3 文件与关键证据边界核查通过" if qa_status.startswith("PASS") else "Stage 3 质量核查未通过",
        "checks": checks,
        "counts": {
            "evidence_rows": len(evidence_rows),
            "search_log_rows": len(search_rows),
            "swiss_target_rows": len(swiss),
            "strict_gutmgene_edges": len(strict_edges),
            "permissive_gutmgene_edges": len(permissive_edges),
        },
    }
    qa_path = QA / "stage3_qa.json"
    qa_path.write_text(json.dumps(qa_payload, ensure_ascii=False, indent=2), encoding="utf-8")

    manifest_targets = [
        PROJECT / "README_ET_STAGE0_zh.md",
        REPORT,
        OUT / "stage3_mechanism_evidence_table.csv",
        OUT / "stage3_metabolite_identity.json",
        OUT / "stage3_literature_search_log.csv",
        OUT / "stage3_swisstargetprediction_all.csv",
        OUT / "stage3_swisstargetprediction_nonzero.csv",
        OUT / "stage3_target_prediction_status.json",
        RAW / "Li_2023_FoodChem_PMID36565551.pdf",
        RAW / "swisstargetprediction_job_1517142978.html",
        RAW / "pubmed_efetch_PMID36565551.xml",
        RAW / "pubchem_CID49831816_properties.json",
        RAW / "pubchem_CID49831816_synonyms.json",
        PROJECT / "scripts" / "run_stage3_literature_audit.py",
        PROJECT / "scripts" / "parse_stage3_target_prediction.py",
        PROJECT / "scripts" / "finalize_stage3_mechanism_audit.py",
        qa_path,
    ]
    manifest_rows = []
    for path in manifest_targets:
        manifest_rows.append(
            {
                "relative_path": path.relative_to(PROJECT).as_posix(),
                "size_bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    write_csv(
        QA / "stage3_artifact_manifest_sha256.csv",
        manifest_rows,
        ["relative_path", "size_bytes", "sha256"],
    )
    print(json.dumps(qa_payload, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
