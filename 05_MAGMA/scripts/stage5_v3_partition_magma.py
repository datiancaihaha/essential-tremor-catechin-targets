from __future__ import annotations

from pathlib import Path


ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
WORK = ROOT / "work" / "v3_parallel_evidence_20260822" / "magma"


def clean_duplicate_pvalues() -> None:
    source = WORK / "et_2024_gwas_magma_pval.tsv"
    duplicate_log = WORK / "predicted105_35up10down_Ntotal.log.suppl"
    output = WORK / "et_2024_gwas_magma_pval_deduplicated.tsv"

    duplicate_ids = {
        line.strip()
        for line in duplicate_log.read_text(encoding="utf-8").splitlines()
        if line.startswith("rs")
    }
    minimum_p: dict[str, float] = {}
    with source.open("r", encoding="utf-8") as src, output.open("w", encoding="utf-8") as dst:
        dst.write(src.readline())
        for line in src:
            rsid, p_text = line.rstrip("\r\n").split("\t")
            if rsid in duplicate_ids:
                p_value = float(p_text)
                minimum_p[rsid] = min(minimum_p.get(rsid, p_value), p_value)
            else:
                dst.write(line)
        for rsid in sorted(minimum_p):
            dst.write(f"{rsid}\t{minimum_p[rsid]:.16g}\n")
    if duplicate_ids != set(minimum_p):
        raise RuntimeError("Duplicate SNP set did not match the MAGMA supplementary log")


def partition_annotation(partitions: int = 4) -> None:
    source = WORK / "et_35up10down.genes.annot.txt"
    headers: list[str] = []
    records: list[tuple[int, str]] = []
    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.split()
            # Greedy balancing uses the number of annotated SNPs as a direct workload proxy.
            records.append((max(1, len(fields) - 2), line))

    bins: list[list[str]] = [[] for _ in range(partitions)]
    loads = [0 for _ in range(partitions)]
    for weight, line in sorted(records, key=lambda item: item[0], reverse=True):
        index = min(range(partitions), key=loads.__getitem__)
        bins[index].append(line)
        loads[index] += weight

    for index, lines in enumerate(bins, start=1):
        path = WORK / f"et_35up10down.part{index}.genes.annot.txt"
        with path.open("w", encoding="utf-8") as handle:
            handle.writelines(headers)
            handle.writelines(lines)

    summary = [f"part\tgenes\tannotated_snp_entries"]
    for index, (lines, load) in enumerate(zip(bins, loads, strict=True), start=1):
        summary.append(f"{index}\t{len(lines)}\t{load}")
    (WORK / "magma_partition_summary.tsv").write_text("\n".join(summary) + "\n", encoding="utf-8")


def main() -> None:
    clean_duplicate_pvalues()
    partition_annotation()
    print((WORK / "magma_partition_summary.tsv").read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
