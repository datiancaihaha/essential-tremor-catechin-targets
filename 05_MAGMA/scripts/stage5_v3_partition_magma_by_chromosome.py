from __future__ import annotations

from collections import defaultdict
from pathlib import Path


ROOT = Path(r"D:\CodexProjects\ET_MR_Stage0_20260821_v1")
WORK = ROOT / "work" / "v3_parallel_evidence_20260822" / "magma"


def main() -> None:
    source = WORK / "et_35up10down.genes.annot.txt"
    headers: list[str] = []
    records_by_chr: dict[int, list[str]] = defaultdict(list)
    loads_by_chr: dict[int, int] = defaultdict(int)

    with source.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                headers.append(line)
                continue
            fields = line.split()
            chromosome = int(fields[1].split(":", 1)[0])
            if chromosome > 22:
                continue
            records_by_chr[chromosome].append(line)
            loads_by_chr[chromosome] += max(1, len(fields) - 2)

    bins: list[list[int]] = [[] for _ in range(4)]
    loads = [0, 0, 0, 0]
    for chromosome in sorted(loads_by_chr, key=loads_by_chr.get, reverse=True):
        index = min(range(4), key=loads.__getitem__)
        bins[index].append(chromosome)
        loads[index] += loads_by_chr[chromosome]

    summary = ["part\tchromosomes\tgenes\tannotated_snp_entries"]
    for index, chromosomes in enumerate(bins, start=1):
        chromosomes.sort()
        lines = [line for chromosome in chromosomes for line in records_by_chr[chromosome]]
        path = WORK / f"et_35up10down.chrpart{index}.genes.annot.txt"
        with path.open("w", encoding="utf-8") as handle:
            handle.writelines(headers)
            handle.writelines(lines)
        summary.append(
            f"{index}\t{','.join(map(str, chromosomes))}\t{len(lines)}\t{loads[index - 1]}"
        )

    output = WORK / "magma_chromosome_partition_summary.tsv"
    output.write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(output.read_text(encoding="utf-8"))


if __name__ == "__main__":
    main()
