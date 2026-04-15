process PURGE_HAPLOTYPIC_DUPLICATES {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    publishDir "results/merged_assemblies", mode: 'copy'
    conda 'bioconda::purge_dups bioconda::minimap2'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_major_merged_purged.fa")

    script:
    """
    set -euo pipefail

    seq_count=\$(python3 - <<'PY'
from pathlib import Path
f = Path("$fasta")
count = 0
with f.open() as handle:
    for line in handle:
        if line.startswith(">"):
            count += 1
print(count)
PY
)

    if [[ "\$seq_count" -lt 2 ]]; then
        cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
        exit 0
    fi

    split_fa $fasta > ${sample_id}.split.fa
    minimap2 -x asm5 -DP ${sample_id}.split.fa ${sample_id}.split.fa > ${sample_id}.self.paf

    pbcstat ${sample_id}.self.paf
    calcuts PB.stat > ${sample_id}.cutoffs
    purge_dups -2 -T ${sample_id}.cutoffs -c PB.base.cov ${sample_id}.self.paf > ${sample_id}.dups.bed

    get_seqs -e ${sample_id}.dups.bed $fasta

    out=""
    for candidate in purged.fa ${fasta}.purged.fa; do
        if [[ -f "\$candidate" ]]; then
            out="\$candidate"
            break
        fi
    done

    if [[ -z "\$out" ]]; then
        echo "ERROR: purge_dups did not produce an expected output FASTA." >&2
        exit 1
    fi

    mv "\$out" "${sample_id}_${assembler}_major_merged_purged.fa"
    """
    
    stub:
    """
    cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
    """
}
