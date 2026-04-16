process PURGE_HAPLOTYPIC_DUPLICATES {
    tag "Purging haplotypic duplicates from merged $sample_id assembly with Purge_dups"
    conda 'bioconda::purge_dups bioconda::minimap2'

    input:
    tuple val(sample_id), val(assembler), path(fasta)

    output:
    tuple val(sample_id), path("${sample_id}_${assembler}_major_merged_purged.fa"), emit: purged_assembly
    tuple val(sample_id), path("${sample_id}_${assembler}_purge_diagnostics.tar.gz"), emit: purge_diagnostics

    script:
    """
    set -euo pipefail

    fallback_reason=""

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
        fallback_reason="fewer_than_two_contigs"
        cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
    else
        split_fa $fasta > ${sample_id}.split.fa
        minimap2 -x asm5 -DP ${sample_id}.split.fa ${sample_id}.split.fa > ${sample_id}.self.paf

        pbcstat ${sample_id}.self.paf
        if ! awk 'BEGIN{ok=0} {if (NF>=2 && \$2+0 > 0) ok=1} END{exit(ok?0:1)}' PB.stat; then
            fallback_reason="degenerate_pb_stat"
            cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
        else
            set +e
            calcuts PB.stat > ${sample_id}.cutoffs
            calcuts_rc=\$?
            set -e

            if [[ "\$calcuts_rc" -ne 0 ]]; then
                fallback_reason="calcuts_failed_\${calcuts_rc}"
                cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
            else
                set +e
                purge_dups -2 -T ${sample_id}.cutoffs -c PB.base.cov ${sample_id}.self.paf > ${sample_id}.dups.bed
                purge_rc=\$?
                set -e
                if [[ "\$purge_rc" -ne 0 ]]; then
                    fallback_reason="purge_dups_failed_\${purge_rc}"
                    cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
                else
                    set +e
                    get_seqs -e ${sample_id}.dups.bed $fasta
                    getseqs_rc=\$?
                    set -e

                    if [[ "\$getseqs_rc" -ne 0 ]]; then
                        fallback_reason="get_seqs_failed_\${getseqs_rc}"
                        cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
                    else
                        out=""
                        for candidate in purged.fa ${fasta}.purged.fa; do
                            if [[ -f "\$candidate" ]]; then
                                out="\$candidate"
                                break
                            fi
                        done
                        if [[ -n "\$out" ]]; then
                            mv "\$out" "${sample_id}_${assembler}_major_merged_purged.fa"
                        else
                            fallback_reason="missing_purged_output"
                            cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
                        fi
                    fi
                fi
            fi
        fi
    fi

    if [[ -n "\$fallback_reason" ]]; then
        printf "fallback_reason\\t%s\\n" "\$fallback_reason" > ${sample_id}_${assembler}_purge_fallback_reason.txt
    else
        printf "fallback_reason\\tnone\\n" > ${sample_id}_${assembler}_purge_fallback_reason.txt
    fi

    tar -czf ${sample_id}_${assembler}_purge_diagnostics.tar.gz \
      PB.stat PB.base.cov PB.cov.wig \
      ${sample_id}.split.fa ${sample_id}.self.paf \
      ${sample_id}.cutoffs ${sample_id}.dups.bed \
      ${sample_id}_${assembler}_purge_fallback_reason.txt \
      2>/dev/null || true
    """
    
    stub:
    """
    cp $fasta "${sample_id}_${assembler}_major_merged_purged.fa"
    printf "fallback_reason\\tstub_mode\\n" > ${sample_id}_${assembler}_purge_fallback_reason.txt
    tar -czf ${sample_id}_${assembler}_purge_diagnostics.tar.gz ${sample_id}_${assembler}_purge_fallback_reason.txt
    """
}
