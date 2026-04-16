process MERGE_ASSEMBLIES_QUICKMERGE {
    tag "Merging $sample_id assemblies with Quickmerge"
    conda 'bioconda::quickmerge'

    input:
    tuple val(sample_id_1), val(assembler_1), val(n50_1), val(num_large_contigs_1), path(fasta_1) //this is flye
    tuple val(sample_id_2), val(assembler_2), val(n50_2), val(num_large_contigs_2), path(fasta_2) //this is nextDenovo

    output:
    tuple val(sample_id), val(assembler), path("${sample_id}_${assembler}_major_merged.fa")

    script:
    sample_id = sample_id_1
    assembler = (n50_1 > n50_2) ? assembler_1 : assembler_2
    
    hybrid = (n50_1 > n50_2) ? fasta_1 : fasta_2
    self = (n50_1 > n50_2) ? fasta_2 : fasta_1

    if (n50_1 > n50_2) {
        println("$assembler_1 assembly more contiguous. Merging $fasta_2 into $fasta_1")
    } else { //it's very unlikely that we'll have a tie.
        println("$assembler_2 assembly more contiguous. Merging $fasta_1 into $fasta_2")
    }

    """
    set -euo pipefail

    cp "$hybrid" hybrid_assembly.fasta
    cp "$self" self_assembly.fasta

    merge_wrapper.py hybrid_assembly.fasta self_assembly.fasta

    out=""
    for candidate in merged.fasta merged.fa quickmerge_merged.fasta quickmerge_merged.fa; do
        if [[ -f "\$candidate" ]]; then
            out="\$candidate"
            break
        fi
    done

    if [[ -z "\$out" ]]; then
        out="\$(python3 - <<'PY'
from pathlib import Path

exclude = {"hybrid_assembly.fasta", "self_assembly.fasta"}
candidates = [
    p for p in Path(".").iterdir()
    if p.is_file() and p.suffix.lower() in {".fa", ".fasta", ".fna"} and p.name not in exclude
]
if not candidates:
    raise SystemExit(1)
print(max(candidates, key=lambda p: p.stat().st_mtime).name)
PY
)"
    fi

    mv "\$out" "${sample_id}_${assembler}_major_merged.fa"
    """
    
    stub:
    sample_id = sample_id_1
    assembler = (n50_1 > n50_2) ? assembler_1 : assembler_2
    
    if (n50_1 > n50_2) {
        println("$assembler_1 assembly more contiguous. Merging $fasta_2 into $fasta_1")
        """
        touch "${sample_id}_${assembler}_major_merged.fa"
        """
    }
    else { //it's very unlikely that we'll have a tie.
        println("$assembler_2 assembly more contiguous. Merging $fasta_1 into $fasta_2")
        """
        touch "${sample_id}_${assembler}_major_merged.fa"
        """ 
    }
}
