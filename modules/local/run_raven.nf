
/*
* RAVEN genome-scale metabolic reconstruction.
*/


/*
* Create a draft from either sequences or annotation list.
*/

process RUN_RAVEN_DRAFT {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(faa)

    output:
    tuple val(sample_id),
          path("${sample_id}-draft.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}-reactions.tsv"),
          path("${sample_id}-kos.tsv"),
          emit: evidence

    tuple val(sample_id),
          path("${sample_id}-metadata.json"),
          emit: metadata

    tuple val(sample_id),
          path("${sample_id}.raven-draft.log"),
          emit: logs

    script:
    def artefact_arg = params.raven_artefact_dir ?
        "--artefact-dir '${params.raven_artefact_dir}'" : ""

    def library_arg = params.raven_library ?
        "--library '${params.raven_library}'" : ""

    """
    set -o pipefail

    python ${projectDir}/bin/run_raven_draft.py \
        --fasta '${faa}' \
        --method '${params.raven_method}' \
        --domain '${params.raven_domain}' \
        --threads ${task.cpus} \
        --model-id '${sample_id}' \
        --output '${sample_id}-draft.xml' \
        --evidence '${sample_id}-reactions.tsv' \
        --ko-evidence '${sample_id}-kos.tsv' \
        --metadata '${sample_id}-metadata.json' \
        ${artefact_arg} \
        ${library_arg} \
        2>&1 | tee '${sample_id}.raven-draft.log'
    """

    stub:
    """
    touch '${sample_id}-draft.xml'
    printf 'reaction_id\\n' > '${sample_id}-reactions.tsv'
    printf 'ko_id\\tgene_id\\n' > '${sample_id}-kos.tsv'
    printf '{"selected_method":"stub"}\\n' \
        > '${sample_id}-metadata.json'

    echo 'Stub RAVEN draft: ${sample_id}' \
        > '${sample_id}.raven-draft.log'
    """
}

/*
process RAVEN_ADAPT {
    // Optional reaction customization
}

process RAVEN_GAPFILL {
    // Gap-filling using defined media
}

process RAVEN_REFINE {
    // Optional post-gapfill refinement
}
*/