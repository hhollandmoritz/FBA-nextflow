/*
 * Run gapseq reconstruction as separate stages.
 *
 * Either "DOALL" or "FIND + DRAFT + ADAPT + MEDIUM + FILL" can be used.
 * 
 * The latter option allows more customization of each step.
 *
 * In the second option, 
 * we have set `gapseq find -l all` so that all pathway databases are used:
 * MetaCyc, KEGG, SEED, and (if provided) a gapseq custom pathway table.
 *
 */


process RUN_GAPSEQ_DOALL {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(faa)

    output:
    tuple val(sample_id),
          path("${sample_id}.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}.RDS"),
          emit: rds_models

    tuple val(sample_id),
          path("${sample_id}-all-Reactions.tbl"),
          path("${sample_id}-all-Pathways.tbl"),
          path("${sample_id}-Transporter.tbl"),
          emit: evidence

    tuple val(sample_id),
          path("${sample_id}.gapseq.log"),
          emit: logs

    script:
    def extra_args = params.gapseq_args?.toString() ?: ""

    def database_arg = params.gapseq_db
        ? "-D '${params.gapseq_db}' -O"
        : ""

    def medium_arg = params.gapseq_medium
        ? "-m '${params.gapseq_medium}'"
        : ""

    """
    set -o pipefail

    gapseq doall \
        -A '${params.gapseq_aligner}' \
        -K ${task.cpus} \
        -t '${params.gapseq_taxonomy}' \
        -f . \
        ${database_arg} \
        ${medium_arg} \
        ${extra_args} \
        '${faa}' \
        2>&1 | tee '${sample_id}.gapseq.log'
    """

    stub:
    """
    touch '${sample_id}.xml'
    touch '${sample_id}.RDS'

    touch '${sample_id}-all-Reactions.tbl'
    touch '${sample_id}-all-Pathways.tbl'
    touch '${sample_id}-Transporter.tbl'

    echo 'Stub gapseq run: ${sample_id}' \
        > '${sample_id}.gapseq.log'
    """
}

process RUN_GAPSEQ_FIND {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(faa)

    output:
    tuple val(sample_id),
          path("${sample_id}-all-Reactions.tbl"),
          path("${sample_id}-all-Pathways.tbl"),
          emit: evidence

    tuple val(sample_id),
          path("${sample_id}.gapseq-find.log"),
          emit: logs

    script:
    /*
     * -D and -O are independent in gapseq doall.
     *
     * gapseq_db:
     *     optional path to the gapseq reference sequence database
     *
     * gapseq_offline:
     *     optionally force offline operation
     */
    def database_arg = params.gapseq_db
        ? "-D '${params.gapseq_db}'"
        : ""

    def offline_arg = params.gapseq_offline
        ? "-O"
        : ""

    """
    set -o pipefail

    gapseq find \
        -v 0 \
        -b 200 \
        -p all \
        -l all \
        -t '${params.gapseq_taxonomy}' \
        -K ${task.cpus} \
        -A '${params.gapseq_aligner}' \
        -f . \
        ${database_arg} \
        ${offline_arg} \
        '${faa}' \
        2>&1 | tee '${sample_id}.gapseq-find.log'
    """

    stub:
    """
    touch '${sample_id}-all-Reactions.tbl'
    touch '${sample_id}-all-Pathways.tbl'

    echo 'Stub gapseq find: ${sample_id}' \
        > '${sample_id}.gapseq-find.log'
    """
}


process RUN_GAPSEQ_TRANSPORT {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(faa)

    output:
    tuple val(sample_id),
          path("${sample_id}-Transporter.tbl"),
          emit: evidence

    tuple val(sample_id),
          path("${sample_id}.gapseq-transport.log"),
          emit: logs

    script:
    """
    set -o pipefail

    gapseq find-transport \
        -v 1 \
        -b 200 \
        -K ${task.cpus} \
        -A '${params.gapseq_aligner}' \
        -f . \
        '${faa}' \
        2>&1 | tee '${sample_id}.gapseq-transport.log'
    """

    stub:
    """
    touch '${sample_id}-Transporter.tbl'

    echo 'Stub gapseq find-transport: ${sample_id}' \
        > '${sample_id}.gapseq-transport.log'
    """
}


process RUN_GAPSEQ_DRAFT {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(reactions),
          path(pathways),
          path(transporters)

    output:
    tuple val(sample_id),
          path("${sample_id}-draft.RDS"),
          path("${sample_id}-draft.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}.gapseq-draft.log"),
          emit: logs

    script:
    """
    set -o pipefail

    gapseq draft \
        -r '${reactions}' \
        -t '${transporters}' \
        -u 200 \
        -l 100 \
        -p '${pathways}' \
        -b '${params.gapseq_taxonomy}' \
        -f . \
        2>&1 | tee '${sample_id}.gapseq-draft.log'
    """

    stub:
    """
    touch '${sample_id}-draft.RDS'
    touch '${sample_id}-draft.xml'

    echo 'Stub gapseq draft: ${sample_id}' \
        > '${sample_id}.gapseq-draft.log'
    """
}


process RUN_GAPSEQ_ADAPT {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(draft_rds),
          path(draft_xml)

    output:
    tuple val(sample_id),
          path("${sample_id}-adapt.RDS"),
          path("${sample_id}-adapt.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}.gapseq-adapt.log"),
          emit: logs

    script:
    def add_ids = params.gapseq_adapt_add?.toString()
    def remove_ids = params.gapseq_adapt_remove?.toString()

    /*
     * gapseq adapt requires at least one operation.
     *
     * If no manual adaptations have been requested, this stage acts
     * as a pass-through. This keeps the downstream workflow identical
     * regardless of whether manual curation is enabled.
     */
    if (add_ids || remove_ids) {

        def add_arg = add_ids
            ? "-a '${add_ids}'"
            : ""

        def remove_arg = remove_ids
            ? "-r '${remove_ids}'"
            : ""

        """
        set -o pipefail

        gapseq adapt \
            -m '${draft_rds}' \
            ${add_arg} \
            ${remove_arg} \
            -f . \
            2>&1 | tee '${sample_id}.gapseq-adapt.log'
        """

    } else {

        """
        cp '${draft_rds}' '${sample_id}-adapt.RDS'
        cp '${draft_xml}' '${sample_id}-adapt.xml'

        echo 'No gapseq adaptations requested; using draft unchanged.' \
            > '${sample_id}.gapseq-adapt.log'
        """
    }

    stub:
    """
    touch '${sample_id}-adapt.RDS'
    touch '${sample_id}-adapt.xml'

    echo 'Stub gapseq adapt: ${sample_id}' \
        > '${sample_id}.gapseq-adapt.log'
    """
}


process RUN_GAPSEQ_MEDIUM {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(model_rds),
          path(pathways)

    output:
    tuple val(sample_id),
          path("${sample_id}-medium.csv"),
          emit: media

    tuple val(sample_id),
          path("${sample_id}.gapseq-medium.log"),
          emit: logs

    script:
    """
    set -o pipefail

    gapseq medium \
        -m '${model_rds}' \
        -p '${pathways}' \
        -o '${sample_id}-medium.csv' \
        -f . \
        2>&1 | tee '${sample_id}.gapseq-medium.log'
    """

    stub:
    """
    echo -e 'compounds\\tname\\tmaxFlux' \
        > '${sample_id}-medium.csv'

    echo 'Stub gapseq medium: ${sample_id}' \
        > '${sample_id}.gapseq-medium.log'
    """
}


process RUN_GAPSEQ_FILL {

    tag "${sample_id}"

    input:
    tuple val(sample_id),
          path(model_rds),
          path(medium)

    output:
    tuple val(sample_id),
          path("${sample_id}.xml"),
          emit: models

    tuple val(sample_id),
          path("${sample_id}.RDS"),
          emit: rds_models

    tuple val(sample_id),
          path("${sample_id}.gapseq-fill.log"),
          emit: logs

    script:
    /*
     * gf.suite.R derives the final output name from the input model.
     *
     * gapseq doall normally gives it:
     *
     *     sample-draft.RDS
     *
     * so it produces:
     *
     *     sample.RDS
     *     sample.xml
     *
     * Our input is sample-adapt.RDS. Copy it back to the normal
     * draft-style filename so final filenames remain identical to
     * those produced by doall.
     */
    """
    set -o pipefail

    cp '${model_rds}' '${sample_id}-draft.RDS'

    gapseq fill \
        -m '${sample_id}-draft.RDS' \
        -n '${medium}' \
        -b 100 \
        -f . \
        2>&1 | tee '${sample_id}.gapseq-fill.log'
    """

    stub:
    """
    touch '${sample_id}.xml'
    touch '${sample_id}.RDS'

    echo 'Stub gapseq fill: ${sample_id}' \
        > '${sample_id}.gapseq-fill.log'
    """
}