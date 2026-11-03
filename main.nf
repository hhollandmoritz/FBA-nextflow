nextflow.enable.dsl = 2

include { 
    RUN_CARVEME_DEFAULT
    RUN_CARVEME_CUSTOM
    RUN_CARVEME_NO_GAPFILL
 } from './modules/local/run_carveme'

include {
    RUN_GAPSEQ_DOALL
    RUN_GAPSEQ_FIND
    RUN_GAPSEQ_TRANSPORT
    RUN_GAPSEQ_DRAFT
    RUN_GAPSEQ_ADAPT
    RUN_GAPSEQ_MEDIUM
    RUN_GAPSEQ_FILL
} from './modules/local/run_gapseq'

include { COBRA_ANALYSIS } from './subworkflows/local/cobra_analysis'

include { BUILD_ESCHER_MAP } from './modules/local/build_escher_map'

workflow {

    log.info "Project directory: ${projectDir}"
    /* 
    * Checks on universe specifications. 
    * 1) only one of universe options specified at a time
    */
    if (params.universe && params.universe_file) {
        error "Specify either --universe or --universe_file, not both."
    }
    /*
     * Check if the universe file exists if specified.
     */ 
    if (params.universe_file) {
        file(params.universe_file, checkIfExists: true)
    }

    /*
     * Emit:
     * tuple(sample_id, faa)
     */
    genomes_ch = channel
        .fromPath(params.input, checkIfExists: true)
        .map { faa ->
            tuple(faa.baseName, faa)
        }

    /*
     * Collect models produced by all model reconstruction tools.
     *
     * Model format:
     * tuple(sample_id, medium, model)
     */
    all_models_ch = Channel.empty()

    /*
    * gapseq reconstruction
    */
    if (params.run_gapseq) {

        /*
        * Search for metabolic pathways/reactions (GAPSEQ_FIND) 
        * and transporters (GAPSEQ_TRANSPORT).
        *
        * Inputs:
        * tuple(sample_id, faa)
        */
        RUN_GAPSEQ_FIND(genomes_ch)
        RUN_GAPSEQ_TRANSPORT(genomes_ch)

        /*
        * RUN_GAPSEQ_FIND emits:
        * tuple(
        *     sample_id,
        *     reactions,
        *     pathways
        * )
        *
        * RUN_GAPSEQ_TRANSPORT emits:
        * tuple(
        *     sample_id,
        *     transporters
        * )
        *
        * Join the two channels by sample_id.
        *
        * Result:
        * tuple(
        *     sample_id,
        *     reactions,
        *     pathways,
        *     transporters
        * )
        */
        gapseq_draft_input_ch = RUN_GAPSEQ_FIND.out.evidence
            .join(
                RUN_GAPSEQ_TRANSPORT.out.evidence,
                by: 0
            )

        /*
        * Construct the initial draft model. 
        * Turn the evidence (reactions, pathways, and transporters) into 
        * a draft genome-scale metabolic model. 
        */
        RUN_GAPSEQ_DRAFT(gapseq_draft_input_ch)

        /*
        * Manual curation of draft model.
        *
        * RUN_GAPSEQ_DRAFT emits:
        * tuple(
        *     sample_id,
        *     draft_rds,
        *     draft_xml
        * )
        *
        * Inputs: tuple(sample_id, draft_rds, draft_xml).
        *
        * If no adaptations are specified in params.gapseq_adaptations, 
        * RUN_GAPSEQ_ADAPT passes the draft model through unchanged.
        */
        RUN_GAPSEQ_ADAPT(
            RUN_GAPSEQ_DRAFT.out.models
        )
    }
    

    if (params.run_carveme) {
        /*
        * CarveMe soft-constraint file.
        */
        carveme_soft_constraints_ch = Channel.value(
            params.carveme_soft_constraints
                ? file(params.carveme_soft_constraints, checkIfExists: true)
                : []
        )

        /*
        * Optional/modifiable arguments 
        * (media, extra carveme arguments)
        */
        media_config = params.media_databases ?: [:]

        /*
        * Run without gap-filling when:
        *
        * 1. No media configuration was supplied, or
        * 2. NO_GAPFILL was explicitly requested.
        */
        run_no_gapfill = media_config.isEmpty() ||
            media_config.containsKey('NO_GAPFILL')

        /*
        * Remove the special non-media condition before processing
        * actual media.
        */
        gapfill_config = media_config.findAll { medium, mediadb ->
            medium.toString() != 'NO_GAPFILL'
        }

        run_gapfill = !gapfill_config.isEmpty()

        /*
        * Non-gap-filled models
        */
        if (run_no_gapfill) {
            RUN_CARVEME_NO_GAPFILL(genomes_ch, 
            carveme_soft_constraints_ch)
        }

        /*
        * Gap-filled models
        */
        if (run_gapfill) {

            /*
            * Built-in media have no custom database.
            */
            default_media = gapfill_config
                .findAll { medium, mediadb -> mediadb == null }
                .keySet()
                .collect { medium ->
                    medium.toString()
                }

            /*
            * Custom media have an associated database file.
            */
            custom_media = gapfill_config
                .findAll { medium, mediadb -> mediadb != null }
                .collect { medium, mediadb ->
                    tuple(
                        medium.toString(),
                        file(mediadb, checkIfExists: true)
                    )
                }

            /*
            * Emit:
            * tuple(sample_id, faa, medium)
            */
            default_jobs_ch = genomes_ch
                .combine(channel.fromList(default_media))
                .map { sample_id, faa, medium ->
                    tuple(sample_id, faa, medium)
                }

            /*
            * Emit:
            * tuple(sample_id, faa, medium, mediadb)
            */
            custom_jobs_ch = genomes_ch
                .combine(channel.fromList(custom_media))
                .map { sample_id, faa, medium, mediadb ->
                    tuple(sample_id, faa, medium, mediadb)
                }

            RUN_CARVEME_DEFAULT(
                default_jobs_ch,
                carveme_soft_constraints_ch
                )
            RUN_CARVEME_CUSTOM(
                custom_jobs_ch,
                carveme_soft_constraints_ch
                )

            gapfill_models_ch = RUN_CARVEME_DEFAULT.out.models
                .mix(RUN_CARVEME_CUSTOM.out.models)

            gapfill_logs_ch = RUN_CARVEME_DEFAULT.out.logs
                .mix(RUN_CARVEME_CUSTOM.out.logs)
        }

        /*
        * Create unified output channels.
        */
        if (run_no_gapfill && run_gapfill) {

            models_ch = RUN_CARVEME_NO_GAPFILL.out.models
                .mix(gapfill_models_ch)

            logs_ch = RUN_CARVEME_NO_GAPFILL.out.logs
                .mix(gapfill_logs_ch)

        } else if (run_no_gapfill) {

            models_ch = RUN_CARVEME_NO_GAPFILL.out.models
            logs_ch   = RUN_CARVEME_NO_GAPFILL.out.logs

        } else {

            models_ch = gapfill_models_ch
            logs_ch   = gapfill_logs_ch
        }
        /*
         * Add CarveMe reconstructions to the
         * model channel for downstream analysis.
         */
        all_models_ch = all_models_ch.mix(models_ch)
    }

    /*
     * COBRApy analysis
     */
    if (params.run_cobra) {
        cobra_medium_ch = Channel.value(
            params.cobra_medium
                ? file(params.cobra_medium, checkIfExists: true)
                : []
        )
        /* 
        * Hannah note; instead of this being an error, can we make it a warning that outputs will not be transformed
        */ 
        if (!params.modelseed_db) {
            error "COBRA summarization requires --modelseed_db."
        }

        modelseed_db_ch = Channel.value(
            file(
                params.modelseed_db,
                checkIfExists: true
            )
        )

        cobra_run_options_ch = Channel.value([
            method                           : params.cobra_method,
            fraction_of_optimum              : params.cobra_fraction_of_optimum,
            objective                        : params.cobra_objective,
            solver                           : params.cobra_solver,
            medium_mode                      : params.cobra_medium_mode,
            ignore_missing_medium_reactions  : params.cobra_ignore_missing_medium_reactions,
            flux_threshold                   : params.cobra_flux_threshold,
            fail_on_nonoptimal               : params.cobra_fail_on_nonoptimal,
            extra_args                       : params.cobra_run_args
        ])
        

        COBRA_ANALYSIS(
            all_models_ch,
            cobra_medium_ch,
            modelseed_db_ch,
            cobra_run_options_ch,
        )
    }

    /*
     * Build Escher maps
     */
    if (params.escher_map) {
        escher_map_ch = Channel.value(
            file(params.escher_map, checkIfExists: true)
        )
        BUILD_ESCHER_MAP(
            all_models_ch,
            escher_map_ch
        )
    }

}
