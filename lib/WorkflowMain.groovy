//
// This file holds several functions specific to the main.nf workflow in the nf-core/variantalker pipeline
//

class WorkflowMain {

    //
    // Retired parameters and their migration messages. Passing one warns and
    // the run continues; the value is ignored.
    //
    public static final Map RETIRED_PARAMS = [
        window_width            : 'UMIs are no longer grouped in user-set genomic windows.',
        min_dist_umi_collapsing : 'the UMI clustering distance is no longer user-set.',
        baseline_error          : 'renamed; use --error_floor instead.',
        min_base_quality        : 'renamed; use --snv_min_base_quality instead (it applies to SNV targets only).',
    ]

    //
    // Warn about retired parameters with a migration message
    //
    public static void warnRetiredParams(params, log) {
        RETIRED_PARAMS.each { name, migration ->
            if (spellings(name).any { params.containsKey(it) }) {
                log.warn "--${name} is retired and its value is ignored: ${migration}"
            }
        }
    }

    //
    // Every name Nextflow may store a snake_case param under (--a_b, --a-b, aB)
    //
    public static List retiredParamSpellings() {
        return RETIRED_PARAMS.keySet().collectMany { spellings(it) }
    }

    private static List spellings(String name) {
        def camel = name.replaceAll(/_(\w)/) { all, c -> c.toUpperCase() }
        return [name, name.replace('_', '-'), camel]
    }

    //
    // Citation string for pipeline
    //
    public static String citation(workflow) {
        return "If you use ${workflow.manifest.name} for your analysis please cite:\n\n" +
            // TODO nf-core: Add Zenodo DOI for pipeline after first release
            //"* The pipeline\n" +
            //"  https://doi.org/10.5281/zenodo.XXXXXXX\n\n" +
            "* The nf-core framework\n" +
            "  https://doi.org/10.1038/s41587-020-0439-x\n\n" +
            "* Software dependencies\n" +
            "  https://github.com/${workflow.manifest.name}/blob/master/CITATIONS.md"
    }

    //
    // Generate help string
    //
    public static String help(workflow, params, log) {
        def command = "nextflow run main.nf -c <config> -profile <profile> -entry LONG_READS"
        def help_string = ''
        help_string += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        help_string += NfcoreSchema.paramsHelp(workflow, params, command)
        help_string += '\n' + citation(workflow) + '\n'
        help_string += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return help_string
    }

    //
    // Generate parameter summary log string
    //
    public static String paramsSummaryLog(workflow, params, log) {
        def summary_log = ''
        summary_log += NfcoreTemplate.logo(workflow, params.monochrome_logs)
        summary_log += NfcoreSchema.paramsSummaryLog(workflow, params)
        summary_log += '\n' + citation(workflow) + '\n'
        summary_log += NfcoreTemplate.dashedLine(params.monochrome_logs)
        return summary_log
    }

    //
    // Validate parameters and print summary to screen
    //
    public static void initialise(workflow, params, log) {
        // Print help to screen if required
        if (params.help) {
            log.info help(workflow, params, log)
            System.exit(0)
        }

        // Print workflow version and exit on --version
        if (params.version) {
            String workflow_version = NfcoreTemplate.version(workflow)
            log.info "${workflow.manifest.name} ${workflow_version}"
            System.exit(0)
        }

        // Print parameter summary log to screen
        log.info paramsSummaryLog(workflow, params, log)

        // Retired parameters warn here, so validation doesn't also report them as unexpected
        warnRetiredParams(params, log)

        // Validate workflow parameters via the JSON schema
        if (params.validate_params) {
            NfcoreSchema.validateParameters(workflow, params, log, retiredParamSpellings())
        }

        // Check that a -profile or Nextflow config has been provided to run the pipeline
        NfcoreTemplate.checkConfigProvided(workflow, log)

        // Check that conda channels are set-up correctly
        if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
            Utils.checkCondaChannels(log)
        }

        // Check AWS batch settings
        NfcoreTemplate.awsBatch(workflow, params)

        // Check input has been provided
        if (!params.input) {
            log.error "Please provide an input samplesheet to the pipeline e.g. '--input samplesheet.csv'"
            System.exit(1)
        }

        // Only the fastq entry point aligns reads, so only it needs the minimap2 index
        if (params.input_type == 'fastq' && !params.minimap2_mindex) {
            log.error "--minimap2_mindex is required with --input_type fastq"
            System.exit(1)
        }
    }
    //
    // Get attribute from genome config file e.g. fasta
    //
    public static Object getGenomeAttribute(params, attribute) {
        if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
            if (params.genomes[ params.genome ].containsKey(attribute)) {
                return params.genomes[ params.genome ][ attribute ]
            }
        }
        return null
    }
}
