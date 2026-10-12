// SCMSEQ_NEW, first version
nextflow.enable.dsl = 2


// Validate input parameters
WorkflowMain.initialise(workflow, params, log)

// Include here
include { long_reads; long_reads_from_tagged_bam } from "./subworkflows/long_reads/main"


//----------------------------------------------------------------------------//
// SCM-seq pipeline entry points
//----------------------------------------------------------------------------//

//

workflow LONG_READS {
    if (params.input_type == 'tagged_bam') {
        // Tagged, aligned BAMs per sample start the run at SPLIT_BAM
        input = extract_tagged_bam_csv(file(params.input))
        long_reads_from_tagged_bam(input)
        long_reads_from_tagged_bam.out.results.view()
    } else {
        // Create a Channel of nanopore fastqs from the input csv file
        input = extract_csv(file(params.input))
        long_reads(input)
        long_reads.out.results.view()
    }

}


// extract channels from input annotation sample sheet
def extract_csv(csv_file) {
    check_samplesheet(csv_file, ["sample", 'lane', 'fastq'])

    Channel.from(csv_file)
        .splitCsv(header: true)
        .map { row -> [row.sample, row.lane, row.fastq]
        }
}


// extract [sample, bam] from a tagged BAM sample sheet, one row per sample
def extract_tagged_bam_csv(csv_file) {
    check_samplesheet(csv_file, ["sample", "bam"])

    def rows = csv_file.splitCsv(header: true)
    rows.countBy { it.sample }.findAll { sample, n -> n > 1 }.each { sample, n ->
        log.error "Sample ${sample} is listed more than once in ${csv_file}. Give one tagged BAM per sample."
        System.exit(1)
    }

    Channel.fromList(rows.collect { row -> [row.sample, file(row.bam, checkIfExists: true)] })
}


def check_samplesheet(csv_file, requiredColumns) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
            if (numberOfLinesInSampleSheet == 1){
                def headerColumns = line.tokenize(',')*.trim()
                if (!requiredColumns.every { headerColumns.contains(it) }) {
                    log.error "Header missing or CSV file does not contain all of the required columns in the header: ${requiredColumns} (--input_type ${params.input_type})"
                    System.exit(1)
                }
            }
        }

        if (numberOfLinesInSampleSheet < 2) {
            log.error "Provided SampleSheet has less than two lines. Provide a samplesheet with header and at least a sample."
            System.exit(1)
        }
    }
}
