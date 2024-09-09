// SCMSEQ_NEW, first version
nextflow.enable.dsl = 2


// Validate input parameters
WorkflowMain.initialise(workflow, params, log)

// Include here
include { long_reads } from "./subworkflows/long_reads/main"


//----------------------------------------------------------------------------//
// SCM-seq pipeline entry points
//----------------------------------------------------------------------------//

//

workflow LONG_READS {
    // Create a Channel of nanopore fastqs from the input csv file
    input = extract_csv(file(params.input))
    long_reads(input)
    long_reads.out.results.view()

}


// extract channels from input annotation sample sheet 
def extract_csv(csv_file) {
    // check that the sample sheet is not 1 line or less, because it'll skip all subsequent checks if so.
    file(csv_file).withReader('UTF-8') { reader ->
        def line, numberOfLinesInSampleSheet = 0;
        while ((line = reader.readLine()) != null) {
            numberOfLinesInSampleSheet++
            if (numberOfLinesInSampleSheet == 1){
                def requiredColumns = ["sample", 'lane', 'fastq']
                def headerColumns = line
                if (!requiredColumns.every { headerColumns.contains(it) }) {
                    log.error "Header missing or CSV file does not contain all of the required columns in the header: ${requiredColumns}"
                    System.exit(1)
                }
            }
        }
        
        if (numberOfLinesInSampleSheet < 2) {
            log.error "Provided SampleSheet has less than two lines. Provide a samplesheet with header and at least a sample."
            System.exit(1)
        }
    }

    Channel.from(csv_file)
        .splitCsv(header: true)
        .map { row -> [row.sample, row.lane, row.fastq]
        }
}