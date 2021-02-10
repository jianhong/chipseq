#!/usr/bin/env nextflow
/*
========================================================================================
                         qiubio-nf-core/chipseq
========================================================================================
 qiubio-nf-core/chipseq Analysis Pipeline.
 #### Homepage / Documentation
 https://github.com/jianhong/chipseq
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

def json_schema = "$projectDir/nextflow_schema.json"

/*
 * Print help message if required
 */
if (params.help) {
    def command = "nextflow run jianhong/chipseq -r dev --input design.csv --genome GRCh37 -profile docker"
    log.info Headers.nf_core(workflow, params.monochrome_logs)
    log.info Schema.params_help(json_schema, command)
    exit 0
}

////////////////////////////////////////////////////
/* --         PRINT PARAMETER SUMMARY          -- */
////////////////////////////////////////////////////
/*
 * Print parameter summary
 */
def summary_params = Schema.params_summary_map(workflow, params, json_schema)
log.info Schema.params_summary_log(workflow, params, json_schema)


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

/*
 * Validate input
 */
 
workflow {
    if (params.input)     { 
        if(["GSE", "SRP", "PRJ"].contains(params.input.substring(0, 3)) & params.input[-3..-1].isNumber()){
            /*
             * SUBWORKFLOW: Get SRA run information for public database ids, download and md5sum check FastQ files, auto-create samplesheet
             */
            include { JO_SRADOWNLOADER } from './modules/local/process/sra_downloader/main' addParams(options: [:])
            id = Channel.of(params.input)
            JO_SRADOWNLOADER(id)
            JO_SRADOWNLOADER.out.design
              .collectFile(keepHeader:true, name:"desginTab.csv", newLine: false)
              .set{ch_input}
        }else{
            ch_input = file(params.input, checkIfExists: true)
        }
    } else { exit 1, 'Samples design file not specified!' }

    /*
     * SUBWORKFLOW: Run main jianhong/chipseq analysis pipeline
     */
    include { CHIPSEQ } from './chipseq' addParams( summary_params: summary_params )
    //CHIPSEQ (ch_input)
}


