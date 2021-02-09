/*
 * Check input samplesheet and get read channels
 */
params.options = [:]

include { JO_SRAINFO           } from './sra_info/main' addParams(options: params.options)
include { JO_FASTQDUMP         } from './fasterqdump/main' addParams(options: params.options)

workflow JO_SRADOWNLOADER {
    take:
    id // GEO accession ID or SRA project ID
    
    main:
    // get srr id and information to create meta values
    JO_SRAINFO(id)

    //dowload data
    JO_FASTQDUMP(JO_SRAINFO.out.info.splitCsv ( header:true, sep:',' ).map{
        meta -> 
            path = "fastq/${meta.ScientificName.replaceAll('\\W+', '_')}/${meta.LibraryLayout}/${meta.LibraryStrategy}/${meta.condition}/${meta.replicate}"
            meta
    })

    emit:
    reads = JO_FASTQDUMP.out.fastq             // channel: [ val(meta), [ reads ] ]
    design = JO_FASTQDUMP.out.design           // channel: path
    entrez_direct = JO_SRAINFO.out.version     // channel: path
    sra_tools = JO_FASTQDUMP.out.version       // channel: path
}
