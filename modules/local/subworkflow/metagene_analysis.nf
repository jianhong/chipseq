/*
 * Metagene analysis
 */

include { JO_MERGE_REP_BAM   } from '../process/mergebam/bam_merge'
include { JO_METAGENE        } from '../process/metagene/metagene'

workflow JO_METAGENE_ANALYSIS {
    take:
    ch_clean_bam               // channel: [ [meta], bam ]
    ch_sigle_bw                // channel: [ [meta], bw  ]
    ch_bed                     // channel: [ bed ]
    bam_merge_options          //     map: options for bam_merge module
    metagene_options           //     map: options for metagene module

    main:
    
    ch_clean_bam
        .map {
            meta, bam ->
                [ meta.control? 
                     meta.control.replaceAll(/_R\d+.*$/, ""):
                     meta.id.replaceAll(/_R\d+.*$/, ""),
                  meta.id.replaceAll(/_R\d+.*$/, ""), meta.single_end, bam ] }
       .groupTuple(by: [0, 1, 2])
       .map{control, name, se, bam -> [control, [name, se, bam]]}
       .set{ch_to_be_merged}
    ch_clean_bam
        .map {
            meta, bam ->
                [ meta.id.replaceAll(/_R\d+.*$/, ""), bam ] }
       .groupTuple(by: [0])
       .set{ch_to_be_merged_input}
    ch_to_be_merged_input.cross(ch_to_be_merged)
       .map{input, chip -> 
             if(input[0]==chip[1][0]){
                [chip[1][0], null, chip[1][1], chip[1][2].unique(), []]
             }else{
                [chip[1][0], input[0], chip[1][1], chip[1][2].unique(), input[1].unique()]
             }
           }
       .set{ch_to_be_merged}
    
    JO_MERGE_REP_BAM(ch_to_be_merged,bam_merge_options)
    JO_MERGE_REP_BAM.out.bw
       .map{ name, bw -> 
               [name, bw.findAll{it.toString().endsWith('.CPM.bw')}] }
       .collect().ifEmpty([[], []])
       .set{ch_cpm_bw}
    JO_MERGE_REP_BAM.out.bw
       .map{ name, bw -> 
               [name, bw.findAll{
                           it.toString()
                             .endsWith('.normByInput.CPM.log2ratio.bw')}] }
       .collect().ifEmpty([[], []])
       .set{ch_log2_bw}
    ch_sigle_bw
       .map{ meta, bw -> [meta.id, bw] }
       .collect().ifEmpty([[], []])
       .set{ch_sigle_bw}
    ch_cpm_bw.concat(ch_log2_bw, ch_sigle_bw).set{ch_bw}
    JO_METAGENE(ch_bw, ch_bed, metagene_options)
    
    emit:
    bam = JO_MERGE_REP_BAM.out.bam  // channel: [ val(name), path(bam), path(bai) ]
    bw  = JO_MERGE_REP_BAM.out.bw   // channel: [ val(name), [bw] ]
}
