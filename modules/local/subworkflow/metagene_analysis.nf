/*
 * Metagene analysis
 */

include { JO_MERGE_REP_BAM   } from '../process/mergebam/bam_merge'
include { JO_METAGENE        } from '../process/metagene/metagene'

workflow JO_METAGENE_ANALYSIS {
    take:
    ch_clean_bam               // channel: [ [meta], bam ]
    ch_single_bw                // channel: [ [meta], bw  ]
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
                  meta.id.replaceAll(/_R\d+.*$/, ""), meta.single_end, meta.antibody, bam ] }
       .groupTuple(by: [0, 1, 2, 3])
       .map{control, name, se, ab, bam -> [control, [name, se, ab, bam]]}
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
                [['id':chip[1][0], 
                  'single_end':chip[1][1],
                  'antibody':chip[1][2], 
                  'control':null],
                 chip[1][3].unique(), []]
             }else{
                [['id':chip[1][0], 
                  'single_end':chip[1][1],
                  'antibody':chip[1][2], 
                  'control':input[0]], 
                 chip[1][3].unique(), input[1].unique()]
             }
           }
       .set{ch_to_be_merged}
    ch_to_be_merged.view()
    
    JO_MERGE_REP_BAM(ch_to_be_merged,bam_merge_options)
    JO_MERGE_REP_BAM.out.bw
       .map{ meta, bw -> 
               ["CPM", meta.id, bw.findAll{it.toString().endsWith('.CPM.bw')}] }
       .set{ch_cpm_bw}
    JO_MERGE_REP_BAM.out.bw
       .map{ meta, bw -> 
               ["log2", meta.id, bw.findAll{
                           it.toString()
                             .endsWith('.normByInput.CPM.log2ratio.bw')}] }
       .set{ch_log2_bw}
    ch_single_bw.map{meta, bw ->
                       ["single", meta.id, bw]}
       .set{ch_single_bw}
    ch_cpm_bw.concat(ch_log2_bw, ch_single_bw)
       .groupTuple(by: [0]).map{[it[1], it[2].flatten()]}
       .set{ch_bw}
    JO_METAGENE(ch_bw, ch_bed, metagene_options)
    
    emit:
    bam = JO_MERGE_REP_BAM.out.bam  // channel: [ val(meta), [bam], [bai] ]
    bw  = JO_MERGE_REP_BAM.out.bw   // channel: [ val(meta), [bw] ]
}
