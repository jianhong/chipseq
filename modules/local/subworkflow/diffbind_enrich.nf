/*
 * Metagene analysis
 */

include { JO_DIFFBIND            } from '../process/diffbind/diffbind'
include { JO_ENRICHMENTANALYSIS  } from '../process/enrichment_analysis/enrich'

workflow JO_DIFFBIND_ENRICHMENT {
    take:
    ch_cbam_peak                        // channel: [ [meta], bam, peaks ]
    gtf                                 // channel: [ gtf ]
    blacklist                           // channel: [ blacklist ]
    diffbind_options                    // map: options for diffbind module

    main:
    JO_DIFFBIND (
            ch_cbam_peak,
            gtf,
            blacklist.ifEmpty([]),
            diffbind_options
        )
    JO_ENRICHMENTANALYSIS (
            JO_DIFFBIND.out.res,
            diffbind_options
        )
    
    emit:
    res     = JO_DIFFBIND.out.res                    // channel: [ 'DiffBind/*' ]
    enrich  = JO_ENRICHMENTANALYSIS.out.enrichment   // channel: [ 'enrichment/*' ]
}
