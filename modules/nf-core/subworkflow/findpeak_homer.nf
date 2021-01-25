/*
 * Picard MarkDuplicates, sort, index BAM file and run samtools stats, flagstat and idxstats
 */

include { HOMER_MAKETAGDIRECTORY } from '../software/homer/maketagdirectory/main'
include { HOMER_FINDPEAKS        } from '../software/homer/findpeaks/main'
include { HOMER_ANNOTATEPEAKS    } from '../software/homer/annotatepeaks/main'
include { HOMER_POS2BED          } from '../software/homer/pos2bed/main'

workflow HOMER_CALLPEAK {
    take:
    ch_bam                        // channel: [ val(meta), [ bam ], [ controlbam ] ]
    fasta                         // channel: path
    gtf                           // channel: path
    maketagdirectory_options      // map: options for maketagdirectory
    findpeaks_options             // map: options for findpeaks
    annotationpeaks_options       // map: options for annotationpeaks
    pos2bed_options               // map: options for pos2bed
    homer_callpeak_options        // map: options


    main:
    HOMER_MAKETAGDIRECTORY(ch_bam, maketagdirectory_options)
    HOMER_FINDPEAKS(HOMER_MAKETAGDIRECTORY.out.tagdir, findpeaks_options, homer_callpeak_options.publish_dir)
    HOMER_POS2BED(HOMER_FINDPEAKS.out.peak, pos2bed_options, homer_callpeak_options.publish_dir)
    HOMER_ANNOTATEPEAKS(HOMER_FINDPEAKS.out.peak, fasta, gtf, annotationpeaks_options, homer_callpeak_options.publish_dir)

    emit:
    peak = HOMER_FINDPEAKS.out.peak                    // channel: [ val(meta), [ peak ] ]
    annopeaks = HOMER_ANNOTATEPEAKS.out.txt            // channel: [ val(meta), [ annotatePeaks ] ]
    bed = HOMER_POS2BED.out.bed                        // channel: [ val(meta), [ bed ] ]
    homer_version = HOMER_MAKETAGDIRECTORY.out.version //    path: *.version.txt
}
