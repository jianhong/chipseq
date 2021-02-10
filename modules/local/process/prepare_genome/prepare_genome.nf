/*
 * Uncompress and prepare reference genome files
*/

params.genome_options       = [:]
params.gffread_options      = [:]

include {
    GUNZIP as GUNZIP_FASTA
    GUNZIP as GUNZIP_GTF
    GUNZIP as GUNZIP_GFF
    GUNZIP as GUNZIP_GENE_BED
    GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../utilities/gunzip'                       addParams( options: params.genome_options       )
include { GTF2BED                     } from '../gtf2bed'                                addParams( options: params.genome_options       )
include { CAT_ADDITIONAL_FASTA        } from '../cat_additional_fasta'                   addParams( options: params.genome_options       )
include { GTF_GENE_FILTER             } from '../gtf_gene_filter'                        addParams( options: params.genome_options       )
include { GET_CHROM_SIZES             } from '../get_chrom_sizes'                        addParams( options: params.genome_options       )
include { MAKE_GENOME_FILTER          } from '../make_genome_filter'                     addParams( options: params.genome_options       )
include { GFFREAD                     } from '../../../nf-core/software/gffread/main'    addParams( options: params.gffread_options      )
include { BWA_INDEX                   } from '../../../nf-core/software/bwa/index/main'  addParams( options: params.modules['bwa_index'])

workflow PREPARE_GENOME {
    take:
    genome    // map: options for genome [genome, fasta, gtf, gene_bed, bwa_index, macs_gsize, deep_gsize, blacklist, readme, species]
    
    main:
    
    /*
     * Uncompress genome fasta file if required
     */
    if (params.fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file(params.fasta)
    }

    /*
     * Uncompress GTF annotation file or create from GFF3 if required
     */
    ch_gffread_version = Channel.empty()
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF ( params.gtf ).gunzip
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff = GUNZIP_GFF ( params.gff ).gunzip
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf = GFFREAD ( ch_gff ).gtf
        ch_gffread_version = GFFREAD.out.version
    }

    /*
     * Uncompress additional fasta file and concatenate with reference fasta and gtf files
     */
    /*if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( params.additional_fasta ).gunzip
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta )
        ch_fasta = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf   = CAT_ADDITIONAL_FASTA.out.gtf
    }*/

    /*
     * Uncompress gene BED annotation file or create from GTF if required
     */
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( params.gene_bed ).gunzip
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf )
    }

    /*
     * Create chromosome sizes file
     */
    ch_chrom_sizes = GET_CHROM_SIZES ( ch_fasta ).sizes
    
    if (params.blacklist) { 
        ch_blacklist = Channel.fromPath(params.blacklist, checkIfExists: true) 
    } else { ch_blacklist = Channel.empty() }

    MAKE_GENOME_FILTER (
        ch_chrom_sizes,
        ch_blacklist.ifEmpty([])
    )

    /*
     * Uncompress bwa index or generate from scratch if required
     */ 
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta ).index

    emit:
    fasta            = ch_fasta            // path: genome.fasta
    gtf              = ch_gtf              // path: genome.gtf
    gene_bed         = ch_gene_bed         // path: gene.bed
    chrom_sizes      = ch_chrom_sizes      // path: genome.sizes
    filter_bed       = MAKE_GENOME_FILTER.out.bed   // path: *.bed
    ch_index         = ch_bwa_index        // path: fasta.*
    blacklist        = ch_blacklist        // path: blacklist.bed
    gffread_version  = ch_gffread_version  // path: *.version.txt
    filter_version   = MAKE_GENOME_FILTER.out.version // path: *.version.txt
}
