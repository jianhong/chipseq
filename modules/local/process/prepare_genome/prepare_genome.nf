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
include { GTF2BED                     } from '../gtf2bed/gtf2bed'                        addParams( options: params.genome_options       )
include { GET_CHROM_SIZES             } from '../get_chrom_sizes'                        addParams( options: params.genome_options       )
include { MAKE_GENOME_FILTER          } from '../make_genome_filter'                     addParams( options: params.genome_options       )
include { GFFREAD                     } from '../../../nf-core/software/gffread/main'    addParams( options: params.gffread_options      )
include { BWA_INDEX                   } from '../../../nf-core/software/bwa/index/main'  addParams( options: params.modules['bwa_index'])

workflow PREPARE_GENOME {
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
    ch_gtf = Channel.empty()
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
        ch_gffread_version = GFFREAD.version
    }

    /*
     * Uncompress gene BED annotation file or create from GTF if required
     */
    ch_gene_bed = Channel.empty()
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
    ch_blacklist = ch_blacklist.ifEmpty([])

    filtered_bed = MAKE_GENOME_FILTER (
        ch_chrom_sizes,
        ch_blacklist
    ).bed

    /*
     * Uncompress bwa index or generate from scratch if required
     */ 
    ch_bwa_index = params.bwa_index ? Channel.value(file(params.bwa_index)) : BWA_INDEX ( ch_fasta ).index

    emit:
    fasta             = ch_fasta                       // path: genome.fasta,
    gtf               = ch_gtf                         // path: genome.gtf,
    gene_bed          = ch_gene_bed                    //           path: gene.bed,
    chrom_sizes       = ch_chrom_sizes                 //           path: genome.sizes,
    blacklist         = ch_blacklist                   //           path: blacklist.bed,
    bed               = filtered_bed                   //           path: *.bed,
    bwa_index         = ch_bwa_index                   //           path: fasta.*]]
    data              = [fasta:ch_fasta,
                         gtf:ch_gtf,
                         gene_bed:ch_gene_bed,
                         chrom_sizes:ch_chrom_sizes,
                         blacklist:ch_blacklist,
                         bed:filtered_bed,
                         bwa:ch_bwa_index] //full data for one genome
    gffread_version   = ch_gffread_version             // path: *.version.txt
    filter_version    = MAKE_GENOME_FILTER.out.version // path: *.version.txt
}
