// Import generic module functions
include { saveFiles; getRealPath } from '../functions'

/*
 * Reformat design file, check validitiy and create IP vs control mappings
 */
process SAMPLESHEET_CHECK {
    tag "$samplesheet"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:options, publish_dir:"pipeline_info", publish_id:'') }

    input:
    path samplesheet
    val options

    output:
    path '*.csv'

    script:  // This script is bundled with the pipeline, in nf-core/chipseq/bin/
    def curr_path = getRealPath()
    """
    ${curr_path}/check_samplesheet/check_samplesheet.py $samplesheet samplesheet.valid.csv
    """
}

// Function to get list of [ meta, [ fastq_1, fastq_2 ] ]
def get_samplesheet_paths(LinkedHashMap row, String seq_center, String genome) {
    def meta = [:]
    meta.id = row.sample
    meta.single_end = row.single_end.toBoolean()
    meta.antibody = row.antibody
    meta.control = row.control
    meta.md5 = [row.md5_1, row.md5_2]
    meta.peaktype = row.peaktype
    meta.ctrl = row.control
    meta.group = row.sample.split('_')[0..-2].join('_')
    meta.techrep = row.sample.replaceAll("^.*_T(.*?)\$", "T\$1")
    if(row.treatment){
        meta.treatment = row.treatment
    }
    meta.sample_uniq_keys = ['md5', 'techrep', 'read_group']
    
    def genomes=["homo_sapiens":"GRCh38", "mus_musculus":"GRCm38", "arabidopsis_thaliana":"TAIR10", "bacillus_subtilis_168":"EB2", "bos_taurus":"UMD3.1", "caenorhabditis_elegans":"WBcel235", "canis_familiaris":"CanFam3.1", "danio_rerio":"GRCz10", "Drosophila_melanogaster":"BDGP6", "equus_caballus":"EquCab2", "gallus_gallus":"Galgal4", "escherichia_coli":"EB1", "gallus_gallus":"Galgal4", "glycine_max":"Gm01", "macaca_mulatta":"Mmul_1", "oryza_sativa_japonica":"IRGSP-1.0", "pan_troglodytes":"CHIMP2.1.4", "rattus_norvegicus":"Rnor_6.0", "saccharomyces_cerevisiae":"R64-1-1", "schizosaccharomyces_pombe":"EF2", "sorghum_bicolor":"Sbi1", "sus_scrofa":"Sscrofa10.2", "zea_mays":"AGPv3"]
    
    if(row.ScientificName){
        meta.genome = genomes[row.ScientificName.toString().replaceAll("\\W+", "_").toLowerCase()]
    }else{
        meta.genome = genome
    }

    def rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\'"
    if (seq_center) {
        rg = "\'@RG\\tID:${meta.id}\\tSM:${meta.id.split('_')[0..-2].join('_')}\\tPL:ILLUMINA\\tLB:${meta.id}\\tPU:1\\tCN:${seq_center}\'"
    }
    meta.read_group = rg

    def array = []
    if (meta.single_end) {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true) ] ]
    } else {
        array = [ meta, [ file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ] ]
    }
    return array
}
