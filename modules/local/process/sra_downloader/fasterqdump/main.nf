// Import generic module functions
include { initOptions; saveFiles } from '../../functions'
params.options = [:]
/*
 * Create index.html
 */
process JO_FASTQDUMP {
    tag "${meta.Run}"
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:"${meta.path}", publish_id:"${meta.Run}") }

    conda (params.conda ? "${params.conda_softwares.sra_tools} ${params.conda_softwares.pigz}" : null)
    
    input:
    val meta

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastq
    path "designTab.csv", emit: design
    path "*.version.txt", emit: version
    
    script:
    def ioptions = initOptions(params.options)
    def genomes=["homo_sapiens":"GRCh38", "mus_musculus":"GRCm38", "arabidopsis_thaliana":"TAIR10", "bacillus_subtilis_168":"EB2", "bos_taurus":"UMD3.1", "caenorhabditis_elegans":"WBcel235", "canis_familiaris":"CanFam3.1", "danio_rerio":"GRCz10", "Drosophila_melanogaster":"BDGP6", "equus_caballus":"EquCab2", "gallus_gallus":"Galgal4", "escherichia_coli":"EB1", "gallus_gallus":"Galgal4", "glycine_max":"Gm01", "macaca_mulatta":"Mmul_1", "oryza_sativa_japonica":"IRGSP-1.0", "pan_troglodytes":"CHIMP2.1.4", "rattus_norvegicus":"Rnor_6.0", "saccharomyces_cerevisiae":"R64-1-1", "schizosaccharomyces_pombe":"EF2", "sorghum_bicolor":"Sbi1", "sus_scrofa":"Sscrofa10.2", "zea_mays":"AGPv3"]
    def genome=genomes[meta.ScientificName.replaceAll("\\W+", "_").toLowerCase()]
    """
    mkdir -p raw
    echo "group,replicate,fastq_1,fastq_2,genome"> designTab.csv

    if [[ "${meta.LibraryLayout.toLowerCase()}" =~ "paired" ]]; then
 	    until fasterq-dump --split-files ${ioptions.args} --outdir raw -e ${task.cpus} ${meta.Run} 2>&1 | grep -q "reads written"
     	do
     		echo "Try again for ${meta.Run}"
     	done
     	pigz -p ${task.cpus} raw/${meta.Run}*
     	cat raw/${meta.Run}*_1.fastq.gz > ${meta.Run}_R1.fastq.gz
     	cat raw/${meta.Run}*_2.fastq.gz > ${meta.Run}_R2.fastq.gz
      echo "${meta.condition},${meta.replicate},${params.outdir}/${meta.path}/${meta.Run}_R1.fastq.gz,${params.outdir}/${meta.path}/${meta.Run}_R2.fastq.gz,${genome}" >> designTab.csv
     else
     	until fasterq-dump --concatenate-reads ${ioptions.args} --outdir raw -e ${task.cpus} ${meta.Run} 2>&1 | grep -q "reads written"
     	do
     		echo "Try again for ${meta.Run}"
     	done
     	pigz -p ${task.cpus} raw/${meta.Run}*
     	cat raw/${meta.Run}*.gz > ${meta.Run}_R1.fastq.gz
      echo "${meta.condition},${meta.replicate},${params.outdir}/${meta.path}/${meta.Run}_R1.fastq.gz,,${genome}" >> designTab.csv
    fi
    
    rm -rf raw
    
    echo \$(fasterq-dump -V 2>&1) | sed 's/^.*fasterq-dump : //; s/ .*//' > sra_tools.version.txt
    """
}
