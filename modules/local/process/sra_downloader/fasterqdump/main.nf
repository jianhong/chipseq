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
    def antibody = meta.antibody?:''
    def control  = meta.input?:meta.control?:''
    """
    mkdir -p raw
    echo "group,replicate,fastq_1,fastq_2,antibody,control,ScientificName"> designTab.csv

    if [[ "${meta.LibraryLayout.toLowerCase()}" =~ "paired" ]]; then
 	    until fasterq-dump --split-files ${ioptions.args} --outdir raw -e ${task.cpus} ${meta.Run} 2>&1 | grep -q "reads written"
     	do
     		echo "Try again for ${meta.Run}"
     	done
     	pigz -p ${task.cpus} raw/${meta.Run}*
     	cat raw/${meta.Run}*_1.fastq.gz > ${meta.Run}_R1.fastq.gz
     	cat raw/${meta.Run}*_2.fastq.gz > ${meta.Run}_R2.fastq.gz
      echo "${meta.condition},${meta.replicate},${params.outdir}/${meta.path}/${meta.Run}_R1.fastq.gz,${params.outdir}/${meta.path}/${meta.Run}_R2.fastq.gz,${antibody},${control},${meta.ScientificName}" >> designTab.csv
     else
     	until fasterq-dump --concatenate-reads ${ioptions.args} --outdir raw -e ${task.cpus} ${meta.Run} 2>&1 | grep -q "reads written"
     	do
     		echo "Try again for ${meta.Run}"
     	done
     	pigz -p ${task.cpus} raw/${meta.Run}*
     	cat raw/${meta.Run}*.gz > ${meta.Run}_R1.fastq.gz
      echo "${meta.condition},${meta.replicate},${params.outdir}/${meta.path}/${meta.Run}_R1.fastq.gz,,${antibody},${control},${meta.ScientificName}" >> designTab.csv
    fi
    
    rm -rf raw
    
    echo \$(fasterq-dump -V 2>&1) | sed 's/^.*fasterq-dump : //; s/ .*//' > sra_tools.version.txt
    """
}
