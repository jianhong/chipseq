# ![qiubio-nf-core/chipseq](assets/chipseqlogo.png)

## Introduction

_qiubio/chipseq_ is a bioinformatics analysis pipeline used for Chromatin ImmunopreciPitation sequencing (ChIP-seq) data based on [nfcore/chipseq](https://nf-co.re/chipseq).

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It comes with docker containers making installation trivial and results highly reproducible.

This pipeline will generate the UCSC genome browser track hub and metagene analysis
resuls in addition to original output of _nfcore/chipseq_ pipeline.

The most important change for this pipeline from _nfcore/chipseq_ pipeline is
that I tried to improve the reproducibility of the pipeline depend on conda but
not docker for following 2 reasons:

1. I can not use docker in our cluster.

2. The memory required for the pipeline is too heavy for personal computer if using docker.

However, conda always throw errors when create environment even I use modules.
I add conda_softwares section in module.conf setting to make the pipeline more flexible
to figure out this issue.
I also changed the R/Biocondactor package installation methods from conda installation
to BiocManager installation, which will be much slower than conda installation.
The reason for that is because lots of package in conda is malfunction. By using
BiocManager to avoid the dependece issues.

## Pipeline summary

1. Raw read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Adapter trimming ([`Trim Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/))
3. Alignment ([`BWA`](https://sourceforge.net/projects/bio-bwa/files/))
4. Mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
5. Merge alignments from multiple libraries of the same sample ([`picard`](https://broadinstitute.github.io/picard/))
    1. Re-mark duplicates ([`picard`](https://broadinstitute.github.io/picard/))
    2. Filtering to remove:
        * reads mapping to blacklisted regions ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/), [`BEDTools`](https://github.com/arq5x/bedtools2/))
        * reads that are marked as duplicates ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that arent marked as primary alignments ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that are unmapped ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads that map to multiple locations ([`SAMtools`](https://sourceforge.net/projects/samtools/files/samtools/))
        * reads containing > 4 mismatches ([`BAMTools`](https://github.com/pezmaster31/bamtools))
        * reads that have an insert size > 2kb ([`BAMTools`](https://github.com/pezmaster31/bamtools); *paired-end only*)
        * reads that map to different chromosomes ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads that arent in FR orientation ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
        * reads where only one read of the pair fails the above criteria ([`Pysam`](http://pysam.readthedocs.io/en/latest/installation.html); *paired-end only*)
    3. Alignment-level QC and estimation of library complexity ([`picard`](https://broadinstitute.github.io/picard/), [`Preseq`](http://smithlabresearch.org/software/preseq/))
    4. Create normalised bigWig files scaled to 1 million mapped reads ([`BEDTools`](https://github.com/arq5x/bedtools2/), [`bedGraphToBigWig`](http://hgdownload.soe.ucsc.edu/admin/exe/))
    5. Generate gene-body meta-profile from bigWig files ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotProfile.html))
    6. Calculate genome-wide IP enrichment relative to control ([`deepTools`](https://deeptools.readthedocs.io/en/develop/content/tools/plotFingerprint.html))
    7. Calculate strand cross-correlation peak and ChIP-seq quality measures including NSC and RSC ([`phantompeakqualtools`](https://github.com/kundajelab/phantompeakqualtools))
6. Call broad/narrow peaks
    1. By ([`MACS2`](https://github.com/taoliu/MACS))
        1. Call Peaks
        2. Differentail binding analysis by shell script
            * Annotate peaks relative to gene features ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
            * Create consensus peakset across all samples and create tabular file to aid in the filtering of the data ([`BEDTools`](https://github.com/arq5x/bedtools2/))
            * Count reads in consensus peaks ([`featureCounts`](http://bioinf.wehi.edu.au/featureCounts/))
            * Differential binding analysis, PCA and clustering ([`R`](https://www.r-project.org/), [`DESeq2`](https://bioconductor.org/packages/DESeq2/))
        3. Differential binding analysis by [`DiffBind`](https://bioconductor.org/packages/DiffBind/)
    2. By ([`HOMER`](http://homer.ucsd.edu/homer/download.html))
        1. Call Peaks
        2. Differential binding analysis by [`DiffBind`](https://bioconductor.org/packages/DiffBind/)
        3. Annotate peaks relative to gene features ([`ChIPpeakAnno`](https://bioconductor.org/packages/ChIPpeakAnno/))
7. Visualisation the tracks.
    1. Create IGV session file containing bigWig tracks, peaks and differential sites for data visualisation ([`IGV`](https://software.broadinstitute.org/software/igv/)).
    2. Create UCSC genome browser track hub for bigWig tracks [trackhub](https://daler.github.io/trackhub/quickstart.html).
8. Present QC for raw read, alignment, peak-calling and differential binding results ([`MultiQC`](http://multiqc.info/), [`R`](https://www.r-project.org/))
9. Create index.html ([`R`](https://www.r-project.org/))

## Installation

```bash
# install nextflow
wget -qO- https://get.nextflow.io | bash
nextflow pull jianhong/chipseq -r dev
nextflow run jianhong/chipseq -profile test,conda -r dev
```

## Update

```bash
nextflow pull jianhong/chipseq -r dev
```

## Remove

```bash
nextflow drop jianhong/chipseq
```

## design table

To make sure DiffBind to be run, antibody must be provided.
The control could be empty. `track_` columes are optional.
save the design table as a csv file. See samples in assets/design*.

| group | replicate | fastq_1 | fastq_2 | antibody | control | track_color | track_group |
|-------|-----------|---------|---------|----------|---------|-------------|-------------|
| WT | 1 | fastq/WT1.fastq.gz| | ANT1 | Input | #E69F00 | SAMPLE |
| WT | 2 | fastq/WT2.fastq.gz| | ANT1 | Input | #E69F00 | SAMPLE |
| KD | 1 | fastq/KD1.fastq.gz| | ANT1 | Input | #0000FF | SAMPLE |
| KD | 2 | fastq/KD2.fastq.gz| | ANT1 | Input | #0000FF | SAMPLE |
| Input | 1 | fastq/KD1.fastq.gz| |  |  | #000000 | SAMPLE |
| Input | 2 | fastq/KD2.fastq.gz| |  |  | #000000 | SAMPLE |

## Metagene analysis

The pipeline can do metagene analysis for predefined bed files.

```bash
nextflow run jianhong/chipseq -profile test -resume --genomicElements beds/*.bed
```

## Example profile

Here is the example profile file named as sample.config for human samples.

```bash
params {
  config_profile_name = 'Full test profile'
  config_profile_description = 'Full test dataset to check pipeline function'

  // Input data
  input = 'path_to_design.csv'

  // Genome references
  genome = 'GRCh38'
}
```

And run the pipeline by:

```bash
nextflow run jianhong/chipseq -c sample.config --conda
```

Or by docker:

```bash
nextflow run jianhong/chipseq -c sample.config --docker
```

## Analysis data from GEO database

Here is the example to download data from GEO database and run analysis

```bash
nextflow run jianhong/chipseq --input GSE90661 --genome R64-1-1
```

## Change parameters for module setting and rerun the pipeline

First create a config file following this format:

```bash
params {
    // change conda software version.
    conda_softwares {
        samtools = "bioconda::samtools=1.09"
        trimgalore = "bioconda::cutadapt=1.18 bioconda::trim-galore=0.6.6"
    }
    // change module parameters, for example for homer findpeaks and macs2
    modules {
        'homer_findpeaks' {
            args       = ['h3k4me1': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k4me3': "-region -nfr",
                          'h3k9me1': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k9me2': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k9me3': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k14ac': "-region -size 1000 -minDist 2500 -C 0",
                          'h4k20me1': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k27me3': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k27ac': "-region -nfr",
                          'h3k36me3': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k79me2': "-region -size 1000 -minDist 2500 -C 0",
                          'h3k79me3': "-region -size 1000 -minDist 2500 -C 0",
                          'others': ""]
            publish_dir   = "homer"
        }
        'macs2_callpeak' {
            args          = "--keep-dup all"
            publish_dir   = "macs2"
        }
    }
}
```

To get more modules setting, please refer to
[modules.config file](https://raw.githubusercontent.com/jianhong/chipseq/master/conf/modules.config)

And then run the pipeline as following command:

```bash
nextflow run jianhong/chipseq -c sample.config -c path/to/your/module/config/file --conda -resume
```

## Get help

Please create an issue to submit your questions.

## Citation

See [citation](CITATIONS.md)
