#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on Nov. 10, 2020 to create UCSC trackhub file from file list
## Copyright (c) 2020 Jianhong Ou (jianhong.ou@gmail.com)
#######################################################################
#######################################################################

import os
import sys
import glob
import errno
import argparse
import trackhub

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

def parse_args(args=None):
    Description = 'Create UCSC trackhub file from a list of files and associated colours - ".bed", ".narrowPeak", ".broadPeak", ".bw", ".bigwig" files currently supported.'
    Epilog = """Example usage: python create_trackhub.py <OUTPUT_FOLDER> <LIST_FILE> <GENOME>"""

    argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    ## REQUIRED PARAMETERS
    argParser.add_argument('LIST_FILE', help="Tab-delimited file containing two columns i.e. samplename\tsignalfile. Header isnt required.")
    argParser.add_argument('GENOME', help="Full path to genome fasta file or shorthand for genome available in UCSC e.g. hg19.")
    argParser.add_argument('CHROM_SIZE', help="Full path to chrom size")
    argParser.add_argument('EMAIL', help="email address")
    argParser.add_argument('DESIGN_FILE', help="design file")

    return argParser.parse_args(args)


############################################
############################################
## HELPER FUNCTIONS
############################################
############################################

def makedir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

############################################
############################################
## MAIN FUNCTION
############################################
############################################
tracktypes = ['bigWig', 'bam', 'bigBed', 'bigInteract', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigBroadPeak',
              'bigMaf', 'bigPsl', 'halSnake']
TrackType = {'bw':'bigWig', 'bb':'bigBed', 'bed':'bigBed',
             'narrowpeak':'bigBed 6+4', 'broadpeak':'bigBed 6+3',
             'bedpe':'bigBed 5+13'}
Visibility = {'bw':'full', 'bb':'dense', 'bed':'dense', 'bedpe':'dense',
             'narrowpeak':'dense', 'broadpeak':'dense', 'bigInteract':'dense'}

bigNarrowPeak=open("narrowPeak.as", "w")
bigNarrowPeak.write('''table bigNarrowPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
    int   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."
)
''')
bigNarrowPeak.close()
bigBroadPeak=open("broadPeak.as", "w")
bigBroadPeak.write('''table bigBroadPeak
"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."
(
    string chrom;        "Reference sequence chromosome or scaffold"
    uint   chromStart;   "Start position in chromosome"
    uint   chromEnd;     "End position in chromosome"
    string name;	 "Name given to a region (preferably unique). Use . if no name is assigned"
    uint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "
    char[1]  strand;     "+ or - or . for unknown"
    float  signalValue;  "Measurement of average enrichment for the region"
    float  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."
    float  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."
)
''')
bigBroadPeak.close()
bigInteracion=open("interact.as", "w")
bigInteracion.write('''table interact
"interaction between two regions"
    (
    string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
    uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region"
    uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
    string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
    uint score;          "Score (0-1000)"
    double value;        "Strength of interaction or other data value. Typically basis for score"
    string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
    string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
    string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
    uint sourceStart;    "Start position in chromosome of source/lower/this region"
    uint sourceEnd;      "End position in chromosome of source/lower/this region"
    string sourceName;   "Identifier of source/lower/this region"
    string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
    string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
    uint targetStart;    "Start position in chromosome of target/upper/this region"
    uint targetEnd;      "End position in chromosome of target/upper/this region"
    string targetName;   "Identifier of target/upper/this region"
    string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"

    )
''')
bigInteracion.close()

for tt in tracktypes:
    TrackType[tt.lower()] = tt
    if tt in ['bam', 'bigBed', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigBroadPeak',
              'bigMaf', 'bigPsl', 'halSnake']:
      Visibility[tt.lower()] = 'dense'
    else:
      Visibility[tt.lower()] = 'full'

def create_trackhub(OutFolder,ListFile,Genome,ChrSize, EMAIL,DesignFile):
    ERROR_STR = 'ERROR: Please check design file'
    HEADER = ['group', 'replicate', 'fastq_1', 'fastq_2']

    makedir(OutFolder)

    dIn = open(DesignFile, 'r')
    header = dIn.readline().strip().split(',')
    if header[:4] != HEADER:
        print("{} header: {} != {}".format(ERROR_STR,','.join(header),','.join(HEADER)))
        sys.exit(1)

    paramColn = {}
    for i in range(len(header)):
      if header[i][:6]=="track_": # header start with track_
        paramColn[header[i][6:]]=i

    sampleDesignDict = {}
    designDict = {}
    if paramColn:
      while True:
        line = dIn.readline()
        if line:
          lspl = [x.strip() for x in line.strip().split(',')]
          lspl[0] = [lspl[0]+Postfix[1], lspl[0]+'_R'+lspl[1]]
          lspl[0] = [trackhub.helpers.sanitize(lspl[0][0].replace(".", "_"), strict=False),trackhub.helpers.sanitize(lspl[0][1].replace(".", "_"), strict=False)]
          sampleDesignDict[lspl[0][0]] = {}
          sampleDesignDict[lspl[0][1]] = {}
          for k in paramColn.keys():
            sampleDesignDict[lspl[0][0]][k]=lspl[paramColn[k]]
            sampleDesignDict[lspl[0][1]][k]=lspl[paramColn[k]]
            if k in designDict:
              designDict[k][lspl[paramColn[k]]] = lspl[paramColn[k]]
            else:
              designDict[k] = {lspl[paramColn[k]]:lspl[paramColn[k]]}
        else:
          break


    dIn.close()

    fileList = []
    fin = open(ListFile,'r')
    while True:
        line = fin.readline()
        if line:
            ifile = [x.strip() for x in line.strip().split('\t')]
            colour = ""
            if sampleDesignDict:
                kfile = trackhub.helpers.sanitize(ifile[0].replace(".", "_"), strict=False)
                if kfile in sampleDesignDict:
                  if "color" in sampleDesignDict[kfile]:
                    h = sampleDesignDict[kfile]["color"].lstrip('#')
                    colour = ','.join(str(x) for x in tuple(int(h[i:i+2], 16) for i in (0, 2, 4)))
            if len(colour.strip()) == 0:
              colour = '0,0,178'
            fileList.append((ifile[1],colour,ifile[0]))
        else:
            break
            fin.close()

    fileList = sorted(fileList, key=lambda x: x[2])

    # Initialize the components of a track hub, already connected together
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name="RNISRS_hub",
        short_label='Regeneromics Shared Resource hub',
        long_label='Regeneration Next Initiative Regeneromics Shared Resource hub',
        genome=Genome,
        email=EMAIL)

    # create compositeTracks
    if sampleDesignDict:
      composite = trackhub.CompositeTrack(
        name = 'composite',
        short_label='singlal'
      )
      # Add those subgroups to the composite track
      subgroups = []
      for k in designDict.keys():
        if k!='color':
          subg = trackhub.SubGroupDefinition(
            name=k,
            label=k,
            mapping=designDict[k]
          )
          subgroups.append(subg)

      composite.add_subgroups(subgroups)

      # Add the composite track to the trackDb
      trackdb.add_tracks(composite)
      signal_view = trackhub.ViewTrack(
          name='signalviewtrack',
          view='signal',
          short_label='Signal')
      composite.add_view(signal_view)
      regions_view = trackhub.ViewTrack(
          name='regionsviewtrack',
          view='regions',
          short_label='Regions')
      composite.add_view(regions_view)

    for ifile,color,id in fileList:
        extension = os.path.splitext(ifile)[1].replace(".", "").lower()
        filename = trackhub.helpers.sanitize(os.path.splitext(os.path.basename(ifile))[0].replace(".", "_"), strict=False)
        if extension in ['bed','broadpeak','narrowpeak','bedpe']:
          # convert bed to bigbed
          # sort bed file
          cmd = "sort -k1,1 -k2,2n " + ifile +" >" + ifile + "_sorted.bed"
          os.system(cmd)
          # bedToBigBed
          os.system("awk '$1 != \"track\" {$5=($5>1000)?1000:$5; print ;}' "+ifile+"_sorted.bed > "+ifile+"_srt.bed")
          if extension == "bed":
            cmd = "bedToBigBed "+ifile+"_srt.bed"+" "+ChrSize+" "+ifile+".bb"
            extension = "bb"
          if extension == "broadpeak":
            cmd = "bedToBigBed -type=bed6+3 -as=broadPeak.as "+ifile+"_srt.bed"+" "+ChrSize+" "+ifile+".bb"
          if extension == "narrowpeak":
            cmd = "bedToBigBed -type=bed6+4 -as=narrowPeak.as "+ifile+"_srt.bed"+" "+ChrSize+" "+ifile+".bb"
          if extension == "bedpe":
            cmd = "bedToBigBed -type=bed5+13 -as=interact.as "+ifile+"_srt.bed"+" "+ChrSize+" "+ifile+".bb"
          os.system(cmd)
          # change ifile to new bigbed file
          ifile = ifile+".bb"
        if extension in TrackType.keys():
          if sampleDesignDict:
            track = trackhub.Track(
              name=filename,
              source=filename,
              short_label=id,
              long_label=filename,
              color=color,
              visibility=Visibility[extension],
              tracktype=TrackType[extension],
              subgroups=sampleDesignDict[filename],
              autoScale='on')
            signal_view.add_tracks(track)
          else:
            track = trackhub.Track(
              name=filename,
              source=filename,
              short_label=id,
              long_label=filename,
              color=color,
              visibility=Visibility[extension],
              tracktype=TrackType[extension],
              autoScale='on')
            trackdb.add_tracks(track)
          linkname=os.path.join(OutFolder, Genome, filename+"."+TrackType[extension].split()[0])
          makedir(os.path.join(OutFolder, Genome))
          os.symlink("../../"+ifile, linkname)
        else:
          pass

    hub.render(staging=OutFolder)

############################################
############################################
## RUN FUNCTION
############################################
############################################
def main(args=None):
    args = parse_args(args)
    create_trackhub(
      OutFolder="trackhub",
      ListFile=args.LIST_FILE,
      Genome=args.GENOME,
      ChrSize=args.CHROM_SIZE,
      EMAIL=args.EMAIL,
      DesignFile=args.DESIGN_FILE)

if __name__ == "__main__":
    sys.exit(main())
