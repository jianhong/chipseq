#!/usr/bin/env python

#######################################################################
#######################################################################
## Created on Nov. 10, 2020 to create UCSC trackhub file from file list
#######################################################################
#######################################################################

import os
import glob
import errno
import argparse
import trackhub

############################################
############################################
## PARSE ARGUMENTS
############################################
############################################

Description = 'Create UCSC trackhub file from a list of files and associated colours - ".bed", ".narrowPeak", ".broadPeak", ".bw", ".bigwig" files currently supported.'
Epilog = """Example usage: python create_trackhub.py <OUTPUT_FOLDER> <LIST_FILE> <GENOME>"""

argParser = argparse.ArgumentParser(description=Description, epilog=Epilog)

## REQUIRED PARAMETERS
argParser.add_argument('OUTPUT_FOLDER', help="Folder for UCSC trackhub files")
argParser.add_argument('LIST_FILE', help="Tab-delimited file containing two columns i.e. file_name\tcolour. Header isnt required.")
argParser.add_argument('GENOME', help="Full path to genome fasta file or shorthand for genome available in UCSC e.g. hg19.")
argParser.add_argument('EMAIL', help="email address")

## OPTIONAL PARAMETERS
argParser.add_argument('-pp', '--path_prefix', type=str, dest="PATH_PREFIX", default='', help="Path prefix to be added at beginning of all files in input list file.")
args = argParser.parse_args()

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
tracktypes = ['bigWig', 'bam', 'bigBed', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigNarrowPeak',
              'bigMaf', 'bigPsl', 'halSnake']
TrackType = {'bw':'bigWig', 'bb':'bigBed'}
Visibility = {'bw':'full', 'bb':'dense'}
for tt in tracktypes:
    TrackType[tt.lower()] = tt
    if tt in ['bam', 'bigBed', 'vcfTabix', 'bigNarrowPeak',
              'bigBarChart', 'bigChain', 'bigGenePred', 'bigNarrowPeak',
              'bigMaf', 'bigPsl', 'halSnake']:
      Visibility[tt.lower()] = 'dense'
    else:
      Visibility[tt.lower()] = 'full'

def create_trackhub(OutFolder,ListFile,Genome,EMAIL,PathPrefix=''):

    makedir(OutFolder)
    
    fileList = []
    fin = open(ListFile,'r')
    while True:
        line = fin.readline()
        if line:
            ifile,colour = line.strip().split('\t')
            if len(colour.strip()) == 0:
                colour = '0,0,178'
            fileList.append((PathPrefix.strip()+ifile,colour))
        else:
            break
            fin.close()

    # Initialize the components of a track hub, already connected together
    hub, genomes_file, genome, trackdb = trackhub.default_hub(
        hub_name="RNISRS_hub",
        short_label='Regeneromics Shared Resource hub',
        long_label='Regeneration Next Initiative Regeneromics Shared Resource hub',
        genome=Genome,
        email=EMAIL)
    
    for ifile,color in fileList:
        extension = os.path.splitext(ifile)[1].replace(".", "").lower()
        if extension in ['bed','broadpeak','narrowpeak']:
          pass
        elif extension in TrackType.keys():
          track = trackhub.Track(
            name=os.path.splitext(os.path.basename(ifile))[0],
            source=ifile,
            color=color,
            visibility=Visibility[extension],
            tracktype=TrackType[extension],
            autoScale='on')
          trackdb.add_tracks(track)
          linkname=os.path.join(OutFolder, Genome, os.path.splitext(os.path.basename(ifile))[0]+"."+TrackType[extension])
          makedir(os.path.join(OutFolder, Genome))
          os.symlink(ifile, linkname)
        else:
          pass
    
    hub.render(staging=OutFolder)

############################################
############################################
## RUN FUNCTION
############################################
############################################

create_trackhub(OutFolder=args.OUTPUT_FOLDER,ListFile=args.LIST_FILE,Genome=args.GENOME,EMAIL=args.EMAIL,PathPrefix=args.PATH_PREFIX)

############################################
############################################
############################################
############################################
