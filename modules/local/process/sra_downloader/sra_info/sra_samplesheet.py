#!/usr/bin/env python

import sys
import csv
import json
import re
import argparse

def parse_args(args=None):
    Description = "Create runinfo samplesheet."
    Epilog = 'Example usage: python get_runinfo.py <TAB_FILE> <JSON_FILE> <FILE_OUT>'

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument('TAB_FILE', help="Comma-separated list of metadata file created by efetch -format runinfo from sra database.")
    parser.add_argument('JSON_FILE', help="JSON of metadata file created by efetch --json from sra database.")
    parser.add_argument('INFO_OUT', help="Output file containing information for fasterq_dump.")
    parser.add_argument('DESIGN_OUT', help="Sample file for nextflow pipeline.")
    parser.add_argument('LIBRARY_FILTER', help="library strategy filter. eg. RNA-seq")
    return parser.parse_args(args)

def list_to_factor(l):
  f=[]
  u=[]
  for x in l:
    if x not in u:
      u.append(x)
  for x in l:
    f.append(u.index(x))
  return f

def get_sra_run_info(runTableFile, runJSONFile, outfile, designfile, libraryFilter):
  runheader = []
  runinfo = []
  runJSON = []
  
  with open(runTableFile, 'r') as f:
    reader = csv.reader(f, delimiter=',')
    runheader = next(reader)
    for row in reader:
      if row:
        runinfo.append(row)
  
  with open(runJSONFile, 'r') as f:
    runJSON = json.load(f)
  
  experimentInfo = []
  attrTAG = ["title"]
  if "title" not in runheader:
    runheader.append("title")
  
  for exp in runJSON["EXPERIMENT_PACKAGE_SET"]["EXPERIMENT_PACKAGE"]:
    accession = exp["SAMPLE"]["accession"]
    title = exp["SAMPLE"]["TITLE"]
    sample = {'accession': accession, 'title': title}
    attributes = exp["SAMPLE"]["SAMPLE_ATTRIBUTES"]["SAMPLE_ATTRIBUTE"]
    for attr in attributes:
      if attr["TAG"] not in runheader:
        runheader.append(attr["TAG"])
      if attr["TAG"] not in attrTAG:
        attrTAG.append(attr["TAG"])
      sample[attr["TAG"]] = attr["VALUE"]
    experimentInfo.append(sample)
  
  for run in runinfo:
    for i in range(len(runheader)-len(run)):
      run.append(None)
  
  #print(experimentInfo)
  #print(attrTAG)
  for exp in experimentInfo:
    for run in runinfo:
      if run[runheader.index("Sample")] == exp["accession"]:
        for item in attrTAG:
          if item in exp.keys() and item in runheader:
            run[runheader.index(item)] = exp[item]
          else:
            run[runheader.index(item)] = ""
        break
  
  selectheader = ["Run", "LibraryStrategy", "LibraryLayout", "ScientificName", 
                  "condition", "replicate"]
  unselectheader = ["Run", "ReleaseDate", "LoadDate", "spots", "bases", 
                    "spots_with_mates", 'avgLength', 'size_MB', 
                    'AssemblyName', 'download_path', 'Experiment', 'Sample',
                    'InsertSize', 'InsertDev', 'SRAStudy', 'BioProject', 
                    'Study_Pubmed_id', 'ProjectID', 'TaxID', 
                    'SampleName', 'g1k_pop_code', 'source', 'g1k_analysis_group', 
                    "Submission", "dbgap_study_accession", "Consent",
                    "BioSample",	"RunHash",	"ReadHash"]

  conditionA = []
  conditionB = []
  replicate = {}
  condition = {}

  for h in runheader:
    if h not in unselectheader:
      if h.lower() == "antibody":
        selectheader.append(h)
        next
      allInfo = []
      runID = []
      runCond = []
      for run in runinfo:
        curr = re.sub('\s+rep[licates ]*\d+$', '', run[runheader.index(h)], flags=re.I)
        runCond.append(curr)
        runID.append(run[runheader.index("Run")])
        if curr not in allInfo:
          allInfo.append(curr)
      if len(allInfo)>1 and len(allInfo)!=len(runCond):
        if conditionA:
          conditionB = list_to_factor(runCond)
          check = False
          for cond in conditionA:
            if cond != conditionB:
              check = True
              break
          if check:
            for i in range(len(conditionB)):
              condition[runID[i]] = condition[runID[i]] + "_" + "".join([c for c in runCond[i] if c.isalpha() or c.isdigit() or c=='_']).rstrip()
            conditionA.append(conditionB)
            selectheader.append(h)
        else:
          conditionA = [list_to_factor(runCond)]
          for i in range(len(runCond)):
              condition[runID[i]] = "".join([c for c in runCond[i] if c.isalpha() or c.isdigit() or c=='_']).rstrip()
          selectheader.append(h)
  
  replicateFlag = True
  for run in runinfo:
    x = re.search(r'rep\d+', run[runheader.index("title")].lower())
    if x:
      replicate[run[runheader.index("Run")]] = x.group().replace("rep", "")
      replicateFlag = False
    else:
      replicate[run[runheader.index("Run")]] = 1
  
  if replicateFlag:
    gp_id={}
    for key, value in condition.items():
      if value in gp_id.keys():
        gp_id[value] = gp_id[value]+1
      else:
        gp_id[value] = 1
      replicate[key] = gp_id[value]
  
  with open(outfile, mode='w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    writer.writerow(selectheader)
    
    for run in runinfo:
      if run[runheader.index('LibraryStrategy')].lower() == libraryFilter.lower():
        outinfo = []
        for h in selectheader:
          if h not in ["condition", "replicate"]:
            outinfo.append(run[runheader.index(h)])
          else:
            if h=="condition":
              outinfo.append(condition[run[runheader.index("Run")]])
            if h=="replicate":
              outinfo.append(replicate[run[runheader.index("Run")]])
        writer.writerow(outinfo)
  
  designTabColname = ['group','replicate','fastq_1','fastq_2', 'ScientificName']
  with open(designfile, mode='w') as f:
    writer = csv.writer(f, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    
    writer.writerow(designTabColname)
    
    for run in runinfo:
      if run[runheader.index('LibraryStrategy')].lower() == libraryFilter.lower():
        folder = "/".join(["fastq", run[runheader.index("ScientificName")],
                          run[runheader.index("LibraryLayout")],
                          run[runheader.index("LibraryStrategy")],
                          condition[run[runheader.index("Run")]],
                          str(replicate[run[runheader.index("Run")]])])
                  
        outinfo = [condition[run[runheader.index("Run")]], 
                   replicate[run[runheader.index("Run")]],
                   folder+"/"+run[runheader.index("Run")]+"_R1.fastq.gz",
                   None,
                   run[runheader.index("ScientificName")]]
        if run[runheader.index("LibraryLayout")].lower() == "paired":
          outinfo[designTabColname.index("fastq_2")] = folder+"/"+run[runheader.index("Run")]+"_R2.fastq.gz"
        writer.writerow(outinfo)

def main(args=None):
    args = parse_args(args)
    get_sra_run_info(args.TAB_FILE, args.JSON_FILE, args.INFO_OUT, args.DESIGN_OUT, args.LIBRARY_FILTER)
    

if __name__ == '__main__':
    sys.exit(main())
