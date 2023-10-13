#----------------------------------------------------------------------------------------
# Objective:
#	Calculate average coverage from alignment SAM file
#
# Input:
#   SAM file
#
# Output:
#   file with the following columns:
#       genome | coverage

# Author:
#	Vitalii Stebliankin (vsteb002@fiu.edu)
#           Florida International University
#           Bioinformatics Research Group
#----------------------------------------------------------------------------------------

import pandas as pd
import os
import argparse

#----------------------------------------------------------------------------------------
# Functions
#----------------------------------------------------------------------------------------

def average_coverage_unit(genome_length, read_length):
    # Add average coverage based on single hit
    # C = NL / G
    #     where:
    #     C is the average coverage.
    #     N is the total number of reads that align to "this" sequence.
    #     L is the length of the read.
    #     G is the sequence length.
    return (int(read_length))/int(genome_length)

import PyIO
import PyPluMA
class AbundanceGenomesPlugin:
    def input(self, inputfile):
        self.parameters=PyIO.readParameters(inputfile)
    def run(self):
        pass
    def output(self, outputfile):

      SAM = PyPluMA.prefix()+"/"+self.parameters["SAM"]
      COV_CUTOFF = float(self.parameters["COV_CUTOFF"])
      GR_COV = float(self.parameters["GR_COV"])
      OUTPUT_ABUNDANCE = outputfile+".abundance.txt"
      OUTPUT_LIST = outputfile+".genes.txt"

      with open(SAM, 'r') as f:
       coverage_dict = {}
       length_dict = {}
       for row in f.readlines():
        if ("@HD" in row) or ("@PG" in row):
            pass
        elif "@SQ" in row:
            row = row.strip("\n")
            row = row.split("\t")
            genome = row[1][3:]
            length = row[2].split(":")[1]
            length_dict[genome] = int(length)
            coverage_dict[genome] = 0
        else:
            row = row.strip("\n")
            row = row.split("\t")
            genome = row[2]
            READ_LENGTH = len(row[9])
            curr_coverage = coverage_dict[genome] + average_coverage_unit(length_dict[genome], READ_LENGTH)

            coverage_dict[genome] = curr_coverage
      coverage_df = pd.DataFrame.from_dict(coverage_dict, orient="index")
      coverage_df["genome"] = coverage_df.index
      coverage_df = coverage_df.rename(index=str, columns={0: "coverage"}).reset_index(drop=True)
      coverage_df = coverage_df[coverage_df["coverage"]>COV_CUTOFF]

      coverage_df = coverage_df.sort_values(by="coverage", ascending=False)
      coverage_df[["genome", "coverage"]].to_csv(OUTPUT_ABUNDANCE, sep="\t", index=False)
      # Getting the list of genomes that pass the coverage cutoff

      gr_df = pd.DataFrame()
      gr_df["genome"] = coverage_df[coverage_df["coverage"]> GR_COV]["genome"]
      gr_df["genome"] = gr_df["genome"].apply(lambda x: x.split("_")[1] + "_" + x.split("_")[2] + "_" + x.split("_")[4] )
      gr_df = gr_df["genome"]
      if len(gr_df)>0:
        gr_df.to_csv(OUTPUT_LIST, index=False)
      pass
