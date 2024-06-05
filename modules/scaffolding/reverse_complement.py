#!/usr/bin/env python.11
"""
    reverse_compliment.py
    30 May 2024

    Copyright (C) 2024-2025 Keiler Collier
    This is a script which takes a fasta file and a CSV of filtered alignments. 
    For each contig in the list, it reverse-compliments them and prints them back into a fasta file.
    It is designed to be used in conjunction with the outputs from get_contact_map.R, 
    and, more generally, to partially automate chromosomal curation.


    
    # RULES:
        If the CSV has a '-', you flip it.
        If two fragments are mapped to the same fragment, fuse them.

"""
import argparse
import itertools
import pandas as pd

#from Bio import SeqIO

#--------------------------------------------
# HELPER FUNCTIONS AND CLASSES:
#--------------------------------------------
class Alignment:
    """This is a class for holding subsets of alignment information from PAF files."""
    def __init__(self, qname = "query", tname = "target", length = 0, strand = None):
        self.qname = qname
        self.tname = tname
        self.length = length
        self.strand = strand
    
    def __str__(self):
      return f"Query name: {self.qname}\nTarget name: {self.tname}\nAlignment length: {self.length}\nPolarity: {self.strand}"
     
    def to_dict(self):
        return {
            'qname': self.qname,
            'tname': self.tname,
            'length': self.length,
            'strand': self.strand
            }
    
class Contig:
    """This is a class for holding contig data/names from fasta files. Yes, this is a worse reinvention of SeqIO, but it's also more portable."""
    def __init__(self, name = None, seq = None, length = 0):
        self.name = name
        self.seq = seq
        self.length = length

    def __str__(self):
      return f"Contig name: {self.name}\nContig length: {self.length}."
  
    def to_dict(self):
        return {
            'name': self.name,
            'seq': self.seq,
            'length': self.length
            }
    
    def reverseComplement(self):
        # reverse complements your string
        complement_dict = {"A": "T", "a": "t", # four standard nucleotides
                           "T": "A", "t": "a", 
                           "C": "G", "c": "g",
                           "G": "C", "g": "c",
                           # I don't think I'll have many ambiguous bases, but...
                           "R": "Y", "r": "y",
                           "Y": "R", "y": "r",
                           "K": "M", "k": "m",
                           "M": "K", "m": "k",
                           "B": "V", "b": "v",
                           "V": "B", "v": "b",
                           "D": "H", "d": "h",
                           "H": "D", "h": "d"
                           }
        complement_trans_tbl = self.seq.maketrans(complement_dict)
        self.seq = self.seq.translate(complement_trans_tbl)[::-1] # inverts and complements string
    
    def dumpSeq(self):
        # When used with print(), it dumps an entire formatted name and seq fasta combo to stdout.
        # This is useful for outputting the fasta.
        return f">{self.name}\n{self.seq}\n"
    
 # FIXME - find some way to deal with your stupid headers in your PAF output
def parse_paf(paf): 
    aln_lst= []
    # now, we read in lines and assign each to an Alignment obj. in aln_list
    with open(paf, 'r') as paf_obj:
        for line in paf_obj:
            aln = line.split()
            Aln = Alignment(aln[0], aln[5], int(aln[10]), aln[4])
            aln_lst.append(Aln)
            # We return the list. The file is closed as a condition of the while loop
        return(aln_lst)
         
# FIXME - find some way to deal with multi-line fastas
# Perl and seqkit oneliners to go from MULTI-TRACK DRIFTING to single track fasta
    #seqkit sort --by-length --reverse ${INPUT_FASTA} | seqkit replace --pattern '.+' --replacement 'Contig_{nr}' > ${OUTPUT_FASTA}
    #perl -pe '/^>/ ? print "\n" : chomp' Falco_peregrinus_best_genome_sorted_renamed.fa | tail -n +2 > new_falcon.fasta
def parse_fasta(fasta): 
    # initialize a list
    ctg_lst= []
    
    # determine if fasta is single or multiline
    
    
    
    # now, we read in lines and assign each to an Alignment obj. in aln_list
    with open(fasta, 'r') as fa_obj:
        for line, line2 in itertools.zip_longest(fa_obj, fa_obj, fillvalue = None):
            Ctg = Contig(name=line.strip().lstrip('>'), seq=line2.strip(), length=len(line2.strip()))
            ctg_lst.append(Ctg)
    return(ctg_lst)
 
def sortSeq_lst(seq):
    # Small function that provides a key to sort list objects
    # Works for Alignment and Contig objects because they both have lengths
    return seq.length
    
#--------------------------------------------
# MAIN:
#--------------------------------------------
def main():
    # PARSE ARGS:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", help = "Supply a .fa file.")
    parser.add_argument("-a", "--paf", help = "Supply a PAF file")
    parser.add_argument("-o", "--output", default = "stdout", help = "Writes output to path provided and suppresses output to stdout.")
    parser.add_argument("-p", "--prefix", default = "output", help = "Prefix for all output files.")
    args = parser.parse_args()

    
    #--------------------------------------------
    # parse a PAF and give us a list of Alignments in descending order.
    #--------------------------------------------
    aln_lst = parse_paf(args.paf)  
    aln_lst.sort(reverse=True, key=sortSeq_lst)
    
    
    #--------------------------------------------
    # parse a fasta and give us a list of Contigs  in descending order
    #--------------------------------------------
    ctg_lst = parse_fasta(args.fasta)
    ctg_lst.sort(reverse=True, key=sortSeq_lst)
    
     
    #--------------------------------------------
    # Iterate through the fasta. If the sequence name matches a sequence tagged as '-' in the PAF, reverse complement it.
    #--------------------------------------------
    # generate data frames for both 
    aln_df = pd.DataFrame.from_records(f.to_dict() for f in aln_lst) 
    ctg_df = pd.DataFrame.from_records(f.to_dict() for f in ctg_lst) #this wasn't actually used
    ctg_df["name"]
        #then subset to only strand = '-'
    # subset the alignment dataframe to only negative alignments over a given length - these are likely to be biologically signiciant
    min_aln_len = 500000
    ctgs_to_reverse = pd.unique(aln_df[(aln_df["strand"] == '-') & (aln_df["length"] >= min_aln_len)]["qname"])  # we could combine this all onto one line, but...
    
    
    # print fa. If the 'name' of ctg is in ctgs_to_reverse, reverse-complement it
    for ctg in ctg_lst:
        if ctg.name in ctgs_to_reverse:
            ctg.reverseComplement()
    
    
    #--------------------------------------------
    # Write the fasta object to a file.
    #--------------------------------------------
    # FIXME - figure out how to write this to a file
    
    if args.output == 'stdout':
        for contig in ctg_lst:
            print(contig.dumpSeq())# writes to stdout
            
        
if __name__ == "__main__":
    main()
        
        
