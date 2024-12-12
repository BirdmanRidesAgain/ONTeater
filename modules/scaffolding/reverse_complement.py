#!/usr/bin/env python3.11
"""
    reverse_compliment.py
    30 May 2024

    Copyright (C) 2024-2025 Keiler Collier
    This is a script which takes a fasta file and a CSV of filtered alignments. 
    For each contig in the list, it reverse-compliments them and prints them back into a fasta file.
    It is designed to be used in conjunction with the outputs from get_contact_map.R, 
    and, more generally, to partially automate chromosomal curation.
    
    It will also output a new PAF file reflecting the contigs that have been flipped to negatives.


    
    # RULES:
        If the CSV has a '-', you flip it.
        If two fragments are mapped to the same fragment, fuse them.

"""
import os
import argparse
import itertools
import pandas as pd

#from Bio import SeqIO

#--------------------------------------------
# HELPER FUNCTIONS AND CLASSES:
#--------------------------------------------
# FIXME - Would be nice to intake the SAM values as well
        
class Alignment:
    """This is a class for holding all non-SAM alignment information from PAF files."""
    def __init__(self, line):
        # Initiating an alignment truncates all SAM key-value entries at the end. We only count up to twelve.
        aln = line.split()
        self.qname = aln[0]
        self.qlen = int(aln[1])
        self.qstart = int(aln[2])
        self.qend = int(aln[3])
        self.strand = aln[4]
        self.tname = aln[5]
        self.tlen = int(aln[6])
        self.tstart = int(aln[7])
        self.tend = int(aln[8])
        self.res_match = int(aln[9])
        self.length = int(aln[10])
        self.mapq = int(aln[11])
        
        self.SAM_tags = '\t'.join(i for i in aln[12:])
    
    def __str__(self):
      return f"Query name: {self.qname}\nTarget name: {self.tname}\nAlignment length: {self.length}\nPolarity: {self.strand}"
     
    def to_dict(self):
        return {
            'qname': self.qname,
            'qlen': self.qlen,
            'qstart': self.qstart,
            'strand': self.strand,
            'tname': self.tname,
            'tlen': self.tlen,
            'tstart': self.tstart,
            'tend': self.tend,
            'res_match': self.res_match,
            'length': self.length,
            'mapq': self.mapq
            }
    def dump(self):
        # When used with print(), it dumps a formatted alignment to args.output - this is our PAF file
        # It is okay for us to print this PAF, despite it not coming from an actual aligner, because we're not actually changing anything
        # The alignment polarities AND the reads are reversed - it's an equivalent alignment
        return f"{self.qname}\t{self.qlen}\t{self.qstart}\t{self.qend}\t{self.strand}\t{self.tname}\t{self.tlen}\t{self.tstart}\t{self.tend}\t{self.res_match}\t{self.length}\t{self.mapq}\t{self.SAM_tags}"
    
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
    
    def dump(self):
        # When used with print(), it dumps an entire formatted name and seq fasta combo to stdout.
        # This is useful for outputting the fasta.
        return f">{self.name}\n{self.seq}\n"
    
def parse_paf(paf): 
    aln_lst= []
    # now, we read in lines and assign each to an Alignment obj. in aln_list
    with open(paf, 'r') as paf_obj:
        for line in paf_obj:
            Aln = Alignment(line)
            aln_lst.append(Aln)
            # We return the list. The file is closed as a condition of the while loop
        return(aln_lst)
         
# FIXME - find some way to deal with multi-line fastas
# Perl and seqkit oneliners to go from MULTI-TRACK DRIFTING to single track fasta
    #seqkit sort --by-length --reverse ${INPUT_FASTA} | seqkit replace --pattern '.+' --replacement 'Contig_{nr}' | > ${OUTPUT_FASTA}
    #perl -pe '/^>/ ? print "\n" : chomp' ${INPUT_FASTA} | tail -n +2 > ${OUTPUT_FASTA}
def parse_fasta(fasta): 
    # initialize a list
    ctg_lst= []    
    
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
    parser.add_argument("-f", "--fasta", default = None, help = "Supply a .fa file")
    parser.add_argument("-a", "--paf", default = None, help = "Supply a PAF file")
    parser.add_argument("-m", "--min_len", default = 1000, help = "Minimum length for alignments to be reverse-complemented.\nDefault is 1000 - should be set shorter for fragmented or taxonomically distant alignments")
    parser.add_argument("-o", "--out_pre", default = "stdout", help = "Writes output to path provided. If unset, output written to stdout")
    args = parser.parse_args()
    
    #--------------------------------------------
    # Sanity-check for fasta and PAF arguments
    #--------------------------------------------

    if (args.fasta == None or args.paf == None):
        print("Required inputs not found:")
        if (args.fasta == None):
            print("\tPlease provide a fasta (-f or --fasta)")
        if (args.paf == None):
            print("\tPlease provide a PAF (-a or --paf)")
        return 1
    if (not os.path.exists(args.fasta) or not os.path.exists(args.paf)):
        print("Required inputs not found:")
        if (not os.path.exists(args.fasta)):
            print("\t'" + args.fasta + "' does not exist")
        if (not os.path.exists(args.paf)):
            print("\t'" + args.paf + "' does not exist")
        return 2

    
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
    #ctg_df = pd.DataFrame.from_records(f.to_dict() for f in ctg_lst) #this wasn't actually used

        #then subset to only strand = '-'
    # subset the alignment dataframe to only negative alignments over a given length - these are likely to be biologically signiciant
    ctgs_to_reverse = pd.unique(aln_df[(aln_df["strand"] == '-') & (aln_df["length"] >= args.min_len)]["qname"])  # we could combine this all onto one line, but...
    
    
    
    
    #--------------------------------------------
    # Write the fasta object to a file (or stout). If the contig is a '-' alignment, reverse complement it
    #--------------------------------------------
    # Print to stdout
    if args.out_pre == 'stdout':
        for i in ctg_lst:       # fasta output
            if i.name in ctgs_to_reverse:
                i.reverseComplement()
                i.name = i.name + "_FLIP" # test this out - not sure how it'll interact with the 'ctgs_to_reverse' check in line 214
            print(i.dump())
            
        for i in aln_lst:       # PAF output
            if i.strand == '-':
                i.strand = '+'
                i.qname = i.qname + "_FLIP"
            print(i.dump())
            
    # Print to files
    else:
        filename = args.out_pre + ".fa" # fasta output
        with open(filename, 'w') as f:
            for i in ctg_lst:
                if i.name in ctgs_to_reverse:
                    i.reverseComplement()
                    i.name = i.name + "_FLIP"
                print(i.dump(), file=f)
            f.close()
            
        filename = args.out_pre + ".paf" # PAF output
        with open(filename, 'w') as f:
            for i in aln_lst:
                if i.strand == '-':
                    i.strand = '+'
                    i.qname = i.qname + "_FLIP"
                print(i.dump(), file=f)
            f.close()
            


    return 0 # End of the main function
        
if __name__ == "__main__":
    main()
        
        
