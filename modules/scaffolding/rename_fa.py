#!/usr/bin/env python3.11
"""
    reverse_compliment.py
    12 Jun 2024

    Copyright (C) 2024-2025 Keiler Collier
    This is a script which takes a fasta file (multiline or oneline) and a TSV.
    The TSV has old contig names on the left and new ones on the right, separated by a tab character.
    This script goes in and renames the contigs in the fasta to fit the new names.
    If a contig is present in the fasta but not in the TSV, it is ignored by default.
        However, if you enable the -r (--remove) flag, it removes those contigs from the assembly.
    
    It will also support pattern renaming - this is mostly useful for our uncurated assemblies.
    

"""
import argparse
import os
import pandas as pd

class Fasta:
    """"This is a class made to hold multiple contigs."""
    def __init__(self, contigs = [], species = None):
        self.contigs = contigs # this is a list of Contig objects
        self.species = species # this is metadata in case we want it
        self.length = len(contigs)
    
    def __str__(self):
      return f"Fasta species: {self.species}\nFasta length: {self.length}."
        
    def dump(self):
        # When used with print(), it dumps an entire formatted name and seq fasta combo to stdout.
        # This is useful for outputting the entire fasta without having to use a list (like with Contig class)
        for contig in self.contigs:
            print(contig.dump())


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
        return f">{self.name}\n{self.seq}"
    

def parse_fasta(fasta): 
    ctg_lst= [] # initialize a list 
    seq = name = ""
    
    with open(fasta, 'r') as f:
        for line in f:
            if line.startswith('>'):
                if seq == "":
                    name = line.strip().lstrip('>')
                else:
                    #initialize new contig and write to list
                    Ctg = Contig(name = name, seq = seq, length = len(seq))
                    ctg_lst.append(Ctg)
                    name = seq = "" #reset the loop for the next round
                    name = line.strip().lstrip('>')
            else:
                seq = seq + line.strip()
        Ctg = Contig(name = name, seq = seq, length = len(seq))
        ctg_lst.append(Ctg)
    return ctg_lst

    
    

#--------------------------------------------
# MAIN:
#--------------------------------------------    
def main():
    # PARSE ARGS:
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required = True, help = "Supply a .fa file")
    parser.add_argument("-n", "--names", default = None, help = "Supply a TSV file of contig names. Should be formatted as <old_name>\t<new_name>. Unnamed seqs will be renamed as per '-s'")
    parser.add_argument("-s", "--sequential", default = None, help = "Supply a prefix with which to rename your contigs. If no name is given via --names, seqs will be renamed as '<prefix>_1', '<prefix>_2'...")    
    parser.add_argument("-r", "--remove", default = False, action= 'store_true', help = "Remove contigs where a new name is not found, instead of ignoring them. False by default; incompatible with -s.")    
    parser.add_argument("-o", "--out_pre", default = "stdout", help = "Writes output to path provided. If unset, output written to stdout")
    args = parser.parse_args()
    # set default of 'sequential' to none. if a person sets the flag with a prefix, it sets a bool in the program to turn on seq. mode
    
    #--------------------------------------------
    # Sanity-check for fasta and TSV arguments
    #--------------------------------------------
    if (args.fasta == None or args.names == None):
        print("Required inputs not found:")
        if (args.fasta == None):
            print("\tPlease provide a fasta (-f or --fasta)")
        if (args.names == None):
            print("\tPlease provide a TSV (-n or --names)")
        return 1
    if (not os.path.exists(args.fasta) or not os.path.exists(args.names)):
        print("Required inputs not found:")
        if (not os.path.exists(args.fasta)):
            print("\t'" + args.fasta + "' does not exist")
        if (not os.path.exists(args.names)):
            print("\t'" + args.names + "' does not exist")
        return 2
    
    
    #--------------------------------------------
    # parse a fasta and give us a list of Contigs. DO NOT SORT.
    #--------------------------------------------
    ctg_lst = parse_fasta(args.fasta)
    #ctg_lst.sort(reverse=True, key=sortSeq_lst)
    
    #--------------------------------------------
    # parse a TSV and give us a list of old/new seqnames
    #--------------------------------------------
    names_df = pd.read_csv(args.names, sep = '\t', names = ["oldname", "newname"])
    names_df["oldname"]
    names_df["newname"]
    
    #--------------------------------------------
    # Replace contigs found in 'oldname' with 'newname' values
    #--------------------------------------------
    if args.sequential != None:
        unnamed_ctr = 0 # this is a counter for any unnamed contigs. Used when -s is enabled

    for i in ctg_lst:
        if i.name in list(names_df["oldname"]):
            # if it's in oldnames, take the corresponding 'oldname' and get its rownum
            idx = (list(names_df["oldname"]).index(i.name))
            newname = names_df["newname"].iloc[idx] # use the rownum to get the 'newname' value
            i.name = newname # rename with newname

        else:
            if (args.sequential != None):
                unnamed_ctr = unnamed_ctr + 1
                i.name = args.sequential + "_" + str(unnamed_ctr)
            elif (args.remove): #FIXME - find a way to get ride of contigs
                i.name = None

    #--------------------------------------------
    # Write the fasta object to a file (or stdout)
    #--------------------------------------------
    # Print to stdout
    if args.out_pre == 'stdout':
        for i in ctg_lst:       # fasta output
            if i.name == None and args.remove == True:
                del(i)
                continue
            print(i.dump())  
    # Print to files
    else:
        filename = args.out_pre + ".fa" # fasta output
        with open(filename, 'w') as f:
            for i in ctg_lst:
                if i.name == None and args.remove == True:
                    del(i)
                    continue
                print(i.dump(), file=f)
            f.close()
    return 0 # End of the main function
    
if __name__ == "__main__":
    main()
