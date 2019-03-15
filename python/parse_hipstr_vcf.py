#!/bin/env python
import sys
import gzip

vcffn = sys.argv[1]

def get_data_from_line(line,dat):
    length = len(dat)
    index_start = line.index(dat) + length
    index_end = line.index(";",index_start)
    return line[index_start:index_end]

def get_period(line):
    return int(get_data_from_line(line,"PERIOD="))

def get_end(line):
    return get_data_from_line(line,"END=")

def get_genotype(line):
    index_end = line.index(":")
    return line[:index_end]
    
def get_motif_counts(ref_seq,alt_seq,period,genotype):
    motif_counts = []
    seqs = [ref_seq] + alt_seq.split(",")
    indices = map(int,genotype.split("|"))
    for index in indices:
        motif_counts.append(len(seqs[index])/period)
    return motif_counts

with open(vcffn,'rb') as vcffh:
    for line in vcffh:
        if line[0] == "#":
            continue
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = line[1]
        tr_id = line[2]
        ref_seq = line[3]
        alt_seq = line[4]
        period = get_period(line[-3])
        genotype = get_genotype(line[-1])
        end = get_end(line[-3])
        if genotype in ["./.","0/0","0|0"]:
            continue
        motif_counts = get_motif_counts(ref_seq,alt_seq,period,genotype)
        out = [chrom,start,end] + motif_counts
        print "\t".join(map(str,out))
