#!/bin/env python
import sys
import vcf
from collections import namedtuple

Number_of_motifs = namedtuple(
    "chrom",
    "start",
    "end",
    "hap_1_number_of_motifs",
    "hap_2_number_of_motifs")

def read_trs_file(filename):
    trs = {}
    with open(filename,'r') as fh:
        for line in fh:
            line = line.rstrip().split('\t')
            site = Number_of_motifs._make(line)
            motif_count = [site.hap_1_number_of_motifs,site.hap_2_number_of_motifs]
            trs[(site.chrom,site.start,site.end)] = motif_count
    return trs

def number_of_motifs_off(f1_h1_num,f1_h2_num,f2_h1_num,f2_h2_num):    
    min_score = min(abs(f1_h1_num-f2_h1_num) + abs(f1_h2_num-f2_h2_num),
                    abs(f1_h1_num-f2_h2_num) + abs(f1_h2_num-f2_h1_num))
    return min_score

def compare_trs(file1_trs,file2_trs):
    output = {}
    for chrom,start,end in file1_trs:
        f1_h1_num = file1_trs[(chrom,start,end)][0]
        f1_h2_num = file1_trs[(chrom,start,end)][1]
        if (chrom,start,end) not in file2_trs:
            f2_h1_num = "."
            f2_h2_num = "."
            diff = "."
        else:
            f2_h1_num = file2_trs[(chrom,start,end)][0]
            f2_h2_num = file2_trs[(chrom,start,end)][1]
            diff = number_of_motifs_off(f1_h1_num,f1_h2_num,f2_h1_num,f2_h2_num)
        output[(chrom,start,end)] = [f1_h1_num,f1_h2_num,f2_h1_num,f2_h2_num,diff]
    for chrom,start,end in file2_trs:
        if (chrom,start,end) in output:
            continue        
        assert (chrom,start,end) not in file1_trs        
        f2_h1_num = file2_trs[(chrom,start,end)][0]
        f2_h2_num = file2_trs[(chrom,start,end)][1]
        f1_h1_num = "."
        f1_h2_num = "."
        diff = "."
        output[(chrom,start,end)] = [f1_h1_num,f1_h2_num,f2_h1_num,f2_h2_num,diff]
    return output

            
file1_trs = read_trs_file(sys.argv[1])
file2_trs = read_trs_file(sys.argv[2])
comparisons = compare_trs(file1_trs,file2_trs)

header = ["chrom","start","end","diff",
          "file_1_hap1_motif_count","file_1_hap2_motif_count",
          "file_2_hap1_motif_count","file_2_hap2_motif_count"]

print "\t".join(header)
for chrom,start,end in comparisons:
    line = [chrom,start,end] + comparisons[(chrom,start,end)]
    print "\t".join(line)

