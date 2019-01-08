#!/bin/env python
import sys
import gzip

vcffn = sys.argv[1]

def parse_info(info,item):
    length = len(item)
    index_start = info.index(item) + length
    if ";" in info[index_start:]: 
        index_end = info.index(";",index_start)
    else:
        index_end = len(info)
    parsed_info = info[index_start:index_end]
    return parsed_info

with gzip.open(vcffn,'rb') as vcffh:
    for line in vcffh:
        if line[0] == "#":
            continue
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = line[1]
        sv_id = line[2]
        ref_seq = line[3]
        sv_seq = line[4]
        qual = line[6]
        info = line[7]
        genotype = line[9].split(":")[0]
        if genotype in ["./.","0/0"]:
            continue
        sv_type = parse_info(info,"SVTYPE=")
        end = parse_info(info,"END=")
        sv_len = None
        if "SVLEN" in info:
            sv_len = parse_info(info,"SVLEN=")
        if sv_len != None:
            if "-" in sv_len:
                sv_len = sv_len[1:]
        out = [chrom,start,end,sv_type,genotype,sv_len,ref_seq,sv_seq,sv_id,qual]
        print "\t".join(map(str,out))

        
        
