#!/bin/env python
from bisect import bisect_left 
import pysam
import sys

def read_in_contig_coords(contig_coordsfn):
    coords = {}
    with open(contig_coordsfn,'r') as contig_coordsfh:
        header = contig_coordsfh.readline()
        for line in contig_coordsfh:
            line = line.rstrip().split('\t')
            ref = line[0]
            start = int(line[1])
            end = int(line[2])
            if ref not in coords:
                coords[ref] = {}
            coords[ref][(start,end)] = [None,None,None,None]
    return coords

def find_position(list_to_look_in,position_to_find,start_to_look_here):
    lb = 0
    if start_to_look_here != None:
        lb = start_to_look_here
    ub = len(list_to_look_in)
    while True:
        if lb == ub:
            return -1
        mid_index = (lb + ub) // 2
        item_at_mid = list_to_look_in[mid_index][0]
        if item_at_mid == position_to_find:
            return mid_index      
        if item_at_mid < position_to_find:
            lb = mid_index + 1    
        else:
            ub = mid_index        

def get_mapping_coords(mappings,contig_coords):
    samfile = pysam.AlignmentFile(mappings)
    for j,read in enumerate(samfile):
        if read.query_name not in contig_coords:
            continue
        chrom = samfile.get_reference_name(read.reference_id)
        ap = list(read.get_aligned_pairs())
        contig_start_ends = sorted(contig_coords[read.query_name].keys(),key=lambda x: x[0])        
        prev_start_position = None
        for i,(contig_start, contig_end) in enumerate(contig_start_ends): #contig_coords[read.query_name]):            
            ref_start = None
            ref_end = None
            query_start = None
            query_end = None
            exact = None
            start_search_position = find_position(ap,contig_start,prev_start_position)
            prev_start_position = start_search_position
            if start_search_position != -1:
                for query_pos, ref_pos in ap[start_search_position:]:
                    if query_pos == None:
                        continue
                    if ref_pos == None:
                        continue
                    if query_pos <= contig_start:
                        ref_start = ref_pos
                        query_start = query_pos
                    if contig_end < query_pos:
                        break
                    ref_end = ref_pos
                    query_end = query_pos
                    exact = False
                    if contig_start == query_start:
                        if contig_end == query_end:
                            exact = True
            contig_coords[read.query_name][(contig_start,contig_end)] = [chrom,ref_start,ref_end,exact]                    
    return contig_coords

contig_coords = read_in_contig_coords(sys.argv[1])
mapping = get_mapping_coords(sys.argv[2],contig_coords)

for contig_name in mapping:
    for contig_start,contig_end in mapping[contig_name]:
        out = [contig_name,contig_start,contig_end] + mapping[contig_name][(contig_start,contig_end)]
        print "\t".join(map(str,out))



