#!/bin/env python
import sys

def read_in_contig_coords(contig_coordsfn):
    coords = {}
    with open(contig_coordsfn,'r') as contig_coordsfh:
        for line in contig_coordsfh:
            line = line.rstrip().split('\t')
            ref = line[0]
            start = line[1]
            end = line[2]
            if ref not in coords:
                coords[ref] = []
            coords[ref].append([start,end,None,None,None,None])
    for ref in coords:
        coords[ref] = sorted(coords[ref],key=lambda x: x[0])
    return coords

def get_mapping_coords(mappings,contig_coords):
    samfile = pysam.AlignmentFile(mappings)
    for read in samfile:
        if read.query_name not in contig_coords:
            continue
        contig_start = contig_coords[read.query_name][0]
        contig_end = contig_coords[read.query_name][1]
        ref_start = None
        ref_end = None
        query_start = None
        query_end = None
        for query_pos, ref_pos in read.get_aligned_pairs():
            if query_pos == None:
                continue
            if ref_pos == None:
                continue
            if query_pos <= contig_start:
                ref_start = ref_pos
                query_start = query_pos
            if contig_end > query_pos:
                break
            ref_end = ref_pos
            query_end = query_pos
        exact = False
        if contig_start == query_start:
            if contig_end == query_end:
                exact = True
        contig_coords[read.query_name][2] = samfile.get_reference_name(read.reference_id)
        contig_coords[read.query_name][3] = ref_start
        contig_coords[read.query_name][4] = ref_end
        contig_coords[read.query_name][5] = exact
    return contig_coords

contig_coords = read_in_contig_coords(sys.argv[1])
mapping = get_mapping_coords(sys.argv[2],contig_coords)

for contig_name in mapping:
    for mapped_coords in mapping[contig_name]:
        out = [contig_name] + mapped_coords
        print "\t".join(out)


