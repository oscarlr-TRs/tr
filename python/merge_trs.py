#!/bin/env python
import sys

tandem_repeatfn = sys.argv[1]
hap1_fn = sys.argv[2]
hap2_fn = sys.argv[3]

output = {}

def score_alignment(ref,query):
    score = 0.0
    alignment = 0.0
    for r,q in zip(ref,query):
        if r.upper() == q.upper():
            score += 1.0
        alignment += 1.0
    score = score/alignment
    return score

def add_haplotype_info(hap_fn):
    with open(hap_fn, 'r') as hap_fh:
        for index,line in enumerate(hap_fh):
            line = line.rstrip().split('\t')
            id_1 = line[1]
            id_2 = line[2]
            motif_seq = line[3]
            number_of_motifs = line[4]
            perfect_motif_locus = line[5]
            query_motif_locus = line[6]
            ref_left_flank = line[7]
            query_left_flank = line[8]
            ref_right_flank = line[7]
            query_right_flank = line[8]

            assert id_1 == id_2

            motif_score = score_alignment(perfect_motif_locus,query_motif_locus)
            left_flank_score = score_alignment(ref_left_flank,query_left_flank)
            right_flank_score = score_alignment(ref_right_flank,query_right_flank)
            
            items_to_add = [number_of_motifs,motif_score,left_flank_score,right_flank_score,
                            perfect_motif_locus,query_motif_locus,
                            ref_left_flank,query_left_flank,
                            ref_right_flank,query_right_flank]

            motif_info = id_1.split("_")
            chrom = motif_info[0]
            start = motif_info[1]
            end = motif_info[2]
            motif_size = motif_info[3]
            number_of_motifs = motif_info[4]
            for item in items_to_add:
                output[(chrom,start,end,motif_size,number_of_motifs)].append(item)

with open(tandem_repeatfn,'r') as tandem_repeatfh:
    for line in tandem_repeatfh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = line[1]
        end = line[2]
        motif_size = line[3]
        number_of_motifs = line[4]
        motif_seq = line[6]
        output[(chrom,start,end,motif_size,number_of_motifs)] = [chrom,start,end,motif_size,number_of_motifs,motif_seq]

add_haplotype_info(hap1_fn)
add_haplotype_info(hap2_fn)

header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq",
          "hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score",
          "hap_1_perfect_motif_locus","hap_1_query_motif_locus",
          "hap_1_ref_left_flank","hap_1_query_left_flank",
          "hap_1_ref_right_flank","hap_1_query_right_flank", 
          "hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score",
          "hap_2_perfect_motif_locus","hap_2_query_motif_locus",
          "hap_2_ref_left_flank","hap_2_query_left_flank",
          "hap_2_ref_right_flank","hap_2_query_right_flank"]         

print "\t".join(header)
for index in output:
    if len(header) != len(output[index]):
        continue
    out_line = output[index]
    print "\t".join(map(str,out_line))


