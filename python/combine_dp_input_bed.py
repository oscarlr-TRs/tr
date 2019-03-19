#!/bin/env python
import sys

tandem_repeatfn = sys.argv[1]
in_fn = sys.argv[2]

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
            for item in items_to_add:
                output[(chrom,start,end)].append(item)

with open(tandem_repeatfn,'r') as tandem_repeatfh:
    for line in tandem_repeatfh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = line[1]
        end = line[2]
        motif_size = line[3]
        number_of_motifs = line[4]
        motif_seq = line[6]
        #index = str(index)
        output[(chrom,start,end)] = [chrom,start,end,motif_size,number_of_motifs,motif_seq]

add_haplotype_info(in_fn)

header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq",
          "in_fn_number_of_motifs","in_fn_motif_score","in_fn_left_flank_score","in_fn_right_flank_score",
          "in_fn_perfect_motif_locus","in_fn_query_motif_locus",
          "in_fn_ref_left_flank","in_fn_query_left_flank",
          "in_fn_ref_right_flank","in_fn_query_right_flank"]

print "\t".join(header)
for index in output:
    if len(header) != len(output[index]):
        continue
    out_line = output[index]
    print "\t".join(map(str,out_line))


