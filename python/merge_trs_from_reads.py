#!/bin/env python
import sys

tandem_repeatfn = sys.argv[1]
hap0_fn = sys.argv[2]
hap1_fn = sys.argv[3]
hap2_fn = sys.argv[4]

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

def add_haplotype_info(hap_fn,hap_key):
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
            if len(line) != 9:
                ref_right_flank = line[9]
                query_right_flank = line[10]
                right_flank_score = score_alignment(ref_right_flank,query_right_flank)
            else:
                right_flank_score = 0
            motif_score = score_alignment(perfect_motif_locus,query_motif_locus)
            left_flank_score = score_alignment(ref_left_flank,query_left_flank)            
            items_to_add = [number_of_motifs,motif_score,left_flank_score,right_flank_score]
            motif_info = id_1.split("_")
            chrom = motif_info[0]
            start = motif_info[1]
            end = motif_info[2]
            output[(chrom,start,end,motif_seq)][hap_key].append(items_to_add)

def output_hap_line(hap_data):
    hap_number_of_motifs_out = []
    hap_motif_score_out = []
    hap_left_flank_score_out = []
    hap_right_flank_score_out = []
    for hap_number_of_motifs,hap_motif_score,hap_left_flank_score,hap_right_flank_score in hap_data:
        hap_number_of_motifs_out.append(hap_number_of_motifs)
        hap_motif_score_out.append(hap_motif_score)
        hap_left_flank_score_out.append(hap_left_flank_score)
        hap_right_flank_score_out.append(hap_right_flank_score)
    if len(hap_number_of_motifs_out) == 0:
        hap_number_of_motifs_out = ["."]
        hap_motif_score_out = ["."]
        hap_left_flank_score_out = ["."]
        hap_right_flank_score_out = ["."]
        hap_number_of_motifs_average = "."
        hap_number_of_motifs_max = "."
    else:
        hap_number_of_motifs_average = sum(map(float,hap_number_of_motifs_out))/float(len(hap_number_of_motifs_out))
        hap_number_of_motifs_max = max(map(float,hap_number_of_motifs_out))
    hap_number_of_motifs_out = ",".join(map(str,hap_number_of_motifs_out))
    hap_motif_score_out = ",".join(map(str,hap_motif_score_out))
    hap_left_flank_score_out = ",".join(map(str,hap_left_flank_score_out))
    hap_right_flank_score_out = ",".join(map(str,hap_right_flank_score_out))
    out = [hap_number_of_motifs_out,hap_motif_score_out,hap_left_flank_score_out,hap_right_flank_score_out,hap_number_of_motifs_average,hap_number_of_motifs_max]
    return out
    
with open(tandem_repeatfn,'r') as tandem_repeatfh:
    for line in tandem_repeatfh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = line[1]
        end = line[2]
        motif_size = line[3]
        number_of_motifs = line[4]
        motif_seq = line[6]
        if "/" in motif_seq:
            motif_seq = motif_seq.split("/")[0]
        #index = str(index)
        output[(chrom,start,end,motif_seq)] = {}
        output[(chrom,start,end,motif_seq)]["TR"] = [chrom,start,end,motif_size,number_of_motifs,motif_seq]
        output[(chrom,start,end,motif_seq)]["hap0"] = []
        output[(chrom,start,end,motif_seq)]["hap1"] = []
        output[(chrom,start,end,motif_seq)]["hap2"] = []

add_haplotype_info(hap0_fn,"hap0")
add_haplotype_info(hap1_fn,"hap1")
add_haplotype_info(hap2_fn,"hap2")

header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq",          
          "hap_0_number_of_motifs","hap_0_motif_score","hap_0_left_flank_score","hap_0_right_flank_score","hap_0_avg_num_of_motifs","hap_0_max_num_of_motifs",
          "hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score","hap_1_avg_num_of_motifs","hap_1_max_num_of_motifs",
          "hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score","hap_2_avg_num_of_motifs","hap_2_max_num_of_motifs"]

print "\t".join(header)
for chrom,start,end,motif_seq in output:
    tr = output[(chrom,start,end,motif_seq)]["TR"]
    hap0 = output_hap_line(output[(chrom,start,end,motif_seq)]["hap0"])
    hap1 = output_hap_line(output[(chrom,start,end,motif_seq)]["hap1"])
    hap2 = output_hap_line(output[(chrom,start,end,motif_seq)]["hap2"])        
    out_line = tr + hap0 + hap1 + hap2
    print "\t".join(map(str,out_line))


