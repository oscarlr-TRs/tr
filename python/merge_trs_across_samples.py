#!/bin/env python
import sys

samples_fofn = sys.argv[1]

header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq",
          "hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score",
          "hap_1_perfect_motif_locus","hap_1_query_motif_locus",
          "hap_1_ref_left_flank","hap_1_query_left_flank",
          "hap_1_ref_right_flank","hap_1_query_right_flank", 
          "hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score",
          "hap_2_perfect_motif_locus","hap_2_query_motif_locus",
          "hap_2_ref_left_flank","hap_2_query_left_flank",
          "hap_2_ref_right_flank","hap_2_query_right_flank"]         



def read_fofn(samples_fofn):
    sample_list = []
    merged_samples = {}
    with open(samples_fofn,'r') as samples_fofh:
        for line in samples_fofh:
            sample, fn = line.rstrip().split('\t')
            if sample not in sample_list:
                sample_list.append(sample)
            with open(fn,'r') as fh:
                fn_header = fh.readline()
                for line in fh:
                    line = line.rstrip().split('\t')
                    tr = "_".join(line[0:6])
                    if tr not in merged_samples:
                        merged_samples[tr] = {}
                    out = line[6:10] + line[16:20]
                    merged_samples[tr][sample] = out
    return (sample_list,merged_samples)

out_header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq"]
for sample in sample_list:
    for column in ["hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score"]:
        out_header.append("%s_%s" % (sample,column))
    for column in ["hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score"]:
        out_header.append("%s_%s" % (sample,column))
out_header += ["number_of_samples","num_alleles","max_expansion","max_contraction","range"]
print "\t".join(out_header)

for tr in merged_samples:
    alleles = set()
    number_of_samples = 0
    out = tr.split("_")
    for sample in sample_list:
        if sample not in merged_samples[tr]:
            out += [".",".",".","."]*2
        else:
            if merged_samples[tr][sample][0] == "hap_1_number_of_motifs":
                continue
            number_of_samples += 1
            if merged_samples[tr][sample][0] != ".":
                alleles.add(float(merged_samples[tr][sample][0]))
            if merged_samples[tr][sample][4] != ".":
                alleles.add(float(merged_samples[tr][sample][4]))                
            out += merged_samples[tr][sample]
    if len(alleles) == 0:
        continue
    if max(alleles) > float(out[4]):
        max_expansion = max(alleles) - float(out[4])
    else:
        max_expansion = 0
    if min(alleles) < float(out[4]):
        max_contraction = float(out[4]) - min(alleles)
    else:
        max_contraction = 0
    range_ = max(alleles) - min(alleles)
    num_alleles = len(alleles)
    out.append(number_of_samples)
    out.append(num_alleles)
    out.append(max_expansion)
    out.append(max_contraction)
    out.append(range_)
    print "\t".join(map(str,out))
