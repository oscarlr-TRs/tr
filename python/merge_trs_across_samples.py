#!/bin/env python
import sys

def read_assembly_genotypes(fn,merged_samples,sample):
    header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq",
              "hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score",
              "hap_1_perfect_motif_locus","hap_1_query_motif_locus",
              "hap_1_ref_left_flank","hap_1_query_left_flank",
              "hap_1_ref_right_flank","hap_1_query_right_flank", 
              "hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score",
              "hap_2_perfect_motif_locus","hap_2_query_motif_locus",
              "hap_2_ref_left_flank","hap_2_query_left_flank",
              "hap_2_ref_right_flank","hap_2_query_right_flank"]         
    with open(fn,'r') as fh:
        fn_header = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            tr = "_".join(line[0:6])
            if tr not in merged_samples:
                merged_samples[tr] = {}
            out = line[6:10] + line[16:20]
            merged_samples[tr][sample] = out
    return merged_samples

def read_subreads_genotypes(fn,merged_samples,sample):
    header = [ "chrom","start","end","motif_size","number_of_motifs","motif_seq",
               "hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score",
               "hap_1_avg_num_of_motifs","hap_1_max_num_of_motifs",
               "hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score",
               "hap_2_avg_num_of_motifs","hap_2_max_num_of_motifs"]
    with open(fn,'r') as fh:
        fn_header = fh.readline()
        for line in fh:
            line = line.rstrip().split('\t')
            tr = "_".join(line[0:6])
            if tr not in merged_samples:
                merged_samples[tr] = {}
            out = []
            if line[6] == "hap_1_number_of_motifs":
                continue
            if line[12] == "hap_2_number_of_motifs":
                continue
            if line[6] != ".":
                for i in range(6,10):                
                    sum_ = sum(map(float,line[i].split(",")))
                    length = len(line[i].split(","))
                    avg = sum_/float(length)
                    out.append(avg)
            else:
                out = line[6:10]
            if line[12] != ".":
                for i in range(12,16):
                    sum_ = sum(map(float,line[i].split(",")))
                    length = len(line[i].split(","))
                    avg = sum_/float(length)
                    out.append(avg)
            else:
                out += line[12:16]
            merged_samples[tr][sample] = out
    return merged_samples

def read_fofn(samples_fofn):
    assembly_sample_list = []
    subreads_sample_list = []
    merged_samples = {}
    with open(samples_fofn,'r') as samples_fofh:
        for line in samples_fofh:
            sample, source, fn = line.rstrip().split('\t')
            if source == "assembly":
                if sample not in assembly_sample_list:
                    assembly_sample_list.append(sample)
                merged_samples = read_assembly_genotypes(fn,merged_samples,sample)
            elif source == "subreads":
                if sample not in subreads_sample_list:
                    subreads_sample_list.append(sample)
                merged_samples = read_subreads_genotypes(fn,merged_samples,sample)
            else:
                sys.exit("Error 1. Need assembly or subreads genotypes")
    return (assembly_sample_list,subreads_sample_list,merged_samples)

def print_header(assembly_sample_list,subreads_sample_list):
    sample_list = assembly_sample_list + subreads_sample_list
    out_header = ["chrom","start","end","motif_size","number_of_motifs","motif_seq"]
    for sample in sample_list:
        for column in ["hap_1_number_of_motifs","hap_1_motif_score","hap_1_left_flank_score","hap_1_right_flank_score"]:
            out_header.append("%s_%s" % (sample,column))
        for column in ["hap_2_number_of_motifs","hap_2_motif_score","hap_2_left_flank_score","hap_2_right_flank_score"]:
            out_header.append("%s_%s" % (sample,column))
    out_header += ["number_of_samples","num_alleles","max_expansion","max_contraction","range","number_of_genotyped_haps"]
    print "\t".join(out_header)

def add_assembly_samples(assembly_sample_list,merged_samples,tr,alleles,out,tr_size):
    num_haps = 0
    for sample in assembly_sample_list:
        if sample not in merged_samples[tr]:
            out += [".",".",".","."]*2
        else:
            if merged_samples[tr][sample][0] == "hap_1_number_of_motifs":
                continue
            if merged_samples[tr][sample][0] != ".":
                alleles.add(float(merged_samples[tr][sample][0]))
                num_haps += 1
            if merged_samples[tr][sample][4] != ".":
                alleles.add(float(merged_samples[tr][sample][4]))                
                num_haps += 1
            out += merged_samples[tr][sample]
    return (out,alleles,num_haps)

def add_subreads_samples(subreads_sample_list,merged_samples,tr,alleles,out,tr_size):
    subread_threshold = 2
    num_haps = 0
    for sample in subreads_sample_list:
        if sample not in merged_samples[tr]:
            out += [".",".",".","."]*2
        else:
            if merged_samples[tr][sample][0] == "hap_1_number_of_motifs":
                continue
            if merged_samples[tr][sample][0] != ".":
                num_haps += 1
                hap_1_allele_size = float(merged_samples[tr][sample][0])
                add_new_allele = True
                if len(alleles) != 0:
                    if (hap_1_allele_size - max(alleles))/float(tr_size) < subread_threshold:
                        add_new_allele = False
                    elif abs(min(alleles)-hap_1_allele_size)/float(tr_size) < subread_threshold:
                        add_new_allele = False
                if add_new_allele:
                    alleles.add(float(merged_samples[tr][sample][0]))
            if merged_samples[tr][sample][4] != ".":
                num_haps += 1
                hap_2_allele_size = float(merged_samples[tr][sample][4])
                add_new_allele = True
                if len(alleles) != 0:
                    if (hap_2_allele_size - max(alleles))/float(tr_size) < subread_threshold:
                        add_new_allele = False
                    elif abs(min(alleles)-hap_2_allele_size)/float(tr_size) < subread_threshold:
                        add_new_allele = False
                if add_new_allele:
                    alleles.add(float(merged_samples[tr][sample][4]))                
            out += merged_samples[tr][sample]            
    return (out,alleles,num_haps)

def print_merged_genotypes(merged_samples,assembly_sample_list,subreads_sample_list):
    for tr in merged_samples:
        alleles = set()
        number_of_samples = len(assembly_sample_list) + len(subreads_sample_list)
        out = tr.split("_")
        if out[3] == "motif":
            continue
        tr_size = int(out[3])
        out,alleles,num_assembly_haps = add_assembly_samples(assembly_sample_list,merged_samples,tr,alleles,out,tr_size)
        out,alleles,num_subreads_haps = add_subreads_samples(subreads_sample_list,merged_samples,tr,alleles,out,tr_size)        
        num_haps = num_assembly_haps + num_subreads_haps
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
        out.append(num_haps)
        print "\t".join(map(str,out))

samples_fofn = sys.argv[1]
assembly_sample_list, subreads_sample_list, merged_samples = read_fofn(samples_fofn)
print_header(assembly_sample_list,subreads_sample_list)
print_merged_genotypes(merged_samples,assembly_sample_list,subreads_sample_list)
