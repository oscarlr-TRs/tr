#!/bin/env python
import sys
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

bamfile = sys.argv[1]
simple_repeatsbed = sys.argv[2]
reffn = sys.argv[3]
output = sys.argv[4]

tr_fasta = "%s/tr.fasta" % output
prefix_fasta = "%s/prefix.fasta" % output
suffix_fasta = "%s/suffix.fasta" % output
query_fasta = "%s/query.fasta" % output

tr_records = []
prefix_records = []
suffix_records = []
query_records = []

ref = pysam.FastaFile(reffn)
samfile = pysam.AlignmentFile(bamfile)

flank = 500

def get_sequence_from(chrom,start,end):    
    sequences = []
    query_start = None
    query_end = None
    read_sequence = None
    for read in samfile.fetch(chrom,start,end):
        num_reads += 1
        ap = read.get_aligned_pairs()
        for q_pos, r_pos in ap:
            if q_pos == None:
                continue
            if r_pos == None:
                continue
            if r_pos <= start:
                query_start = q_pos
            query_end = q_pos
            if r_pos > end:
                break
        read_sequence = read.query_sequence[query_start:query_end]
        sequences.append((read_sequence,read.query_name))
    return sequences

with open(simple_repeatsbed,'r') as fh:
    header = fh.readline()
    for line in fh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        motif_length = line[3]
        num_of_motif = line[4]
        tr_seq = line[6]
        
        if "/" in tr_seq:
            tr_seq = tr_seq.split("/")[0]

        read_seqs = get_sequence_from(chrom,start-flank,end+flank)
        if len(read_seqs) == 0:
            continue
        
        for read_seq, read_name in read_seqs:
            out_name = "%s_%s_%s_%s_%s_%s" % (chrom,start,end,motif_length,num_of_motif,read_name)
            
            query_record = SeqRecord(Seq(read_seq),id=out_name,name=out_name,description="")
            query_records.append(query_record)
            
            tr_record = SeqRecord(Seq(tr_seq),id=out_name,name=out_name,description="")
            tr_records.append(tr_record)
            
            prefix_seq = ref.fetch(chrom,start-500,start)
            prefix_record = SeqRecord(Seq(prefix_seq),id=out_name,name=out_name,description="")
            prefix_records.append(prefix_record)
            
            suffix_seq = ref.fetch(chrom,end,end+500)
            suffix_record = SeqRecord(Seq(suffix_seq),id=out_name,name=out_name,description="")
            suffix_records.append(suffix_record)

SeqIO.write(tr_records, tr_fasta, "fasta")
SeqIO.write(prefix_records, prefix_fasta, "fasta")
SeqIO.write(suffix_records, suffix_fasta, "fasta")
SeqIO.write(query_records, query_fasta, "fasta")
