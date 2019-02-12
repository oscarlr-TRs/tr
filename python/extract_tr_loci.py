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
    sequence = []
    num_contigs = 0
    query_start = None
    query_end = None
    contig_sequence = None
    for contig in samfile.fetch(chrom,start,end):
        num_contigs += 1
        ap = contig.get_aligned_pairs()
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
        contig_sequence = contig.query_sequence[query_start:query_end]
    if num_contigs > 1:
        contig_sequence = None
    return contig_sequence

with open(simple_repeatsbed,'r') as fh:
    header = fh.readline()
    for i,line in enumerate(fh):
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        tr_seq = line[6]
        
        contig_seq = get_sequence_from(chrom,start-flank,end+flank)
        if contig_seq == None:
            continue
        i = str(i)
        query_record = SeqRecord(Seq(contig_seq),id=i,name=i,description="")
        query_records.append(query_record)

        tr_record = SeqRecord(Seq(tr_seq),id=i,name=i,description="")
        tr_records.append(tr_record)
        
        prefix_seq = ref.fetch(chrom,start-500,start)
        prefix_record = SeqRecord(Seq(prefix_seq),id=i,name=i,description="")
        prefix_records.append(prefix_record)

        suffix_seq = ref.fetch(chrom,end,end+500)
        suffix_record = SeqRecord(Seq(suffix_seq),id=i,name=i,description="")
        suffix_records.append(suffix_record)

SeqIO.write(tr_records, tr_fasta, "fasta")
SeqIO.write(prefix_records, prefix_fasta, "fasta")
SeqIO.write(suffix_records, suffix_fasta, "fasta")
SeqIO.write(query_records, query_fasta, "fasta")
