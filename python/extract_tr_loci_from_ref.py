#!/bin/env python
import sys

simple_repeatsbed = sys.argv[1]
fastafile = sys.argv[2]
output = sys.argv[3]

tr_fasta = "%s/tr.fasta" % output
prefix_fasta = "%s/prefix.fasta" % output
suffix_fasta = "%s/suffix.fasta" % output
query_fasta = "%s/query.fasta" % output

ref = pysam.FastaFile(fastafile)

flank = 500

tr_records = []
prefix_records = []
suffix_records = []
query_records = []

with open(simple_repeatsbed,'r') as fh:
    header = fh.readline()
    for line in fh:
        line = line.rstrip().split('\t')
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        motif = line[3]
        
        out_name = "%s_%s_%s_%s" % (chrom,start,end,motif)
        
        tr_record = SeqRecord(Seq(motif),id=out_name,name=out_name,description="")
        tr_records.append(tr_record)
        
        prefix_seq = ref.fetch(chrom,start-flank,start)
        prefix_record = SeqRecord(Seq(prefix_seq),id=out_name,name=out_name,description="")
        prefix_records.append(prefix_record)
        
        suffix_seq = ref.fetch(chrom,end,end+flank)
        suffix_record = SeqRecord(Seq(suffix_seq),id=out_name,name=out_name,description="")
        suffix_records.append(suffix_record)
        
        contig_seq = ref.fetch(chrom,start-flank,end+flank)
        query_record = SeqRecord(Seq(contig_seq),id=out_name,name=out_name,description="")
        query_records.append(query_record)
        
SeqIO.write(tr_records, tr_fasta, "fasta")
SeqIO.write(prefix_records, prefix_fasta, "fasta")
SeqIO.write(suffix_records, suffix_fasta, "fasta")
SeqIO.write(query_records, query_fasta, "fasta")

