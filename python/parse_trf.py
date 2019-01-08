#!/bin/env python
import sys

trffn = sys.argv[1]

# @pbsv.DEL.3.REF
# 0 1 
# 1 22 
# 2 2 
# 3 11.0 
# 4 2 
# 5 100 
# 6 0 
# 7 44 
# 8 50 
# 9 0 
# 10 0
# 11 50 
# 12 1.00 
# 13 AT 
# 14 ATATATATATATATATATATAT 
# 15 . 
# 16 .

out = ["chrom","chrom_start","chrom_end","sv_id","start","end","repeat_motif_size","copies","percent_rep","repeat_motif","period_size","percent_match","repeat_seq_length","repeat_seq"]
print "\t".join(out)
with open(trffn, 'r') as trffh:
    for line in trffh:
        if line[0] == "@":
            line = line.rstrip()
            line = line[1:].rstrip().split("_")
            sv_id = line[0]
            chrom = line[1]
            chrom_start = line[2]
            chrom_end = line[3]
        else:
            line = line.rstrip().split()
            start = line[0]
            end = line[1]
            period_size = line[2]
            copies = line[3]
            consensus_period_size = line[4]
            percent_match = line[5]
            repeat_motif = line[13]
            repeat_seq = line[14]
            seq_length = 0.0
            if line[14] != ".":
                seq_length += len(line[14])
            if line[15] != ".":
                seq_length += len(line[15])
            if line[16] != ".":
                seq_length += len(line[16])
            percent_seq = len(repeat_seq)/seq_length
            repeat_motif_size = len(repeat_motif)
            if not chrom_start.isdigit():
                continue
            out = [chrom,chrom_start,chrom_end,sv_id,start,end,repeat_motif_size,copies,percent_seq,repeat_motif,period_size,percent_match,len(line[14]),repeat_seq]
            print "\t".join(map(str,out))
