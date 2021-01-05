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

seq = None

header = ["chrom","start","end","period_size","copies","percent_match","motif","sequence"]    

print "\t".join(header)

with open(trffn, 'r') as fh:
    for line in fh:
        if "Sequence" in line:
            seq = line[10:].rstrip()
        line = line.rstrip().split()
        if len(line) != 15:
            continue
        start = line[0]
        end = line[1]
        period_size = line[2] # Period size of the repeat.
        copies = line[3]
        consensus_period_size = line[4] # Size of consensus pattern (may differ slightly from the period size).
        percent_match = line[5] # Percent of matches between adjacent copies overall.
        alignment_score = line[7] # Alignment score.
        motif = line[13]
        sequence = line[14]
        
        out = [seq,start,end,period_size,copies,percent_match,motif,sequence]
        print "\t".join(out)
