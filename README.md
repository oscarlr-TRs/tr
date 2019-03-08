# tr - code for tandem repeat projects

[Comparison code](#comparison-code)  

# Comparison code
### compare_trs_between_two_beds.py  
```
python python/compare_trs_between_two_beds.py  <bed1> <bed2>
```
`compare_trs_between_two_beds.py` takes in two bed files and outputs a single bed file. The format of the input bed file is `chrom start end number_of_motifs_in_hap1 number_of_motifs_in_hap2`. It outputs a concatenation of both files with the number of motifs different between both call sets added to each line.
