# tr - code for tandem repeat projects

[Comparison code](#comparison-code)  
[Parsing code](#parsing-code)  

# Comparison code
### compare_trs_between_two_beds.py  
```
python python/compare_trs_between_two_beds.py  <bed1> <bed2>
```
`compare_trs_between_two_beds.py` takes in two bed files and outputs a single bed file. The format of the input bed file is `chrom start end number_of_motifs_in_hap1 number_of_motifs_in_hap2`. It outputs a concatenation of both files with the number of motifs different between both call sets added to each line.

# Parsing code
```
python python/parse_hipstr_vcf.py <vcf>
```
`parse_hipstr_vcf.py` takes in a VCF file created from HipSTR and outputs a bed file with 5 columns. The 5 columns are chrom, start, end, number of repeat units in hap1 and number of repeat units in hap2.
