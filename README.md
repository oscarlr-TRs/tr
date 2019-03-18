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
### parse_hipstr_vcf.py
```
python python/parse_hipstr_vcf.py <vcf>
```
`parse_hipstr_vcf.py` takes in a VCF file created from HipSTR and outputs a bed file with 5 columns. The 5 columns are chrom, start, end, number of repeat units in hap1 and number of repeat units in hap2.

### extract_tr_loci_from_read_alignment.py
```
python python/extract_tr_loci_from_read_alignment.py <bam> <bed with TRs> <ref> <output dir>
```
`extract_tr_loci_from_read_alignment.py` creates 4 fasta files. These fasta files are used as input for Ummat et al's DP TR detecting algorithm. For every TR coordinate in the `<bed with TRs>`, the sequenced before (prefix) the TR coordinate and after (suffix) the TR coordinate are extracted from the `<ref>`. The TR seq from the `<bed with TRs>` is put into a fasta file. The reads overlapping the TR coordinates are also put into a another fasta file. Every entry in a fasta file has the same entry name in all the other fasta files. The 4 fasta files are:
1. `query.fasta`
2. `prefix.fasta`
3. `suffix.fasta`
4. `tr.fasta`

### extract_tr_loci_from_ref.py
```
python python/extract_tr_loci_from_ref.py <bed> <ref> <dir>
```
This is the same as `extract_tr_loci_from_read_alignment.py` but extracts sequence just from the reference. It is used to detect the number of repeat units in the reference. The `query.fasta` is the reference. The bed file has to have 4 columns: chrom,start,end and the motif seq. The query and prefix is the sequence before and after the input coords, and the `tr.fasta` has the motif seq from the input bed file.
