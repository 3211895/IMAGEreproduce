# Data processing steps
1.Trimmed the adaptors in the fastq files using TrimGalore!
```
./trim_galore -a <adapter> --gzip --rrbs --length 15 --stringency 4 FILE
```
2.Mapped to the baboon reference genome using BSseeker2
```
python bs_seeker2-align.py -i FILE.fq --aligner=bowtie2 -o FILE.bam -g reference.fa -r --bt2-N 0 --bt2-L 25
python bs_seeker2-call_methylation.py -i FILE.bam -o FILE --db reference 
```
3.Called SNPs using CGmaptools with the BayesWC strategy
```
cgmaptools  snv -i FILE.ATCGmap.gz -m bayes -v FILE.vcf -o FILE.snv  –bayes-e=0.01 –bayes-dynamicP -a
```
4.Merged and filtered VCF files using VCFtools
```
vcftools --vcf FILE.vcf --remove-filtered-all --recode --recode-INFO-all --out FILE
bgzip FILE.vcf
tabix -p vcf FILE.vcf.gz
vcf-merge FILE1.vcf.gz FILE2.vcf.gz ... | bgzip -c > FILE.vcf.gz 
vcftools --vcf FILE.vcf –minDP 3 --max-missing 0.5 --maf 0.05 --recode --recode-INFO-all  --out FILE
vcftools --vcf FILE.vcf --indv n --recode --recode-INFO-all --out FILE
vcftools --vcf FILE.vcf --max-missing 1 --recode --recode-INFO-all --out FILE
```
5.Called allele-specific methylation using CGmaptools
```
cgmaptools asm -r reference.fa -b FILE.bam -l FILE.vcf > FILE.asm
```
