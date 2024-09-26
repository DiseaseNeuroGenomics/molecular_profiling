####makes files for mmqtl

##### generate mmqtl script per chromosome
for i in {1..22}
do
cat > gaba_qtl${i}.sh << EOT
#!/bin/bash  
#BSUB -J gaba_qtl${i}
#BSUB -q premium
#BSUB -P acc_
#BSUB -n 1
#BSUB -W 06:00
#BSUB -L /bin/bash
#BSUB -o %J.stdout
#BSUB -eo %J.stderr


for x in \$(cat genenames${i}.txt); do ../MMQTL-master/MMQTL26 -b  -P  pheno_file.txt   -Z  geno_file${i}.txt   -R grm_file.txt -a GeneInfo.txt  -A random   -T \${x}; done                                                                  
EOT
done


##### split geneinfo bed by chromosome
for i in {1..22}; do awk -v chrom=chr${i} '{ if($1 == chrom) {print $4}}' GeneInfo.txt > genenames${i}.txt; done 



## generate "genofile" per chromosome for mmQTL 
## its just the path of each plink file per chr
##also
for i in {1..22};  do
cat > geno_file${i}.txt << EOT
merged_rna_filtered${i}
EOT
done
