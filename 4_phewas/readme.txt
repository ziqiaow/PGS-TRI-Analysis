#readme

#step 1: 
Download data from the phewas website
https://www.omicspred.org/

#step 2: (can skip if we are only using the SNP names)
Overlift the ASD data from hg38 to hg37/hg19, given that the weights and SNPs are in hg19 format

#step 3:
Given that the data in the website has already removed ambiguous snps, let's directly calculate the PRS scores in PLINK. Combine the original imputed genotype data files into one


#step 4:
put the names into a list, remove "_model"
for f in *.txt; do mv "$f" "${f%%_*.txt}.txt"; done
save names to a list
ls -1 | sed -e 's/\.txt$//' > /dcs04/nilanjan/data/zwang/phewas/data/Metabolon.txt

#no need if saved the txt file in another directory
remove first line
#sed '1d' SomaScan.txt > tmpfile; mv tmpfile SomaScan.txt


#add rsid to OC data, generated from bim files in 1000g (oc data was imputed based on 1000g)
awk '{ print $1":"$4, $2 }' chr_all.bim > 1kg_rsid_hg19.txt

awk '{print $2":"$3,$4,$6}' OPGS002837.txt > test.txt