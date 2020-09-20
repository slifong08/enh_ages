#!/usr/bin/env python

"""
SNPassoc2LDpartners.py - David Rinker 2019

 Change log:

 -----------------------------------------------------------------------------

USAGE:  assoc2LDSNPs.py POP SNP_2_assoc.bed

-POP	1KG POPULATION ("EAS", "EUR", "SAS")
-SNP_2_assoc.bed	input SNP-to-association bed file of tab delimited format: [CHR, POS-1, POS, ASSOCgeneORtrait, pvalue, notes]

IMPORTANT: must use python 3.5 or higher

NOTE: the 4th field of the input file can be anything you want: GWAS study ID, rsID of SNP, effect size, etc (ie. it can be a string or a float)

"""

import os
import sys
import pymysql

POP = sys.argv[1]
# POP = "EUR"

SNP_FILE = sys.argv[2]
# SNP_FILE = "test.SNP_2_assoc.bed"

minLDrequirment= 1.0	#must be 1.0 because r2 LD is not transitive

###FUNCTIONS####

def bed_2_ChrPos (in_file):
	out_list=[]
	for line in open(in_file):
		if not line.startswith("#"):
			seg = line.rstrip().split('\t')
			CHR = seg[0].split('r')[1]
			POS = int(seg[2])
			TRAIT = "na"
			PVAL = "na"
			ES = "na"
			out_list += [(CHR,POS,TRAIT,PVAL,ES)]
	return out_list

def msqlquery (MYSQLTABLE, RSQ, POS): 
	ldSNPs = []
	cursor = cnx.cursor()
	
	query = ("SELECT POS1 FROM "+ MYSQLTABLE + " WHERE (CHROM = %s) AND (Rsquared >= %s) AND (POS2 = %s)")
	cursor.execute(query, (CHR,RSQ,POS))
	for row in cursor:
		ldSNPs.append(str(row['POS1']))
		
	query = ("SELECT POS2 FROM "+ MYSQLTABLE + " WHERE (CHROM = %s) AND (Rsquared >= %s) AND (POS1 = %s)")
	cursor.execute(query, (CHR,RSQ,POS))
	for row in cursor:
		ldSNPs.append(str(row['POS2']))
		
	return ldSNPs

#################################
##_Get LD partners for SNPs in Input SNP from pre-computed 1000 Genomes LD database
##_IMPORTANT NOTE: The LD database was constructed using +/- 500kb windows around queried SNP and only SNPs with Rsquared > 0.5 were retained
##for reference, vcftools command line was: vcftools --gzvcf AFR.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz --hap-r2 --min-r2 0.5 --ld-window-bp 500000 --out AFR.chr1.phase3.ld0.5.collection
	
if POP=='EUR':
	DATABASETABLE='SNP_LD_1KG_PHASE3v5'
else:
	DATABASETABLE='SNP_LD_1KG_PHASE3v5_' + POP

snp_input = bed_2_ChrPos(SNP_FILE)

#make list of unique positions
SNP2assoc = {}	#dict keyed on (CHR,POS) of input asociated SNPs
SNP2LDSNPs = {}	#dict keyed on (CHR,POS) with values of LD snps (CHR,POS)

for CHR,POS,ASSOCgeneORtrait,pvalue,effectsize in snp_input:
	if (CHR,POS) not in SNP2assoc:
		SNP2assoc[(CHR,POS)] = [(CHR,POS,ASSOCgeneORtrait,pvalue,effectsize)]
	else:
		SNP2assoc[(CHR,POS)].append((CHR,POS,ASSOCgeneORtrait,pvalue,effectsize))

snp_input_list = list(dict.fromkeys(SNP2assoc)) #remove duplicate positions by creating a dictionary, using the List items as keys.

##check output
# for keys,values in SNP2assoc.items():
	# print(keys)
	# print(values)


cnx = pymysql.connect(host='vgi01.accre.vanderbilt.edu',
							 user='rinkerd',
							 password='wretched-calculator',
							 db='1kg_snpld_phase3',
							 charset='utf8mb4',
							 cursorclass=pymysql.cursors.DictCursor)
for snp in snp_input_list:
	CHR=snp[0]
	POS=snp[1]
	ldsnps = []
	
	ldsnps = msqlquery(DATABASETABLE, minLDrequirment, POS)
	
	# print '[%s]' % ', '.join(map(str, ldsnps))
	
	if not ldsnps:
		continue
	else:
		for ldsnp in ldsnps:
			if (CHR,POS) not in SNP2LDSNPs:
				SNP2LDSNPs[(CHR,POS)] = [(CHR,int(ldsnp))]
			else:
				SNP2LDSNPs[(CHR,POS)].append((CHR,int(ldsnp)))
cnx.close()

# check output
# for keys,values in SNP2LDSNPs.items():
	# print(keys)
	# print(values)

for snp in snp_input_list:
	CHR=snp[0]
	POS=snp[1]
	if (CHR,POS) not in SNP2LDSNPs:
		for assoc in SNP2assoc[(CHR,POS)]:
			print("chr" + CHR, POS-1, POS, assoc[2], assoc[3], assoc[4], "chr"+str(assoc[0])+":"+str(assoc[1]), sep = '\t')	#print input SNP with all associations
		continue
	else:
		for assoc in SNP2assoc[(CHR,POS)]:
			print("chr" + CHR, POS-1, POS, assoc[2], assoc[3], assoc[4], "chr"+str(assoc[0])+":"+str(assoc[1]), sep = '\t')	#print input SNP with all associations
			for ldsnp in SNP2LDSNPs[(CHR,POS)]:
				print("chr" + ldsnp[0],ldsnp[1]-1, ldsnp[1], assoc[2], assoc[3], assoc[4], "chr"+str(assoc[0])+":"+str(assoc[1]), sep = '\t')	#print LD SNPs with all LD-associated associations
