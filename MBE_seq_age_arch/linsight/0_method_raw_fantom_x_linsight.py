#!/usr/bin/env python
# coding: utf-8

# intersect fantom enhancers w/ LINSIGHT
# In[1]:


import glob
import os, sys

import datetime
last_run = datetime.datetime.now()
today = (datetime.date.today())
print("last run", datetime.datetime.now())

AP = "/home/fongsl/scp/%s/"% today
AP_cmd = "mkdir %s" %AP
os.system(AP_cmd)


# In[2]:


# linsight intersection
linsight_dc_path = "/dors/capra_lab/data/evolutionary_conservation/linsight/"

linsight_f = "/dors/capra_lab/data/evolutionary_conservation/linsight/LINSIGHT.bed"

linsight_path = "/dors/capra_lab/projects/enhancer_ages/linsight/data/"


fantom_path = "/dors/capra_lab/projects/enhancer_ages/fantom/data/"

fantom_fs = ["%sall_fantom_erna_112_tissue.bed" % fantom_path]

fantom_fs = ["%sall_unique_fantom_erna_112_tissue.bed" % fantom_path]


# In[3]:


#split up the linsight_files

os.chdir(linsight_dc_path)

split_linsight = "awk '{print >$1\"_linsight.bed\"}' %s" %linsight_f

print(split_linsight)

#os.system(split_linsight)


# In[4]:


# make a dictionary of chr*_linsight.bed

linsight_chrs = glob.glob("%schr*_linsight.bed" % (linsight_dc_path)) # get all the linsight chromosomes

linsight_dict = {} # {chrX: chrX_linsight.bed}

for lin_chr in linsight_chrs:

    chr_num = (lin_chr.split("/")[-1]).split("_")[0]

    linsight_dict[chr_num] = lin_chr


# In[5]:


# score fantom enhancers with linsight scores

for fantom_f in fantom_fs:

    arch_id = (fantom_f.split("/")[-1]).split(".")[0] # get the architecture

    os.chdir(fantom_path)

    split_fantom = "awk '{print >$1\"_%s.bed\"}' %s" % (arch_id, fantom_f)

    print(split_fantom)

    os.system(split_fantom)

    fantom_chrs = glob.glob("%schr*_%s.bed" % (fantom_path, arch_id)) # get all the linsight chromosomes

    for fantom_chr in fantom_chrs:

        chr_num = (fantom_chr.split("/")[-1]).split("_")[0]
        if chr_num != "chrX":

            linsight_chr = linsight_dict[chr_num]

            bed_out = "%s%s_%s_linsight.bed" %(linsight_path, chr_num, arch_id)

            bed_cmd = "bedtools intersect -a %s -b %s -wao > %s" %(fantom_chr, linsight_chr, bed_out) # 10% of enhancer must overlap linsight score
            print(bed_cmd)

            os.system(bed_cmd)

    # cat chromosomes
    cat_out = "%s%s_linsight.bed" %(linsight_path, arch_id)
    cat_cmd = "cat %schr*_%s_linsight.bed > %s" %(linsight_path, arch_id, cat_out)
    os.system(cat_cmd)

    # clean up
    cleanup_cmd = "rm %schr*_%s_linsight.bed" %(linsight_path, arch_id)
    os.system(cleanup_cmd)


# In[7]:


fantom_chr


# In[8]:


cat_out = "%s%s_linsight.bed" %(linsight_path, arch_id)
cat_cmd = "cat %schr*_%s_linsight.bed > %s" %(linsight_path, arch_id, cat_out)
os.system(cat_cmd)

# clean up
cleanup_cmd = "rm %schr*_%s_linsight.bed" %(linsight_path, arch_id)
os.system(cleanup_cmd)


# In[ ]:
