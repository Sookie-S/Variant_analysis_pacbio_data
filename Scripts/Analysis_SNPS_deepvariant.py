#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 14:34:36 2022

@author: sookie
"""

import os

os.chdir("/home/sookie/Documents/BP_project/LAST_METHODOLOGY/quickstart-output/deepvariant_res/")

#%%

cbs8639_gff = {}
cbs8639_gff2 = {}

with open("../../quickstart-testdata/CBS8639.gff","r") as file:
    for line in file:
        if line.startswith("#") == False:
            line=line.split("\n")[0]
            line=line.split("\t")
            if line[2] == "CDS":
#                print(line)
                geneName = line[8].split(";Name=")[0]
                geneName = geneName.split("ID=")[1]
                geneStart = line[3]
                geneStop = line[4]
                print(geneName,line[0],geneStart,geneStop)
                cbs8639_gff[str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)] = geneName
                cbs8639_gff2[geneName] = str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)
cbs8638_gff = {}
cbs8638_gff2 = {}

with open("../../quickstart-testdata/CBS8638.gff","r") as file:
    for line in file:
        if line.startswith("#") == False:
            line=line.split("\n")[0]
            line=line.split("\t")
            if line[2] == "CDS":
                geneName = line[8].split(";Name=")[0]
                geneName = geneName.split("ID=")[1]
                geneStart = line[3]
                geneStop = line[4]
                print(geneName,line[0],geneStart,geneStop)
                cbs8638_gff[str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)] = geneName
                cbs8638_gff2[geneName] = str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)
nrrl_gff = {}
nrrl_gff2 = {}

with open("../../quickstart-testdata/NRRL.gff","r") as file:
    for line in file:
        if line.startswith("#") == False:
            line=line.split("\n")[0]
            line=line.split("\t")
            if line[2] == "CDS":
                geneName = line[8].split(";Name=")[0]
                geneName = geneName.split("ID=")[1]
                geneStart = line[3]
                geneStop = line[4]
                print(geneName,line[0],geneStart,geneStop)
                nrrl_gff[str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)] = geneName
                nrrl_gff2[geneName] = str(line[0])+"-"+str(geneStart)+"-"+str(geneStop)
#%%
nb_genes_3strains = 0
nb_genes_3strains_same_length = 0
dico_genes_same_length = {}
with open("../../quickstart-testdata/CBS8639all_functional_annotation_species_3strains.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
        cbs8639_id=line[0]
        cbs8638_id=line[11]
        nrrl_id=line[12]
        
        if cbs8638_id != "" and cbs8639_id != "" and nrrl_id != "":
            nb_genes_3strains+=1
            print("\n#### ",cbs8639_id,cbs8638_id,nrrl_id)
            for k,v in cbs8639_gff.items():
                if v == cbs8639_id:
                    kbul39 = v
                    elm = k.split("-")
                    chrm = elm[0]
                    start = int(elm[1])
                    stop = int(elm[2])
                    length_cbs8639 = stop-start
                    print("-",v,k,"gene length",stop-start,"at chromosome",chrm)
            for k,v in cbs8638_gff.items():
                if v == cbs8638_id:
                    kbul38 = v
                    elm = k.split("-")
                    chrm = elm[0]
                    start = int(elm[1])
                    stop = int(elm[2])
                    length_cbs8638 = stop-start
                    print("-",v,k,"gene length",stop-start,"at chromosome",chrm)
            for k,v in nrrl_gff.items():
                if v == nrrl_id:
                    kbulNR = v
                    elm = k.split("-")
                    chrm = elm[0]
                    start = int(elm[1])
                    stop = int(elm[2])
                    length_nrrl = stop-start
                    print("-",v,k,"gene length",stop-start,"at chromosome",chrm)
            if length_cbs8638 == length_cbs8639 == length_nrrl:
                nb_genes_3strains_same_length+=1
                dico_genes_same_length[kbul39] = str(kbul38)+"-"+str(kbulNR)
#        break
#%%
print("Total number of genes shared between the 3 strains:",nb_genes_3strains)
print("Number of genes with same length:",nb_genes_3strains_same_length)
print("Number of genes with different length:",nb_genes_3strains-nb_genes_3strains_same_length)        
#%%
x=0
dico_snp_CDS_39 = {}
with open("CBS8639_vs_CBS8639reads_fastq_withoutRefCall_CDS.vcf") as file:
    for line in file:
        if line.startswith("#") == False:
#            print(line)
            line=line.split("\n")[0]
            line=line.split("\t")
            variant_pos = line[1]
            snp_c_alt=str(line[3])+"-"+str(line[4])
            
            if len(line[3]) == 1 and len(line[4]) == 1:
                x+=1
#                print(line[3],line[4],len(line[3]),len(line[4]))
                dico_snp_CDS_39[str(variant_pos)+"-"+line[0]]=snp_c_alt
print("There are",x,"SNP (no complex case)")
#%%

genes_with_snps_39 = []
for snp_info, snp in dico_snp_CDS_39.items():
    for k,v in cbs8639_gff.items():
        snp_pos = snp_info.split("-")[0]
        snp_chrm = snp_info.split("-")[1]
        elm = k.split("-")
        chrm = elm[0]
        start = int(elm[1])
        stop = int(elm[2])
        newPos = int(snp_pos)-int(start)+1
        if int(start) <= int(snp_pos) <= int(stop) and snp_chrm == chrm:
#            print("\n-Gene",v,k)
            if v not in genes_with_snps_39:
                genes_with_snps_39.append(v)
#            print(v,"has a snp",snp,"at the position",snp_pos)
            try:
                globals()["dico_"+str(v)][newPos] = str(snp)
            except:
                globals()["dico_"+str(v)] = {}
                globals()["dico_"+str(v)][newPos] = str(snp)
            break
#%%
x=0
dico_snp_CDS_38 = {}
with open("CBS8638_vs_CBS8638reads_fastq_withoutRefCall_CDS.vcf") as file:
    for line in file:
        if line.startswith("#") == False:
#            print(line)
            line=line.split("\n")[0]
            line=line.split("\t")
            variant_pos = line[1]
            snp_c_alt=str(line[3])+"-"+str(line[4])
            if len(line[3]) == 1 and len(line[4]) == 1:
                x+=1
#                print(line[3],line[4],len(line[3]),len(line[4]))
                dico_snp_CDS_38[str(variant_pos)+"-"+line[0]]=snp_c_alt
print("There are",x,"SNP (no complex case)")

#%%

genes_with_snps_38 = []
for snp_info, snp in dico_snp_CDS_38.items():
    for k,v in cbs8638_gff.items():
        snp_pos = snp_info.split("-")[0]
        snp_chrm = snp_info.split("-")[1]
        elm = k.split("-")
        chrm = elm[0]
        start = int(elm[1])
        stop = int(elm[2])
        newPos = int(snp_pos)-int(start)+1
        if int(start) <= int(snp_pos) <= int(stop) and snp_chrm == chrm:
#            print("\n-Gene",v,k)
            if v not in genes_with_snps_38:
                genes_with_snps_38.append(v)
#            print(v,"has a snp",snp,"at the position",snp_pos)
            try:
                globals()["dico_"+str(v)][newPos] = str(snp)
            except:
                globals()["dico_"+str(v)] = {}
                globals()["dico_"+str(v)][newPos] = str(snp)
            break
#%%
x=0
dico_snp_CDS_NR = {}
with open("NRRL_vs_NRRLreads_fastq_withoutRefCall_CDS.vcf") as file:
    for line in file:
        if line.startswith("#") == False:
#            print(line)
            line=line.split("\n")[0]
            line=line.split("\t")
            variant_pos = line[1]
            snp_c_alt=str(line[3])+"-"+str(line[4])
            if len(line[3]) == 1 and len(line[4]) == 1:
                x+=1
#                print(line[3],line[4],len(line[3]),len(line[4]))
                dico_snp_CDS_NR[str(variant_pos)+"-"+line[0]]=snp_c_alt
print("There are",x,"SNP (no complex case)")
#%%
genes_with_snps_NR = []
for snp_info, snp in dico_snp_CDS_NR.items():
    for k,v in nrrl_gff.items():
        snp_pos = snp_info.split("-")[0]
        snp_chrm = snp_info.split("-")[1]
        elm = k.split("-")
        chrm = elm[0]
        start = int(elm[1])
        stop = int(elm[2])
        newPos = int(snp_pos)-int(start)+1
        if int(start) <= int(snp_pos) <= int(stop) and snp_chrm == chrm:
#            print("\n-Gene",v,k)
            if v not in genes_with_snps_NR:
                genes_with_snps_NR.append(v)
#            print(v,"has a snp",snp,"at the position",snp_pos)
            try:
                globals()["dico_"+str(v)][newPos] = str(snp)
            except:
                globals()["dico_"+str(v)] = {}
                globals()["dico_"+str(v)][newPos] = str(snp)
                break
#%%
for gene in genes_with_snps_39:
    print("\n-- Gene",gene,globals()["dico_"+str(gene)])
    break
#%%
for gene in genes_with_snps_NR:
    print("\n-- Gene",gene,globals()["dico_"+str(gene)])
    break
#%%
for gene in genes_with_snps_38:
    print("\n-- Gene",gene,globals()["dico_"+str(gene)])
    break
#%%
nr_gene_without_snp = 0
cbs38_gene_without_snp = 0
cbs39_gene_without_snp = 0

nr_gene_without_snp_list = []
cbs38_gene_without_snp_list = []
cbs39_gene_without_snp_list = []

## Case 1 38 vs 39 vs NR
SNPpos_sameVariant_38_39_NR = 0
SNPpos_differentVariant_38_39_NR = 0
SNPpos_sameVariant_38_39_NR_bis = 0

## Case 2 39 vs 38
SNPpos_sameVariant_38_39 = 0
SNPpos_differentVariant_38_39 = 0
SNPpos_sameVariant_38_39_bis = 0

## Case 3 39 vs NR
SNPpos_sameVariant_39_NR = 0
SNPpos_differentVariant_39_NR = 0
SNPpos_sameVariant_39_NR_bis = 0

## Case 3 TEST NR vs 39
SNPpos_sameVariant_NR_39 = 0
SNPpos_differentVariant_NR_39 = 0
SNPpos_sameVariant_NR_39_bis = 0

## Case 4 38 vs NR
SNPpos_sameVariant_38_NR = 0
SNPpos_differentVariant_38_NR = 0
SNPpos_sameVariant_38_NR_bis = 0

## Case 4 TEST NR vs 38
SNPpos_sameVariant_NR_38 = 0
SNPpos_differentVariant_NR_38 = 0
SNPpos_sameVariant_NR_38_bis = 0


for k,v in dico_genes_same_length.items():
    cbs8639gene = k
    cbs8638gene = v.split("-")[0]
    NRRLY27gene = v.split("-")[1]
    
    try:
        dico_cbs8639gene = globals()['dico_'+str(cbs8639gene)]
    except:
        cbs39_gene_without_snp+=1
        cbs39_gene_without_snp_list.append(cbs8639gene)
    try:
        dico_cbs8638gene = globals()['dico_'+str(cbs8638gene)]
    except:
        cbs38_gene_without_snp+=1
        cbs38_gene_without_snp_list.append(cbs8638gene)
    try:
        dico_NRRLY27gene = globals()['dico_'+str(NRRLY27gene)]
    except:
        nr_gene_without_snp+=1
        nr_gene_without_snp_list.append(NRRLY27gene)
    
    ####  CASE 1: 38 VS 39 VS NR
    try:
        for k39,v39 in dico_cbs8639gene.items():
            if (k39 in dico_cbs8638gene.keys()) and (k39 in dico_NRRLY27gene.keys()):
                k38 = k39
                v38 = dico_cbs8638gene[k38]
                kNR = k39
                vNR = dico_NRRLY27gene[kNR]
#                print("\n---")
#                print(cbs8639gene,k39,v39,"and",cbs8638gene,k38,v38,"and",NRRLY27gene,kNR,vNR)
                if v39 == v38 == vNR:
#                    print("same variants")
                    SNPpos_sameVariant_38_39_NR+=1
                else:
                    tmp = v38.split("-")
                    v38_2 = str(tmp[1])+"-"+str(tmp[0])
                    tmp2 = vNR.split("-")
                    vNR_2 = str(tmp2[1])+"-"+str(tmp2[0])
                    if (v39 == v38_2 == vNR) or (v39 == v38_2 == vNR_2) or (v39 == v38 == vNR_2):
                        SNPpos_sameVariant_38_39_NR_bis+=1
#                        print("same variants but way around!!")
                    else:
#                        print("different variants")
                        SNPpos_differentVariant_38_39_NR+=1

    except:
        print("ERROR CASE 1")

    ####  CASE 2: 39 VS 38
    try:
        for k39,v39 in dico_cbs8639gene.items():
            if (k39 in dico_cbs8638gene.keys()):
                k38 = k39
                v38 = dico_cbs8638gene[k38]
                
#                print("\n---")
#                print(cbs8639gene,k39,v39,"and",cbs8638gene,k38,v38)
                if v39 == v38:
#                    print("same variants")
                    SNPpos_sameVariant_38_39+=1
                else:
                    tmp = v38.split("-")
                    v38_2 = str(tmp[1])+"-"+str(tmp[0])
                    if (v39 == v38_2):
                        SNPpos_sameVariant_38_39_bis+=1
#                        print("same variants but way around!!")
                    else:
#                        print("different variants")
                        SNPpos_differentVariant_38_39+=1

    except:
        print("")
        
    ####  CASE 3: 39 VS NR
    try:
        for k39,v39 in dico_cbs8639gene.items():
            if (k39 in dico_NRRLY27gene.keys()):
                kNR = k39
                vNR = dico_NRRLY27gene[kNR]
#                print("\n---")
#                print(cbs8639gene,k39,v39,"and",NRRLY27gene,kNR,vNR)
                if v39 == vNR:
#                    print("same variants")
                    SNPpos_sameVariant_39_NR+=1
                else:
                    tmp2 = vNR.split("-")
                    vNR_2 = str(tmp2[1])+"-"+str(tmp2[0])
                    if (v39 == vNR_2):
                        SNPpos_sameVariant_39_NR_bis+=1
#                        print("same variants but way around!!")
                    else:
#                        print("different variants")
                        SNPpos_differentVariant_39_NR+=1

    except:
        print("issue 3")

    ####  CASE 3 TEST: NR VS 39
    try:
#        print(dico_NRRLY27gene)
        for kNR,vNR in dico_NRRLY27gene.items():
#            print(kNR,vNR)
            if kNR in dico_cbs8639gene.keys():
                print("OK")
                k39 = kNR
                v39 = dico_cbs8639gene[k39]
                print("\n---")
                print(cbs8639gene,k39,v39,"and",NRRLY27gene,kNR,vNR)
                if v39 == vNR:
                    print("same variants")
                    SNPpos_sameVariant_NR_39+=1
                else:
                    tmp = v39.split("-")
                    v39_2 = str(tmp[1])+"-"+str(tmp[0])
                    tmp2 = vNR.split("-")
                    vNR_2 = str(tmp2[1])+"-"+str(tmp2[0])
                    if (vNR == v39_2):
                        SNPpos_sameVariant_NR_39_bis+=1
                        print("same variants but way around!!")
                    else:
                        print("different variants")
                        SNPpos_differentVariant_NR_39+=1

    except:
        print("issue 3 test")
        
        
    ####  CASE 4: 38 VS NR
    try:
        for k38,v38 in dico_cbs8638gene.items():
            if (k38 in dico_NRRLY27gene.keys()):
                v38 = dico_cbs8638gene[k38]
                kNR = k38
                vNR = dico_NRRLY27gene[kNR]
#                print("\n---")
#                print(cbs8638gene,k38,v38,"and",NRRLY27gene,kNR,vNR)
                if v38 == vNR:
#                    print("same variants")
                    SNPpos_sameVariant_38_NR+=1
                else:
                    tmp = v38.split("-")
                    v38_2 = str(tmp[1])+"-"+str(tmp[0])
                    tmp2 = vNR.split("-")
                    vNR_2 = str(tmp2[1])+"-"+str(tmp2[0])
                    if (v38_2 == vNR) or (v38 == vNR_2):
                        SNPpos_sameVariant_38_NR_bis+=1
#                        print("same variants but way around!!")
                    else:
#                        print("different variants")
                        SNPpos_differentVariant_38_NR+=1

    except:
        print("issue 4")
    
    ####  CASE 4 test: NR VS 38
    try:
#        print(dico_NRRLY27gene)
        for kNR,vNR in dico_NRRLY27gene.items():
#            print(kNR,vNR)
            if (kNR in dico_cbs8638gene.keys()):
                
                k38 = kNR
                v38 = dico_cbs8638gene[k38]
#                print("\n---")
#                print(cbs8638gene,k38,v38,"and",NRRLY27gene,kNR,vNR)
                if v38 == vNR:
#                    print("same variants")
                    SNPpos_sameVariant_NR_38+=1
                else:
                    tmp = v38.split("-")
                    v38_2 = str(tmp[1])+"-"+str(tmp[0])
                    if (v38_2 == vNR):
                        SNPpos_sameVariant_NR_38_bis+=1
                       
#                        print("same variants but way around!!")
                    else:
#                        print("different variants")
                        SNPpos_differentVariant_NR_38+=1
    except:
        print("issue 4 test")

#%%
print("\n####  CASE 1: 38 VS 39 VS NR")
print("Number of SNPs that are present in the 3 strains at the same position in the gene:", SNPpos_sameVariant_38_39_NR+SNPpos_differentVariant_38_39_NR+SNPpos_sameVariant_38_39_NR_bis)
print("Number of SNPs that are the same in the 3 strains:", SNPpos_sameVariant_38_39_NR)
print("Number of SNPs that are the same in the 3 strains but consensus/alternative inverted:", SNPpos_sameVariant_38_39_NR_bis)
print("Number of SNPs that are different in the 3 strains:", SNPpos_differentVariant_38_39_NR)
#%%
print("\n####  CASE 2: 39 VS 38")
print("Number of SNPs that are present in the 2 strains at the same position in the gene:", SNPpos_sameVariant_38_39+SNPpos_differentVariant_38_39+SNPpos_sameVariant_38_39_bis)
print("Number of SNPs that are the same in the 2 strains:", SNPpos_sameVariant_38_39)
print("Number of SNPs that are the same in the 2 strains but consensus/alternative inverted:", SNPpos_sameVariant_38_39_bis)
print("Number of SNPs that are different in the 2 strains:", SNPpos_differentVariant_38_39)

#%%
print("\n####  CASE 3: 39 VS NR")
print("Number of SNPs that are present in the 2 strains at the same position in the gene:", SNPpos_sameVariant_39_NR+SNPpos_differentVariant_39_NR+SNPpos_sameVariant_39_NR_bis)
print("Number of SNPs that are the same in the 2 strains:", SNPpos_sameVariant_39_NR)
print("Number of SNPs that are the same in the 2 strains but consensus/alternative inverted:", SNPpos_sameVariant_39_NR_bis)
print("Number of SNPs that are different in the 2 strains:", SNPpos_differentVariant_39_NR)

#%%
print("\n####  CASE 3 TEST: NR vs 39")
print("Number of SNPs that are present in the 2 strains at the same position in the gene:", SNPpos_sameVariant_NR_39+SNPpos_sameVariant_NR_39_bis+SNPpos_differentVariant_NR_39)
print("Number of SNPs that are the same in the 2 strains:", SNPpos_sameVariant_NR_39)
print("Number of SNPs that are the same in the 2 strains but consensus/alternative inverted:", SNPpos_sameVariant_NR_39_bis)
print("Number of SNPs that are different in the 2 strains:", SNPpos_differentVariant_NR_39)

#%%
print("\n####  CASE 4: 38 VS NR")
print("Number of SNPs that are present in the 2 strains at the same position in the gene:", SNPpos_sameVariant_38_NR+SNPpos_differentVariant_38_NR+SNPpos_sameVariant_38_NR_bis)
print("Number of SNPs that are the same in the 2 strains:", SNPpos_sameVariant_38_NR)
print("Number of SNPs that are the same in the 2 strains but consensus/alternative inverted:", SNPpos_sameVariant_38_NR_bis)
print("Number of SNPs that are different in the 2 strains:", SNPpos_differentVariant_38_NR)
#%%
print("\n####  CASE 4 TEST: NR VS 38")
print("Number of SNPs that are present in the 2 strains at the same position in the gene:", SNPpos_sameVariant_NR_38+SNPpos_differentVariant_NR_38+SNPpos_sameVariant_NR_38_bis)
print("Number of SNPs that are the same in the 2 strains:", SNPpos_sameVariant_NR_38)
print("Number of SNPs that are the same in the 2 strains but consensus/alternative inverted:", SNPpos_sameVariant_NR_38_bis)
print("Number of SNPs that are different in the 2 strains:", SNPpos_differentVariant_NR_38)

#%%

print("Number of genes with same length present in the 3 strains:", len(dico_genes_same_length.keys()))
print("Total CBS 8639 genes without SNP:",cbs39_gene_without_snp)
print("Total CBS 8638 genes without SNP:",cbs38_gene_without_snp)
print("Total NRRL Y27205 genes without SNP:",nr_gene_without_snp)

#%%


######### HOMOZYGOTE GENES ANALYSIS



#%%
homozygote_genes_in_three_strains=0
homozygote_genes_in_38_39_only=0
homozygote_genes_in_38_NR_only=0
homozygote_genes_in_39_NR_only=0
homozygote_genes_in_38_only=0
homozygote_genes_in_39_only=0
homozygote_genes_in_NR_only=0

homozygote_genes_in_three_strains_list=[]
homozygote_genes_in_38_39_only_list=[]
homozygote_genes_in_38_NR_only_list=[]
homozygote_genes_in_39_NR_only_list=[]
homozygote_genes_in_38_only_list=[]
homozygote_genes_in_39_only_list=[]
homozygote_genes_in_NR_only_list=[]

#%%
for k,v in dico_genes_same_length.items():
    gene39 = k
    gene38 = v.split("-")[0]
    geneNR = v.split("-")[1]
#    print(gene39,gene38,geneNR)
    genes = str(k)+"-"+str(v)
    
    ## case 1 homozygote gene in the three strains
    if (gene39 in cbs39_gene_without_snp_list) and (gene38 in cbs38_gene_without_snp_list) and (geneNR in nr_gene_without_snp_list):
        homozygote_genes_in_three_strains+=1
        homozygote_genes_in_three_strains_list.append(genes)

    ## case 2 homozygote_genes_in_38_39_only
    if (gene39 in cbs39_gene_without_snp_list) and (gene38 in cbs38_gene_without_snp_list) and (geneNR not in nr_gene_without_snp_list):
        homozygote_genes_in_38_39_only+=1
        homozygote_genes_in_38_39_only_list.append(genes)
        
    ## case 3 homozygote_genes_in_38_NR_only
    if (gene39 not in cbs39_gene_without_snp_list) and (gene38 in cbs38_gene_without_snp_list) and (geneNR in nr_gene_without_snp_list):
        homozygote_genes_in_38_NR_only+=1
        homozygote_genes_in_38_NR_only_list.append(genes)
        
    ## case 4 homozygote_genes_in_39_NR_only
    if (gene39 in cbs39_gene_without_snp_list) and (gene38 not in cbs38_gene_without_snp_list) and (geneNR in nr_gene_without_snp_list):
        homozygote_genes_in_39_NR_only+=1
        homozygote_genes_in_39_NR_only_list.append(genes)
        
    ## case 5 homozygote_genes_in_38_only
    if (gene39 not in cbs39_gene_without_snp_list) and (gene38 in cbs38_gene_without_snp_list) and (geneNR not in nr_gene_without_snp_list):
        homozygote_genes_in_38_only+=1
        homozygote_genes_in_38_only_list.append(genes)
        
    ## case 6 homozygote_genes_in_39_only
    if (gene39 in cbs39_gene_without_snp_list) and (gene38 not in cbs38_gene_without_snp_list) and (geneNR not in nr_gene_without_snp_list):
        homozygote_genes_in_39_only+=1
        homozygote_genes_in_39_only_list.append(genes)
        
    ## case 7 homozygote_genes_in_NR_only
    if (gene39 not in cbs39_gene_without_snp_list) and (gene38 not in cbs38_gene_without_snp_list) and (geneNR in nr_gene_without_snp_list):
        homozygote_genes_in_NR_only+=1
        homozygote_genes_in_NR_only_list.append(genes)
#%%

print("Total number of cases where a gene is homozygote in CBS 8638:",len(cbs38_gene_without_snp_list))
print("Total number of cases where a gene is homozygote in CBS 8639:",len(cbs39_gene_without_snp_list))
print("Total number of cases where a gene is homozygote in CBS 8639:",len(nr_gene_without_snp_list))


print("Total number of cases where a gene is homozygote in the three strains:",homozygote_genes_in_three_strains)

print("Total number of cases where a gene is homozygote in CBS8638 and CBS8639 but heterozygote in NRRL:",homozygote_genes_in_38_39_only)

print("Total number of cases where a gene is homozygote in CBS8638 and NRRL but heterozygote in CBS8639:",homozygote_genes_in_38_NR_only)

print("Total number of cases where a gene is homozygote in CBS8639 and NRRL but heterozygote in CBS8638:",homozygote_genes_in_39_NR_only)

print("Total number of cases where a gene is homozygote in CBS8638 but heterozygote in CBS8639 and NRRL:",homozygote_genes_in_38_only)

print("Total number of cases where a gene is homozygote in CBS8639 but heterozygote in CBS8638 and NRRL:",homozygote_genes_in_39_only)

print("Total number of cases where a gene is homozygote in NRRL but heterozygote in CBS8638 and CBS8639:",homozygote_genes_in_NR_only)

#%%

alias_dir = ""

CBS8638_alias_dico ={}
with open(str(alias_dir)+"CBS8638all_functional_annotation_species.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\t")
        kbul = line[0]
        scerID = line[1]
        alias = line[2]
        print(kbul,scerID, alias)
        CBS8638_alias_dico[kbul]=str(scerID)+";"+str(alias)


CBS8639_alias_dico ={}
with open(str(alias_dir)+"CBS8639all_functional_annotation_species.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\t")
        kbul = line[0]
        scerID = line[1]
        alias = line[2]
        print(kbul,scerID, alias)
        CBS8639_alias_dico[kbul]=str(scerID)+";"+str(alias)
        

NRRL_alias_dico ={}
with open(str(alias_dir)+"NRRLY27205all_functional_annotation_species.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\t")
        kbul = line[0]
        scerID = line[1]
        alias = line[2]
        print(kbul,scerID, alias)
        NRRL_alias_dico[kbul]=str(scerID)+";"+str(alias)

#%%
kbulspe_dir = ""


list_cbs8638_strain_spe = []
with open(str(kbulspe_dir)+"CBS8638intersect_all_functional_annotation_species_for_upsetPlot.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line[10])
        if str(line[10])=="1":
            list_cbs8638_strain_spe.append(line[0])

list_cbs8639_strain_spe = []
with open(str(kbulspe_dir)+"CBS8639intersect_all_functional_annotation_species_for_upsetPlot.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line[10])
        if str(line[10])=="1":
            list_cbs8639_strain_spe.append(line[0])

list_nrrl_strain_spe = []
with open(str(kbulspe_dir)+"NRRLY27205intersect_all_functional_annotation_species_for_upsetPlot.csv","r") as file:
    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line[10])
        if str(line[10])=="1":
            list_nrrl_strain_spe.append(line[0])
            

list_genus_spe = []
list_pHtolerantSpecies_spe = []

with open(str(kbulspe_dir)+"CBS8639intersect_all_functional_annotation_species_for_upsetPlot_3strains.csv","r") as file:
#    next(file)
    for line in file:
        line=line.split("\n")[0]
        line=line.split("\t")
#        print(line[13],line[14])
        if str(line[13])=="1":
            list_genus_spe.append(line[0])
        if str(line[14])=="1":
            list_pHtolerantSpecies_spe.append(line[0])

#%%

print("list_cbs8639_strain_spe",len(list_cbs8639_strain_spe))
print("list_cbs8638_strain_spe",len(list_cbs8638_strain_spe))
print("list_nrrl_strain_spe",len(list_nrrl_strain_spe))
print("list_genus_spe",len(list_genus_spe))
print("list_pHtolerantSpecies_spe",len(list_pHtolerantSpecies_spe))

#%%
x_38 = 0
x_38_pos = 0
path_dir = ""
with open(str(path_dir)+"homozygote_genes_CBS8638.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\n")
    for gene in cbs38_gene_without_snp_list:
        x_38+=1
        infos = cbs8638_gff2[gene]
        infos = infos.split("-")
        chrm = infos[0]
        start = infos[1]
        stop = infos[2]
        aliases = CBS8638_alias_dico[gene]
        aliases=aliases.split(";")
        if aliases[0] != "":
            x_38_pos+=1
        out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

#%%
x_39 = 0
x_39_pos = 0
with open(str(path_dir)+"homozygote_genes_CBS8639.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\n")
    for gene in cbs39_gene_without_snp_list:
        x_39+=1
        infos = cbs8639_gff2[gene]
        infos = infos.split("-")
        chrm = infos[0]
        start = infos[1]
        stop = infos[2]
        aliases = CBS8639_alias_dico[gene]
        aliases=aliases.split(";")
        if aliases[0] != "":
            x_39_pos+=1
        out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

#%%
x_nr = 0
x_nr_pos = 0
with open(str(path_dir)+"homozygote_genes_NRRLY27205.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\n")
    for gene in nr_gene_without_snp_list:
        x_nr+=1
        infos = nrrl_gff2[gene]
        infos = infos.split("-")
        chrm = infos[0]
        start = infos[1]
        stop = infos[2]
        aliases = NRRL_alias_dico[gene]
        aliases=aliases.split(";")
        if aliases[0] != "":
            x_nr_pos+=1
        out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

#%%
print("Out of",x_38,"genes of CBS8638,",x_38_pos," have an S. cerevisiae ortholog")
print("Out of",x_39,"genes of CBS8639,",x_39_pos," have an S. cerevisiae ortholog")
print("Out of",x_nr,"genes of NRRL,",x_nr_pos," have an S. cerevisiae ortholog")
#%%
print("CBS8639 has",len(cbs8639_gff2.keys()),"CDS of which",str(len(cbs8639_gff2.keys())-len(genes_with_snps_39)),"are homozygote and",str(len(genes_with_snps_39)),"are heterozygote")
print("CBS8638 has",len(cbs8638_gff2.keys()),"CDS of which",str(len(cbs8638_gff2.keys())-len(genes_with_snps_38)),"are homozygote and",str(len(genes_with_snps_38)),"are heterozygote")
print("NRRL has",len(nrrl_gff2.keys()),"CDS of which",str(len(nrrl_gff2.keys())-len(genes_with_snps_NR)),"are homozygote and",str(len(genes_with_snps_NR)),"are heterozygote")

#%%
#### here isolate all homozygote genes including the ones not shared by three strains

x_38all = 0
x_38all_pos = 0
path_dir = ""
with open(str(path_dir)+"homozygote_genes_CBS8638_all.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\n")
    for gene in cbs8638_gff2.keys():
        if gene not in genes_with_snps_38:
            x_38all+=1
            infos = cbs8638_gff2[gene]
            infos = infos.split("-")
            chrm = infos[0]
            start = infos[1]
            stop = infos[2]
            aliases = CBS8638_alias_dico[gene]
            aliases=aliases.split(";")
            if aliases[0] != "":
                x_38all_pos+=1
            
            if gene in list_cbs8638_strain_spe:
                strain_spe="yes"
            else:
                strain_spe="no"
            out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\n")

x_39all = 0
x_39all_pos = 0
path_dir = ""
with open(str(path_dir)+"homozygote_genes_CBS8639_all.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    for gene in cbs8639_gff2.keys():
        if gene not in genes_with_snps_39:
            x_39all+=1
            infos = cbs8639_gff2[gene]
            infos = infos.split("-")
            chrm = infos[0]
            start = infos[1]
            stop = infos[2]
            aliases = CBS8639_alias_dico[gene]
            aliases=aliases.split(";")
            if aliases[0] != "":
                x_39all_pos+=1
            if gene in list_cbs8639_strain_spe:
                strain_spe="yes"
            else:
                strain_spe="no"
            if gene in list_genus_spe:
                genus_spe="yes"
            else:
                genus_spe="no"
            if gene in list_pHtolerantSpecies_spe:
                phTolspe="yes"
            else:
                phTolspe="no"
            out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")



x_nrall = 0
x_nrall_pos = 0
path_dir = ""
with open(str(path_dir)+"homozygote_genes_NRRLY27205_all.csv","w") as out:
    out.write("Gene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\n")
    for gene in nrrl_gff2.keys():
        if gene not in genes_with_snps_NR:
            x_nrall+=1
            infos = nrrl_gff2[gene]
            infos = infos.split("-")
            chrm = infos[0]
            start = infos[1]
            stop = infos[2]
            aliases = NRRL_alias_dico[gene]
            aliases=aliases.split(";")
            if aliases[0] != "":
                x_nrall_pos+=1
            if gene in list_nrrl_strain_spe:
                strain_spe="yes"
            else:
                strain_spe="no"
            out.write(str(gene)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str()+"\n")


#%%
print("Out of",x_38all,"genes of CBS8638,",x_38all_pos," have an S. cerevisiae ortholog")
print("Out of",x_39all,"genes of CBS8639,",x_39all_pos," have an S. cerevisiae ortholog")
print("Out of",x_nrall,"genes of NRRL,",x_nrall_pos," have an S. cerevisiae ortholog")

#%%
print(homozygote_genes_in_three_strains_list)

#%%

#%%
################## CREATE SPREADSHEETS
        
#%%
print("Total number of cases where a gene is homozygote in the three strains:",homozygote_genes_in_three_strains)

y = 0
x_38_39_nr = 0
x_38_39_nr_pos = 0
with open(str(path_dir)+"homozygote_genes_in_three_strains.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8639 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8638 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_three_strains_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39_nr+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")

#%%
print("Total number of cases where a gene is homozygote in CBS8638 and CBS8639 but heterozygote in NRRL:",homozygote_genes_in_38_39_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_CBS8639_CBS8638_but_heterozygote_in_NRRLY27205.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8639 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8638 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_38_39_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")



#%%
print("Total number of cases where a gene is homozygote in CBS8638 and NRRL but heterozygote in CBS8639:",homozygote_genes_in_38_NR_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_CBS8638_NRRL_but_heterozygote_in_CBS8639.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8638 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8639 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_38_NR_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")



#%%
print("Total number of cases where a gene is homozygote in CBS8639 and NRRL but heterozygote in CBS8638:",homozygote_genes_in_39_NR_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_CBS8639_NRRLY27205_but_heterozygote_in_CBS8638.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8639 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8638 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_39_NR_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")



#%%
print("Total number of cases where a gene is homozygote in CBS8638 but heterozygote in CBS8639 and NRRL:",homozygote_genes_in_38_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_CBS8638_but_heterozygote_in_CBS8639_NRRLY27205.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8638 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8639 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_38_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")



#%%
print("Total number of cases where a gene is homozygote in CBS8639 but heterozygote in CBS8638 and NRRL:",homozygote_genes_in_39_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_CBS8639_but_heterozygote_in_CBS8638_NRRLY27205.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi CBS8639 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8638 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi NRRLY-27205 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_39_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")



#%%
print("Total number of cases where a gene is homozygote in NRRL but heterozygote in CBS8638 and CBS8639:",homozygote_genes_in_NR_only)

y = 0
x_38_39 = 0
with open(str(path_dir)+"genes_homozygote_in_NRRLY27205_but_heterozygote_in_CBS8639_CBS8638.csv","w") as out:
#    out.write("Group number\tGene name\tChromosome\tStart position\tStop position\tS. cerevisiae ortholog\tGene alias\tStrain specific?\tGenus specific?\tpH tolerant specific?\n")
    out.write("K. bulderi NRRLY-27205 gene"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8639 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"K. bulderi CBS8638 gene (heterozygote!)"+"\t"+"Chromosome"+"\t"+"start-stop position"+"\t"+"S. cerevisiae ortholog ID"+"\t"+"Gene alias"+"\n")

    for geneslist in homozygote_genes_in_NR_only_list:
        y+=1
        geneslist=geneslist.split("-")
        gene39 = geneslist[0]
        gene38 = geneslist[1]
        geneNR = geneslist[2]
        x_38_39+=1

        infos = nrrl_gff2[geneNR]
        infos = infos.split("-")
        chrmNR = infos[0]
        startNR = infos[1]
        stopNR = infos[2]
        aliasesNR = NRRL_alias_dico[geneNR]
        aliasesNR=aliasesNR.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        
#        out.write(str(y)+"\t"+str(geneNR)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8638_gff2[gene38]
        infos = infos.split("-")
        chrm38 = infos[0]
        start38 = infos[1]
        stop38 = infos[2]
        aliases38 = CBS8638_alias_dico[gene38]
        aliases38=aliases38.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
#        out.write(str(y)+"\t"+str(gene38)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\n")

        infos = cbs8639_gff2[gene39]
        infos = infos.split("-")
        chrm39 = infos[0]
        start39 = infos[1]
        stop39 = infos[2]
        aliases39 = CBS8639_alias_dico[gene39]
        aliases39=aliases39.split(";")
#        if aliases[0] != "":
#            x_nr_pos+=1
        if gene in list_cbs8639_strain_spe:
            strain_spe="yes"
        else:
            strain_spe="no"
        if gene in list_genus_spe:
            genus_spe="yes"
        else:
            genus_spe="no"
        if gene in list_pHtolerantSpecies_spe:
            phTolspe="yes"
        else:
            phTolspe="no"
#        out.write(str(y)+"\t"+str(gene39)+"\t"+str(chrm)+"\t"+str(start)+"\t"+str(stop)+"\t"+str(aliases[0])+"\t"+str(aliases[1])+"\t"+str(strain_spe)+"\t"+str(genus_spe)+"\t"+str(phTolspe)+"\t"+str()+"\n")

        out.write(str(geneNR)+"\t"+str(chrmNR)+"\t"+str(startNR)+"-"+str(stopNR)+"\t"+str(gene39)+"\t"+str(chrm39)+"\t"+str(start39)+"-"+str(stop39)+"\t"+str(gene38)+"\t"+str(chrm38)+"\t"+str(start38)+"-"+str(stop38)+"\t"+str(aliases39[0])+"\t"+str(aliases39[1])+"\n")


#%%