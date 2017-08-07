#!/usr/bin/env python3

import sys
import re

def get_region(info):
    region_list = re.findall(";Func.refGene=(.+?);", info)
    if region_list[0] == "splicing":
        return region_list[0]
    elif region_list[0] == "exonic":
        region = re.findall(';ExonicFunc.refGene=(.+?);', info)
        return region[0]

def get_af(info, population_in):
    if population_in == "EAS":
        exac_eas_list = re.findall(';ExAC_EAS=(.+?);', info)
        if exac_eas_list[0] != ".":
            exac_eas = exac_eas_list[0]
        else:
            exac_eas =  0
        eas_1000g_list = re.findall(';1000g2015aug_eas=(.+?);', info)
        if eas_1000g_list[0] != ".":
            eas_1000g = eas_1000g_list[0]
        else:
            eas_1000g = 0
        return exac_eas, eas_1000g
    elif population_in == "EUR":
        exac_af_list_1 = re.findall(';ExAC_NFE=(.+?);', info)
        if exac_af_list_1[0] != ".":
            exac_nfe = exac_af_list_1[0]
        else:
            exac_nfe = 0
        exac_af_list_2 = re.findall(';ExAC_FIN=(.+?);', info)
        if exac_af_list_2[0] != ".":
            exac_fin = exac_af_list_1[0]
        else:
            exac_fin = 0
        exac_eur = float(exac_nfe) + float(exac_fin)
        eur_1000g_list = re.findall(';1000g2015aug_eur=(.+?);', info)
        if eur_1000g_list[0] != ".":
            eur_1000g = eur_1000g_list[0]
        else:
            eur_1000g = 0
        return exac_eur, eur_1000g
    elif population_in == "AMR":
        exac_amr_list = re.findall(';ExAC_AMR=(.+?);', info)
        if exac_amr_list[0] != ".":
            exac_amr = exac_amr_list[0]
        else:
            exac_amr = 0
        amr_1000g_list = re.findall(';1000g2015aug_amr=(.+?);', info)
        if amr_1000g_list[0] != ".":
            amr_1000g = amr_1000g_list[0]
        else:
            amr_1000g = 0
        return exac_amr, amr_1000g
    elif population_in == "AFR":
        exac_afr_list = re.findall(';ExAC_AFR=(.+?);', info)
        if exac_afr_list[0] != ".":
            exac_afr = exac_afr_list[0]
        else:
            exac_afr = 0
        afr_1000g_list = re.findall(';1000g2015aug_afr=(.+?);', info)
        if afr_1000g_list[0] != ".":
            afr_1000g = afr_1000g_list[0]
        else:
            afr_1000g = 0
        return exac_afr, afr_1000g
    elif population_in == "SAS":
        exac_sas_list = re.findall(';ExAC_SAS=(.+?);', info)
        if exac_sas_list[0] != ".":
            exac_sas = exac_sas_list[0]
        else: 
            exac_sas = 0
        sas_1000g_list = re.findall(';1000g2015aug_sas=(.+?);', info)
        if sas_1000g_list[0] != ".":
            sas_1000g = sas_1000g_list[0]
        else:
            sas_1000g = 0
        return exac_sas, sas_1000g
    elif population_in == "ALL":
        exac_all_list = re.findall(';ExAC_ALL=(.+?);', info)
        if exac_all_list[0] != ".":
            exac_all = exac_all_list[0]
        else:
            exac_all = 0
        all_1000g_list = re.findall(';1000g2015aug_all=(.+?);', info)
        if all_1000g_list[0] != ".":
            all_1000g = all_1000g_list[0]
        else:
            all_1000g = 0
        return exac_all, all_1000g


def get_prediction(info):
    gerp_list = re.findall(';GERP\+\+_RS=(.+?);', info)
    if gerp_list[0] != ".":
        gerp = gerp_list[0]
    else:
        gerp = 0 
    phylop_list = re.findall("phyloP20way_mammalian=(.+?);", info)
    if phylop_list[0] != ".":
        phylop = phylop_list[0]
    else:
        phylop = 0
    vest_list = re.findall(";VEST3_score=(.+?);", info)
    if vest_list[0] != ".":
        vest = vest_list[0]
    else:
        vest = 0
    cadd_list = re.findall(';CADD_phred=([-+]?\d*\.\d+|\d+);', info)
    if cadd_list != []:
        cadd = cadd_list[0]
    else:
        cadd = 0
    sift_list = re.findall(';SIFT_pred=([A-Z]);', info)
    if sift_list != []:
        sift = sift_list[0]
    else:
        sift = "NA"
    polyphen_list = re.findall(';Polyphen2_HVAR_pred=([A-Z]);', info)
    if polyphen_list != []:
        polyphen = polyphen_list[0]
    else:
        polyphen = "NA"
    return gerp, phylop, vest, cadd, sift, polyphen     

in_file = sys.argv[1]
population = sys.argv[2]

if population not in ["ALL", "EAS", "EUR", "AMR", "AFR", "SAS"]:
    print("\n***invalid population input***\n") 
else:
    for line in open(in_file):
        #skip header
        if line[0] == "#":
            continue
        data = line.strip().split("\t")
        gene_region = get_region(data[7])
        exac_pop, pop_1000g = get_af(data[7], population)
        gerp, phylop, vest, cadd, sift, polyphen = get_prediction(data[7])
        #keep the rare variants
        if float(exac_pop) >= 0.01 or float(pop_1000g) >= 0.01:
            continue
        #keep the lof variants 
        if gene_region == "splicing" or gene_region ==  "frameshift_insertion" or gene_region == "frameshift_deletion" or gene_region == "stopgain" or gene_region == "stoploss":
            print(line.strip())
        #use the prediction algorithms to keep the lgd missenese muations,
        #require the nonsynonymous_SNV occured at the conserved nucliotide and has the cadd>15 with at least one of the other algrithms predicted as deleterious
        elif gene_region == "nonsynonymous_SNV":
            if float(gerp) > 2.0 and float(phylop) > 0.0:
                if float(cadd) > 15.0:
                    if  float(vest) > 0.5 or sift == "D" or polyphen == "P" or polyphen == "D":
                        print(line.strip())


    
        



