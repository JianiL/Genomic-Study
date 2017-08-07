#!/usr/bin/env python3
import sys
import re
#the allele frequency used in the map file is from the ExAC
#different population are all considered here
#currently only limited to one population, if the cases is mixed,
#I will use the overall freqency in the population
def get_af(info, population_in):
    if population_in == "EAS":
        exac_eas_list = re.findall(';ExAC_EAS=(.+?);', info)
        if exac_eas_list[0] != ".":
            exac_eas = exac_eas_list[0]
        else:
            exac_eas =  0
        return exac_eas
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
        return exac_eur
    elif population_in == "AMR":
        exac_amr_list = re.findall(';ExAC_AMR=(.+?);', info)
        if exac_amr_list[0] != ".":
            exac_amr = exac_amr_list[0]
        else:
            exac_amr = 0
        return exac_amr
    elif population_in == "AFR":
        exac_afr_list = re.findall(';ExAC_AFR=(.+?);', info)
        if exac_afr_list[0] != ".":
            exac_afr = exac_afr_list[0]
        else:
            exac_afr = 0 
        return exac_afr
    elif population_in == "SAS":
        exac_sas_list = re.findall(';ExAC_SAS=(.+?);', info)
        if exac_sas_list[0] != ".":
            exac_sas = exac_sas_list[0]
        else: 
            exac_sas = 0 
        return exac_sas
    elif  population_in == "ALL":
        exac_all_list = re.findall(';ExAC_ALL=(.+?);', info)
        if exac_all_list[0] != ".":
            exac_all = exac_all_list[0]
        else:
            exac_all = 0
        return exac_all


#the order of input and out put file
input = sys.argv[1]
population = sys.argv[2]
out_map = sys.argv[3]
out_tped = sys.argv[4]

if population not in ["ALL", "EAS", "EUR", "AMR", "AFR", "SAS"]:
    print("\n***invalid population input***\n") 
else: 

    map_file = open(out_map, 'w')
    tped_file = open(out_tped, 'w')

    for line in open(input):
        genotype = []
        if line[0] == "#" or line[0] == "X" or line[0] == "Y" or  line[0] == "c":
            continue
        data = line.strip().split("\t")
        if data[2] != ".":
            var_id = data[2]
        else:
            var_id = ":".join([data[0], data[1]])
        gene = re.findall(';Gene.refGene=(.+?);', data[7].strip())
        af = get_af(data[7], population)
        for i in range(9, len(data)):
            if ":lowGQ:" not in data[i]:
                gt = data[i].strip().split(":")
                if gt[0] == "./.":
                    genotype.append("-9")
                    genotype.append("-9")
                elif gt[0] == "0/1":
                    genotype.append("0")
                    genotype.append("1")
                elif gt[0] == "1/1":
                    genotype.append("1")
                    genotype.append("1")
                elif gt[0] == "0/0":
                    genotype.append("0")
                    genotype.append("0")
                else:
                    haplo = gt[0].strip().split("|")
                    genotype.append(haplo[0])
                    genotype.append(haplo[1])
            else:
                genotype.append("-9")
                genotype.append("-9")

        if gene[0] != "." or gene[0] != "unknown":       
            map_file.write(" ".join([gene[0], var_id, str(af)]) + "\n")
            tped_file.write(var_id + " " + " ".join(genotype) + "\n")

    map_file.close()
    tped_file.close()

