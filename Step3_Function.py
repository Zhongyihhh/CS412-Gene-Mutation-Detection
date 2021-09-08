# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 17:38:50 2019

"""
import numpy as np
import os
from matplotlib import pyplot as plt


def cal_KL(prob_list1, prob_list2):
    temp = 0
    pseudo = 0.0001
    for i in range(len(prob_list1)):
        temp += prob_list2[i] * np.log2((prob_list2[i]+pseudo)/(prob_list1[i]+pseudo))
    return temp

def read_Motif(loc):
    f = open(loc, 'r')
    line = f.readline()
    line = f.readline()
    Motif = {}
    flag = 0
    temp = []
    while(flag < 10):
        if (str(line.split()[0]) != "<"):
            temp += [float(i) for i in line.split()]
            line = f.readline()
        elif (str(line.split()[0]) == "<"):
            Motif[flag] = temp.copy()
            temp = []
            flag += 1
            line = f.readline()
            line = f.readline()
    return Motif

def main_Fun(folder_loc):
    Predicted_Motif = read_Motif(folder_loc + '/predictedmotif.txt')
    Given_Motif = read_Motif(folder_loc + '/motif.txt')
    All = []
    for i in range(10):
        All.append(cal_KL(Predicted_Motif[i], Given_Motif[i]))
    average_KL = sum(All)/10
    return (average_KL, All)

def site_eval(file_org,file_pred): # input: path of origin and prediction file


    f1 = open(file_org,'r')
    f2 = open(file_pred,'r')

    ovlp_pst = [] # total number of overlapping positions for each trial
    ovlp_sites = [] # total number of overlapping sites for each trial

    for i in range (10): # 10 trials

        # Get original and prediction information
        location_org = []
        site_org = []
        info = f1.readline()
        while (info != '\n'):
            info = info.split(', ')

            location_org += [int(info[0])]
            site_org += [info[1].rstrip('\n')]
            info = f1.readline()
            if info == '':
                break
        SC = len(location_org)
        ML = len(site_org[1])


        location_pred = []
        site_pred = []
        for j in range(SC):
            info = f2.readline()

            info = info.split(', ')
            location_pred += [int(info[0])]
            site_pred += [info[1].rstrip('\n')]
        info = f2.readline()

        ovlp_pst += [sum([max(ML-abs(x-y),0)  for x,y in zip(location_org,location_pred)])]

        sites_count = 0
        for j in range(SC):
            ovlp_sites_ct = 0
            for k in range(ML):
                if site_org[j][k] == site_pred[j][k]:
                    ovlp_sites_ct += 1
            if ovlp_sites_ct >= ML/2:
                sites_count += 1
        ovlp_sites += [sites_count]

    return(ovlp_pst,ovlp_sites)

def run_Evaluation():
#    print (os.getcwd())
    up_dirs = [d for d in os.listdir(os.getcwd()) if os.path.isdir(os.path.join(os.getcwd(), d))]
    Default = [2.0, 8, 10, 500]
    Default_str = 'ICPC: ' + str(2.0) + ', ML: ' + str(8) + ', SC: ' + str(10) + ', SL: 500'
    Combinations = []
    Combinations_sort = []
    record_All_KL = {}
    record_avg_KL = {}
    record_std_KL = {}
    record_avg_pos = {}
    record_std_pos = {}
    record_avg_sites = {}
    record_std_sites = {}
    record_overlap_pos = {}
    record_overlap_sites = {}
    for i in up_dirs:
        i_trim = str(i).split("_")
        try:
            temp_ICPC = i_trim[0]
            temp_ML = i_trim[1]
            temp_SC = i_trim[2]
            temp = 'ICPC: ' + str(float(temp_ICPC)) + ', ML: ' + str(temp_ML) + ', SC: ' + str(temp_SC) + ', SL: 500'
            if float(temp_ICPC) != Default[0] or int(temp_ML) != Default[1] or int(temp_SC) != Default[2]:
                Combinations_sort.append([float(temp_ICPC), int(temp_ML), int(temp_SC), 500])
            KL, All_KL = main_Fun(i)
            record_avg_KL[temp] = KL
            record_std_KL[temp] = np.std(All_KL.copy())
            record_All_KL[temp] = All_KL.copy()
            file_pred = i + '/predictedsites.txt' # file name of predicted sites
            file_org = i + '/sites.txt'# file name of original sites
            overlap_pos, overlap_sites = site_eval(file_org,file_pred)
            record_avg_pos[temp] = sum(overlap_pos)/len(overlap_pos)
            record_std_pos[temp] = np.std(overlap_pos)
            record_avg_sites[temp] = sum(overlap_sites)/len(overlap_sites)
            record_std_sites[temp] = np.std(overlap_sites)
            record_overlap_pos[temp] = [i/(int(temp_SC) * int(temp_ML)) for i in overlap_pos]
            record_overlap_sites[temp] = [i/(int(temp_SC)) for i in overlap_sites]
        except:
            pass


    Combinations_sort.sort()
    for j in Combinations_sort:
        Combinations.append('ICPC: ' + str(float(j[0])) + ', ML: ' + str(j[1]) + ', SC: ' + str(j[2]) + ', SL: 500')
    #Combinations.append(Default_str)

    ICPC_Combinations = [Combinations[0], Combinations[1], Default_str]
    ML_Combinations = [Combinations[2], Combinations[3], Default_str]
    SC_Combinations = [Combinations[4], Default_str, Combinations[5]]
    All_Combinations = ICPC_Combinations + ML_Combinations + SC_Combinations
    #print (record_std_KL[All_Combinations[0]])
    All = {}
    All['ICPC Variation'] = ICPC_Combinations
    All['ML Variation'] = ML_Combinations
    All['SC Variation'] = SC_Combinations

    # Plot the results
    i = 1
    for j in list(All.keys()):
        plt.figure(i)
        fid, ax = plt.subplots()
        str_title = 'Results: KL Divergence of ' + j
        ax.set_title(str_title)
        dataList = []
        legendList = []
        count = 1
        for k in All[j]:
            dataList.append(record_All_KL[k])
            legendList.append(str(count) + ":" + k)
            count += 1
        bp = ax.boxplot(dataList)
        ax.legend([bp["boxes"][0]]*3, legendList, loc=9, bbox_to_anchor=(0.5, -0.1))

        k = ''
        plt.figure(i+1)
        fid, ax = plt.subplots()
        ax.set_title('Results: Overlapped Position of ' + j)
        dataList = []
        legendList = []
        count = 1
        for k in All[j]:
            dataList.append(record_overlap_pos[k])
            legendList.append(str(count) + ":" + k)
            count += 1
        bp = ax.boxplot(dataList)
        ax.legend([bp["boxes"][0]]*3, legendList, loc=9, bbox_to_anchor=(0.5, -0.1))

        k = ''
        plt.figure(i+2)
        fid, ax = plt.subplots()
        ax.set_title('Results: Overlapped Sites of ' + j)
        dataList = []
        legendList = []
        count = 1
        for k in All[j]:
            dataList.append(record_overlap_sites[k])
            legendList.append(str(count) + ":" + k)
            count += 1
        bp = ax.boxplot(dataList)
        ax.legend([bp["boxes"][0]]*3, legendList, loc=9, bbox_to_anchor=(0.5, -0.1))

        i = i + 3

    file_wt = open("Results.txt",'wt')
    file_wt.writelines('Combinations;Avg_KL;Std_KL;Avg_Pos;Std_Pos;Avg_Sites;Std_Sites\n')
    i = 0
    j = 0
    for t in All_Combinations:
        if i%3 == 0:
            file_wt.writelines('\n')
            file_wt.writelines(list(All.keys())[j] + '\n')
            j += 1
        file_wt.writelines(t + ";" + str(record_avg_KL[t]) + ";" + str(record_std_KL[t]) + ";" + str(record_avg_pos[t]) + ";" + str(record_std_pos[t]) + ";" + str(record_avg_sites[t]) + ";" + str(record_std_sites[t]) + "\n")
        i += 1
    file_wt.close()
