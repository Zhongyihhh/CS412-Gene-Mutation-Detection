# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 00:55:08 2019

"""

# This project use Gibbs Sampling to find the Motif

import numpy as np
from Bio import motifs
# Tutorial: https://biopython-cn.readthedocs.io/zh_CN/latest/en/chr14.html#sec:links
from Bio.Seq import Seq
import os
import time

# Process the Sequence.fa file to read sequences

def input_SQ(loc_sq):
    f = open(loc_sq,'r')
    SQ_all = [] # all sequences of 10 trials
    for i in range(10):
        SQ_one = []
        f1 = f.readline()
        while (f1 != '\n'):
            if len(f1) != 500:
                f1 = f.readline()
            SQ_one += [f1.rstrip('\n')]
            f1 = f.readline()
            if f1 == '':
                break

        SQ_all += [SQ_one]
    SC = len(SQ_all[0])
    return SQ_all, SC

def input_ML(loc_ml):
    f = open(loc_ml,'r')
    ML = int(f.readline())
    return ML


# Calculate the Information Content given one probility
def cal_IC(pos):
    if pos >= 0.0001:
        y = pos*np.log2(pos/0.25)
    else:
        y = 0
    return y


def cal_total_IC(starting, SC, ML, SQ):
    ICsum = 0
    instances = []
    for i in range(SC):
        instances += [Seq(SQ[i][starting[i]:starting[i]+ML])]
    temp = (motifs.create(instances)).counts
    mat = temp.normalize(pseudocounts=0)
    for j in list(mat.values()):
        for jj in j:
            ICsum += cal_IC(jj)
    return ICsum

### Main function to find the motif
def prediction(SQ, ML, SC): # SQ: sequences of DNA in list, ML: Motif Length

    flag = 0
    still_ite = 0
    best_pst = []
    best_IC = 0

    num_record = SC  # Number of sequences to calculate the Residue Sum of Square

    SC = len(SQ) # squence count: number of sequences
    SL = len((' '.join(SQ[0])).split()) # squence length

    pst = [] # following determines the start points of each sequence
    for i in range(SC):
       pst.append(np.random.choice(SL-ML+1)) # set the initial position for each sequence

    ## iteration to update the starting point until convergence

    ite = 1 # iteration count
    sqnum = -1 # The selected sequence to change the starting point

    phase_change = 40
    SSR = 10000 # Initial Residue Sum of Square
    record = {}
    dist_mat = {}
    while (SSR >= 1 and ite <= 20000):

        r = 10 # Parameters to search the potential jumpping point

        # Phase Shifts: every 40 steps, change whole starting point once "Jumping"
        record_IC = {}
        if (ite%(phase_change) == 0):

            # Reocrd the best pst
            temp_IC = cal_total_IC(pst, SC, ML, SQ)
            if temp_IC > best_IC:
                best_IC = temp_IC
                best_pst = pst.copy()

            for a in range(0, r+1):
                instances = []
                movement = a - r/2
                sum_IC = 0
                temp = [int(i + movement) for i in pst]
                if min(temp) >= 0 and max(temp) <= SL - ML: # Examine if this movement is feasible
                    for ii in range(SC):
                        instances += [Seq(SQ[ii][temp[ii]:temp[ii]+ML])]
                    temp_motif = (motifs.create(instances)).counts
                    weight_mat_motif = temp_motif.normalize(pseudocounts=0.1)
                    for j in list(weight_mat_motif.values()):
                        for jj in j:
                            sum_IC += cal_IC(jj) # Calculate the Information Content for each potential new starting points
                    record_IC[movement] = sum_IC
                else:
                    record_IC[movement] = 0
                    continue
            list_IC = list(record_IC.values())
            temp_pst = list(record_IC.keys())
            IC_norm = [item/sum(list_IC) for item in list_IC] # Normalize IC to get the probability
            selected_movement = int(np.random.choice(temp_pst,1,p = IC_norm)) # Sample new starting points according to the probabilities derived from IC
            pst = [int(i + selected_movement) for i in pst] # Update the starting points of all sequences

        sqnum_new = np.random.choice(range(SC)) # sample one sequence to predict startpoint position and the others are used to generate current motif table
        sqnum = sqnum_new
        instances = []
        for i in range(SC):
            if i == sqnum:
                selected_sq = (' '.join(SQ[i])).split()
            else:
                instances += [Seq(SQ[i][pst[i]:pst[i]+ML])]
        m = (motifs.create(instances)).counts
        pwm = m.normalize(pseudocounts=0.1) # Obtain position-weight matrix

        Ax = [] #startpoint weight
        for i in range(SL-ML+1):
            temp = 1
            j = 0
            motif_try = selected_sq[i:i+ML]
            for item in motif_try:
                temp = temp*pwm[item,j]/0.25
                j += 1
            Ax.append(temp)
        Ax_norm = [item/sum(Ax) for item in Ax] # Normalize the starting point weight matrix
        temp1 = int(np.random.choice(range(SL-ML+1),1,p = Ax_norm)) # Get a new start point position for the selected sequence
        pst[sqnum] = temp1

        record[ite] = pst.copy()

        # Calculate SSR of "r" consecutive iterations based on the starting points
        if ite <= 4000 + still_ite:  # Force the algorithm to run at least 1000 times
            SSR = 1000000
        elif flag == 0:    # Conduct test to see how many iterations still need to run according to the IC
            current_IC = best_IC;
            still_ite = 1000*((2*ML - current_IC)/2)**2
            flag = 1
        else:
            del record[list(record.keys())[0]]
            new_dist_mat = {}
            new_dist_list = []
            for m in range(num_record-1):
                mm = list(record.keys())[m]
                for n in range(m + 1, num_record):
                    nn = list(record.keys())[n]
                    try:
                        new_dist_mat[str(mm) + str(nn)] = dist_mat[str(mm) + str(nn)]
                        new_dist_mat[str(nn) + str(mm)] = dist_mat[str(mm) + str(nn)]
                        new_dist_list.append(dist_mat[str(mm) + str(nn)])
                    except:
                        temp = sum([(a - b) ** 2 for a, b in zip(record[mm], record[nn])])
                        new_dist_mat[str(mm) + str(nn)] = temp
                        new_dist_mat[str(nn) + str(mm)] = temp
                        new_dist_list.append(temp)
            dist_mat = new_dist_mat
            SSR = sum(new_dist_list)
        ite += 1


    # conduct final test to see which pst has the largest IC
    temp_IC = cal_total_IC(pst, SC, ML, SQ)
    if temp_IC < best_IC:
        pst = best_pst.copy()


    # At the final step, do one last movement to find ture optimal

    record_IC = {}
    for a in range(0, r+1):
        instances = []
        movement = a - r/2
        sum_IC = 0
        temp = [int(i + movement) for i in pst]
        if min(temp) >= 0 and max(temp) <= SL - ML:
            for ii in range(SC):
                instances += [Seq(SQ[ii][temp[ii]:temp[ii]+ML])]
            temp_motif = (motifs.create(instances)).counts
            weight_mat_motif = temp_motif.normalize(pseudocounts=0.1)
            for j in list(weight_mat_motif.values()):
                for jj in j:
                    sum_IC += cal_IC(jj)
            record_IC[movement] = sum_IC
        else:
            record_IC[movement] = 0
            continue
    record_IC_sorted = sorted(record_IC.items(), key = lambda x:x[1], reverse = True)
    selected_movement = record_IC_sorted[0][0]
    pst = [int(i + selected_movement) for i in pst]



    instances = []
    for i in range(SC):
        instances += [Seq(SQ[i][pst[i]:pst[i]+ML])]
    temp_motif = (motifs.create(instances)).counts
    pwm = temp_motif.normalize(pseudocounts=0)

    return pwm,pst,instances

# create predictedsites.txt, predictedmotif.txt files for future steps
def createFile(randMotif, locations, sites, ML, motifName, folder_ext, seq):

    tempChoice = ["A", "C", "G", "T"]

    directory = "data set"
    if folder_ext != "":
        directory = folder_ext

    if not os.path.exists(directory):
        os.makedirs(directory)

    # saves predictedsites.txt file
    if seq > 0:
        file = open(directory + '/predictedsites.txt', 'a')
    else:
        file = open(directory + '/predictedsites.txt', 'w')

    if seq > 0:
        file.write("\n")
    for i in range(len(sites)):
        file.write(str(locations[i]) + ", " + sites[i] + "\n")
    file.close

    # saves motif.txt file
    if seq > 0:
        file = open(directory + '/predictedmotif.txt', 'a')
    else:
        file = open(directory + '/predictedmotif.txt', 'w')
    if seq > 0:
        file.write("\n")
    file.write(">" + motifName + "_" + str(seq) + "\t" + str(ML) + "\n")
    for i in range(ML):
        itemStr = ""
        for item in tempChoice:
            itemStr += str(randMotif[item, i]) + "\t"

        file.write(itemStr.strip() + "\n")
    file.write("<")
    file.close




############################## Main Script#############################
#ICPC = 2
#ML = 7
#SC = 10

def run_Step2(in_ICPC, in_ML, in_SC):

    folder_ext = str(in_ICPC) + "_" + str(in_ML) + "_" + str(in_SC)
    loc_sq = 'sequences.fa'
    SQ_all, SC= input_SQ(folder_ext + "/" + loc_sq)
    loc_ml = 'motiflength.txt'
    ML = input_ML(folder_ext + "/" + loc_ml)
    motifName = "MOTIF1"

    i = 0
    runtimes = []
    for SQ in SQ_all:
        start = time.time()
        pwm,pst,instances = prediction(SQ, ML, SC)
        runtimes.append(time.time() - start)
        predictedSites = []
        for site in instances:
            predictedSites.append(str(site))

        createFile(pwm, pst, predictedSites, ML, motifName + "_" + folder_ext, folder_ext, i)
        i += 1
    print(folder_ext, sum(runtimes)/10)
