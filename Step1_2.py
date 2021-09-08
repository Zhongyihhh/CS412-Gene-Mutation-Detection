"""
Created on Wed Apr 17 15:12:32 2019

"""

import random
from random import randint
from random import choices
from Step1_3 import motif_results
import os

tempChoice = ["A", "C", "G", "T"]

# step 1.2 function to create sequence
def createSequences(SC, SL):

    tempChoice = ["A", "C", "G", "T"]

    sequeuces = []

    textCntMax = int(SL / 4) + 5

    for i in range(SC):
        tempStr = ""
        backMotif = [0.25] * 4
        textCount = [0] * 4
        for j in range(SL):
            textChosen = choices(tempChoice, backMotif)[0]
            tempStr += textChosen
            textCount[tempChoice.index(textChosen)] += 1

            if textCount[tempChoice.index(textChosen)] >= textCntMax:
                backMotif = [0] * 4
                textCnt = 0
                for cnt in textCount:
                    if cnt < textCntMax:
                        textCnt += 1
                for k in range(4):
                    if textCount[k] < textCntMax:
                        backMotif[k] = 1.0 / textCnt

        sequeuces.append(tempStr)
        random.shuffle(tempChoice)

    return sequeuces

# create site using random motif
def createSite(SC, ML, randMotif):

    tempChoice = ["A", "C", "G", "T"]

    sites = []
    for i in range(SC):
        tempStr = ""
        for j in range(ML):
            tempStr += choices(tempChoice, randMotif[j])[0]
        sites.append(tempStr)
    return sites

# plant each site into each sequence
def plant(SC, ML, RS, sites):
    sequences = []
    locations = []
    for i in range(SC):
        sequence = RS[i]
        location = randint(0, len(sequence) - ML)
        sequence = sequence[0: location] + sites[i] + sequence[location + ML:]
        sequences.append(sequence)
        locations.append(location)
    return sequences, locations

# create sequences.fa, sites.txt, motif.txt, motiflength.txt files for future steps
def createFile(sequences, randMotif, sites, locations, ML, motifName, folder_ext, index):
    directory = "data set"
    if folder_ext != "":
        directory = folder_ext

    if not os.path.exists(directory):
        os.makedirs(directory)

    # saves sequences.fa file
    if index > 0:
        file = open(directory + '/sequences.fa','a')
    else:
        file = open(directory + '/sequences.fa','w')
    i = 0

    if index > 0:
        file.write("\n")

    for sequence in sequences:
        file.write('>seq' + str(i) + "\n")
        file.write(sequence + "\n")
        i += 1
    file.close()

    # saves sites.txt file
    if index > 0:
        file = open(directory + '/sites.txt', 'a')
    else:
        file = open(directory + '/sites.txt', 'w')

    if index > 0:
        file.write("\n")

    for i in range(len(sites)):
        file.write(str(locations[i]) + ", " + sites[i] + "\n")
    file.close

    # saves motif.txt file
    if index > 0:
        file = open(directory + '/motif.txt', 'a')
    else:
        file = open(directory + '/motif.txt', 'w')


    if index > 0:
        file.write("\n")

    file.write(">" + motifName + "\t" + str(ML) + "\n")
    for motif in randMotif:
        itemStr = ""
        for items in motif:
            itemStr += str(items) + "\t"

        file.write(itemStr.strip() + "\n")
    file.write("<")
    file.close

    # saves motiflength.txt file
    if index > 0:
        file = open(directory + '/motiflength.txt', 'a')
    else:
        file = open(directory + '/motiflength.txt', 'w')

    if index > 0:
        file.write("\n\n")

    file.write(str(ML))
    file.close


def generate_data():
    ICPC = 2
    ML = 8
    SL = 500
    SC = 10

    for i in range(10):
       if i == 1:
           isAppend = True
       motifName = "MOTIF1"

       RS = createSequences(SC, SL)
       randMotif = motif_results(ML, ICPC)
       sites = createSite(SC, ML, randMotif)
       sequences, locations = plant(SC, ML, RS, sites)
       file_ext = str(ICPC) + "_" + str(ML) + "_" + str(SC)
       createFile(sequences, randMotif, sites, locations, ML, motifName + "_" + str(i), file_ext, i)

       ICPCSet = [1, 1.5]
       MLSet = [6, 7]
       SCSet = [5, 20]

       # create 2 files for ICPC changes
       for ICPC in ICPCSet:
           file_ext = str(ICPC) + "_" + str(ML) + "_" + str(SC)

           RS = createSequences(SC, SL)
           randMotif = motif_results(ML, ICPC)
           sites = createSite(SC, ML, randMotif)
           sequences, locations = plant(SC, ML, RS, sites)
           createFile(sequences, randMotif, sites, locations, ML, motifName + "_" + file_ext + "_" + str(i), file_ext, i)

       ICPC = 2
       # create 2 files for ML changes
       for ML in MLSet:
           file_ext = str(ICPC) + "_" + str(ML) + "_" + str(SC)

           RS = createSequences(SC, SL)
           randMotif = motif_results(ML, ICPC)
           sites = createSite(SC, ML, randMotif)
           sequences, locations = plant(SC, ML, RS, sites)
           createFile(sequences, randMotif, sites, locations, ML, motifName + "_" + file_ext + "_" + str(i), file_ext, i)

       ML = 8
       # create 2 files for SC changes
       for SC in SCSet:
           motifName = "MOTIF1"
           file_ext = str(ICPC) + "_" + str(ML) + "_" + str(SC)

           RS = createSequences(SC, SL)
           randMotif = motif_results(ML, ICPC)
           sites = createSite(SC, ML, randMotif)
           sequences, locations = plant(SC, ML, RS, sites)
           createFile(sequences, randMotif, sites, locations, ML, motifName + "_" + file_ext + "_" + str(i), file_ext, i)
       SC = 10

generate_data()
