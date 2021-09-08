# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 14:19:35 2019

@author: yinin
"""

from Step1_2 import generate_data
from Step2_Function import run_Step2
from Step3_Function import run_Evaluation

# Generate dataset
generate_data()
print ('Data Generation is finished')
all_Combinations = [[1,8,10],[1.5,8,10],[2,8,10],[2,6,10],[2,7,10],[2,8,5],[2,8,20]]

# For each dataset, find motif
for i in all_Combinations:
    temp_ICPC = i[0]
    temp_ML = i[1]
    temp_SC = i[2]
    run_Step2(temp_ICPC, temp_ML, temp_SC)  ## These values are not used to find motif, only used to read the folders
    print (str(temp_ICPC) + '_' + str(temp_ML) + '_' + str(temp_SC) + ' is finished')

# Do evaluation
run_Evaluation()