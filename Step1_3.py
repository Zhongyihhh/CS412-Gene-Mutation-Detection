# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 16:37:39 2019

"""


import numpy as np

def cal_ICPC(pos):
    if pos >= 0.0001:
        y = pos*np.log2(pos/0.25)
    else:
        y = 0
    return y

def find_possible(ICPC):
    possible_range = []
    for ii in range(0,10001):
        i = ii/10000
        if ii == 10000:
            i = 1
        if ii == 0:
            i = 0
        max_ICPC = cal_ICPC(i)+cal_ICPC(1-i)
        min_ICPC = cal_ICPC(i)+cal_ICPC((1-i)/3)*3
        if max_ICPC >= ICPC - 0.00001 and min_ICPC <= ICPC + 0.00001:
            possible_range.append(i)
    first = np.random.choice(possible_range)

    possible_range_2nd = []
    for iii in range(0, 10001 - int(np.round(first*10000))):
        i1 = iii/10000;
        max_ICPC2 = cal_ICPC(i1) + cal_ICPC(1-first-i1) + cal_ICPC(first)
        min_ICPC2 = cal_ICPC(i1)+cal_ICPC((1-first-i1)/2)*2 + cal_ICPC(first)
        if max_ICPC2 >= ICPC - 0.0001 and min_ICPC2 <= ICPC + 0.0001:
            possible_range_2nd.append(i1)
    second = np.random.choice(possible_range_2nd)

    if second == 0.0:
        second = 0
    if second == 1.0:
        second = 1

    remain_total = 1 - first - second
    if remain_total == 1.0:
        remain_total = 1
    if remain_total == 0.0:
        remain_total = 0
    remain_ICPC = ICPC - cal_ICPC(first) - cal_ICPC(second)

    for kk in range(0, int(np.round(remain_total*10000)+1)):
        k = kk/10000
        if k == 0.0:
            k = 0
        if k == 1.0:
            k = 1
        if cal_ICPC(k) + cal_ICPC(remain_total - k) <= remain_ICPC + 0.001 and cal_ICPC(k) + cal_ICPC(remain_total - k) >= remain_ICPC - 0.001:
            third = k
            forth = np.round(10000*(remain_total - third))/10000
            if forth == 0.0:
                forth = 0
            if forth == 1.0:
                forth = 1
            break
    switch = [third,forth]
    t = np.random.choice([0,1])
    third = switch[t]
    forth = switch[int(1-t)]
    return first,second,third,forth


def motif_results(ML, ICPC):
    motif_table = []
    for i in range(0, ML):
        this = find_possible(ICPC)
        motif_table.append(this)
    return motif_table
