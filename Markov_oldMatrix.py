# -*- coding: utf-8 -*-
"""
Created on Thu Nov 29 11:05:56 2018

@author: ilanaarbisser
"""

#==============================================================================
# Using old code
#==============================================================================

import os
import math
import numpy as np
import pandas as pd
from scipy.stats import rv_discrete
import time


def nCr(n,r):
    f = math.factorial
    if n >= r:
        return f(n) / f(r) / f(n-r)
    else: return 0

### Need to mod this 

def isInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False

def isFloat(s):
    try: 
        float(s)
        return True
    except ValueError:
        return False

def isNum(s):
    return isInt(s) or isFloat(s)


#def str_to_eq(eqString):
#    while(eqString.find('nCr') != -1):
#        ncr_start = eqString.find('nCr(')
#        ncr_end = eqString.find(")")
#        [f1, f2] = eqString[(ncr_start+len('nCr(')):ncr_end].split(",")
#        f1_num = str_to_eq(f1)
#        f2_num = str_to_eq(f2)
#        ncr_out = nCr(f1_num, f2_num)
#        eqString = eqString.replace(eqString[(ncr_start):(ncr_end+1)], str(ncr_out))
#    Output = None
#    div_vec = eqString.split('/')
#    if len(div_vec) == 1:
#        terms_vec = eqString.split("+")
#        Sum = 0
#        for term in terms_vec:
#            Prod = 1
#            for factor in term.split("*"):
#                if isNum(factor):
#                    Prod *= float(factor)
#                else: 
#                    print "Problem with "+factor
#                    break
#            Sum += Prod
#        Output = float(Sum)
#    elif len(div_vec) == 2:
#        Output = float(str_to_eq(div_vec[0]))/float(str_to_eq(div_vec[1]))
#    else:
#        print 'Problem with input: more than one divisor'
#        return None
#    return Output

def str_to_eq(eqString):
    while(eqString.find('nCr') != -1):
        ncr_start = eqString.find('nCr(')
        ncr_end = eqString.find(")")
        [f1, f2] = eqString[(ncr_start+len('nCr(')):ncr_end].split(",")
        f1_num = str_to_eq(f1)
        f2_num = str_to_eq(f2)
        ncr_out = nCr(f1_num, f2_num)
        eqString = eqString.replace(eqString[(ncr_start):(ncr_end+1)], str(ncr_out))
    return float(pd.eval(eqString))

### Make new matrices

pmat_R01M01 = pd.read_csv('Pmat_R01_M01.csv', header = 0, index_col = 0)
pmat_R1M01 = pd.read_csv('Pmat_R1_M01.csv', header = 0, index_col = 0)
pmat_R10M01 = pd.read_csv('Pmat_R10_M01.csv', header = 0, index_col = 0)

pmat_R01M1 = pd.read_csv('Pmat_R01_M1.csv', header = 0, index_col = 0)
pmat_R1M1 = pd.read_csv('Pmat_R1_M1.csv', header = 0, index_col = 0)
pmat_R10M1 = pd.read_csv('Pmat_R10_M1.csv', header = 0, index_col = 0)

tmat_R01M01 = np.vectorize(str_to_eq)(pmat_R01M01.values)
tmat_R1M01 = np.vectorize(str_to_eq)(pmat_R1M01.values)
tmat_R10M01 = np.vectorize(str_to_eq)(pmat_R10M01.values)

tmat_R01M1 = np.vectorize(str_to_eq)(pmat_R01M1.values)
tmat_R1M1 = np.vectorize(str_to_eq)(pmat_R1M1.values)
tmat_R10M1 = np.vectorize(str_to_eq)(pmat_R10M1.values)

#tmat = tmat_R10M1
#for i in range(np.shape(tmat)[0]):
#    if not np.isclose(sum(tmat[i,]), 1):
#        print "At row "+str(i)+", sum is "+str(sum(tmat[i,]))

pmat = pmat_R1M1

old_eqht_states = ['p.10', 'p.1A', 'p.1B', 'p.81', 'p.82', 'p.8A', 'p.8B']
old_uneqht_states = ['p.41', 'p.42', 'p.4A', 'p.4B', 'p.51', 'p.52', 'p.5A', 'p.5B']

old_eqht_num = np.zeros(len(old_eqht_states))
for i in range(len(old_eqht_states)):
    state = old_eqht_states[i]
    if state in list(pmat): old_eqht_num[i] = list(pmat).index(state)
    else: 
        print "Problem with " + state
        break

old_uneqht_num = np.zeros(len(old_uneqht_states))
for i in range(len(old_uneqht_states)):
    state = old_uneqht_states[i]
    if state in list(pmat): old_uneqht_num[i] = list(pmat).index(state)
    else: 
        print "Problem with " + state
        break

old_eqht_num = old_eqht_num.astype(int)
old_uneqht_num = old_uneqht_num.astype(int)

list(pmat)[0]
list(pmat)[44]
list(pmat)[45]

def run_process(tmat):
    start = 0
    end = [44, 45]
    steps = 0
    total_steps = 1000
    i = start
    prev = None
    while i not in end and steps < total_steps:
        row = tmat[i,:]
        if not np.isclose(sum(tmat[i, :]), 1.0):
            print "Not equal to 1 at " + i
            break
        else:
            j = rv_discrete(values=(range(len(row)),row)).rvs(size=1)[0]
            prev = i
            i = j
            steps = steps + 1    
    return(prev)

def run_fullProcess(tmat):
    start = 0
    end = [44, 45]
    steps = 0
    total_steps = 1000
    path = np.zeros(total_steps)
    i = start
    path[0] = i
    while i not in end and steps < total_steps:
        row = tmat[i,:]
        if not np.isclose(sum(tmat[i, :]), 1.0):
            print "Not equal to 1 at " + i
            break
        else:
            j = rv_discrete(values=(range(len(row)),row)).rvs(size=1)[0]
            i = j
            steps = steps + 1    
            path[steps] = int(i)
    return path[0:(steps+1)].astype(int)

def eqht_fullpath(path):
    ## m_tag is boolean flag for if migration included or not
    pathlen = len(path)
    if path[len(path)-1] not in [44, 45]:
        print "Problem: last state was " + str(path[len(path)-1])
        return None
    elif path[pathlen-2] in old_eqht_num: return True
    elif path[pathlen-2] in old_uneqht_num: return False
    else:
        print "Problem: penultimate state was "+str(path[len(path)-2])
        return None

def equal_heights(num):
    if num in old_eqht_num: 
        return True
    elif num in old_uneqht_num:
        return False
    else: 
        print "Problem with index: "+str(num)        
        return None

R_vals = [0.1, 1, 10]

#old_sims_M1 = np.zeros(3)
old_sims_M01 = np.zeros(3)

R_index = 2
tmat = tmat_R10M01

num_sims = 10000
sims = 0
time1 = time.time()
for i in range(num_sims):
    sims += equal_heights(run_process(tmat))
old_sims_M01[R_index] = float(sims)/num_sims
time.time()-time1
time1 = time.time()

# old_sims_M1[0] = 0.79339999999999999 # M = 1, R = 0.1 1.31680881579717
# old_sims_M1[1] = 0.5645 # M = 1, R = 1, 2.8037676334381105
# old_sims_M1[2] = 0.90910000000000002 # M = 1, R = 10, 0.79339999999999999 # 4.513932398955027

# sims_m1[2] = 0.77910000000000001 # 78.71453285217285
# sims_m1[3] = 0.55740000000000001 #
# sims_m1[4] = 0.90200000000000002 # 4.4523774822553

# old_sims_M01[0] = 0.74180000000000001 # 0.7692142963409424
# old_sims_M01[1] = 0.51739999999999997 # 2.659664018948873
# old_sims_M01[2] = 0.90190000000000003 # 5.8513829986254375

#sims_m01[2] = 0.73919999999999997 # 0.8481764674186707
#sims_m01[3] = 0.51119999999999999 # 2.961593035856883
#sims_m01[4] = 0.90590000000000004 # 4.991067500909169

tmat = tmat_R10M01
tmat_n = tmat ** 20
sum(sum(tmat_n[:, old_eqht_num]))
sum(sum(tmat_n[:, old_uneqht_num]))








