# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 10:10:31 2018

@author: ilanaarbisser
"""

#==============================================================================
# Basic functions
#==============================================================================


## Define state space

import os
import math
import numpy as np
import pandas as pd
from scipy.stats import rv_discrete
import time

#os.chdir("/Users/ilanaarbisser/Dropbox/Markov_ReMig")

# n samples
# m loci
# p populations

## Define states
## 1 leq i + k leq n
## 1 leq j + k leq n
## 0 leq i, j, k

def nCr(n,r):
    f = math.factorial
    if n >= r:
        return f(n) / f(r) / f(n-r)
    else: return 0

#def isinbounds(vec, n):
#    i = vec[0]
#    j = vec[1]
#    k = vec[2]
#    return 0 < (i + k) and (i + k) <= n and 0 < (j + k) and (j + k) <= n

def vectoString(vec):
    if vec is not None:
        if len(vec) == 3:
            [i, j, k] = vec
            return str(i) + str(j)+ str(k)
        elif len(vec) == 6:
            [iA, jA, kA, iB, jB, kB] = vec
            return '{'+str(iA)+str(jA)+str(kA)+ '}{'+str(iB)+str(jB)+str(kB)+'}'
    else: return None
    
def stringtoVec(state):
    if state is not None:
        if len(state) == 3:
            return [int(state[0]),int(state[1]),int(state[2])]
        elif len(state) == 10:
            return [int(state[1]),int(state[2]),int(state[3]),int(state[6]),int(state[7]),int(state[8])]
    else: return None
    
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


def string_add(s1, s2):
    addsum = ""
    if isInt(s1) and isInt(s2):
        addsum = str(int(s1)+int(s2))
    elif not isInt(s1) :
        if s1 == s2:
            addsum = str(2)+s1        
        elif s2 == "0":
                addsum = s1
        else: addsum = s2 + " + " + s1
    elif s1 == "0":
        addsum = s2
    else: addsum = s1 + " + " + s2
    return addsum

def string_times(s1, s2):
    product = ""
    if isInt(s1) and isInt(s2):
        product = str(int(s1)*int(s2))
    elif not isInt(s1) :
        if isInt(s2):
            if s2 == "0":
                product = "0"
            elif s2 == 1:
                product = s1
            else: product = s2 + s1
        else: product = s1 + "*" + s2
    elif int(s1) == 0:
        product = "0"
    elif int(s1) == 1:
        product = s2
    return product
    
def find_state_space_size(n):
    return n*(2*(n ** 2) + 9*n + 1)/6

def find_states(n):
    num_states = find_state_space_size(n)
    states = np.reshape(np.zeros(num_states*3, dtype = int), (num_states, 3))
    np.shape(states)
    s_idx = 0
    for i in range(n+1):
        for j in range(n+1):
            for k in range(n+1):
                if 0 < (i + k) and (i + k) <= n and 0 < (j + k) and (j + k) <= n:
                    states[s_idx, 0] = i
                    states[s_idx, 1] = j
                    states[s_idx, 2] = k
                    s_idx = s_idx + 1
    return states

def find_stateNames(states):
    state_names = [""]*np.shape(states)[0]
    for i in range(np.shape(states)[0]):
        state_names[i] = str(int(states[i, 0])) + str(int(states[i, 1])) + str(int(states[i, 2]))
    return state_names


#"Re": [1, 1, -1]
#"C1_ell": [-1, 0, 0]
#"C1_r": [0, -1, 0]
#"C2": [-1, -1, 1]
#"C3_ell": [-1, 0, 0]
#"C3_r": [0, -1, 0]
#"C4": [0, 0, -1]}

#if isinbounds(end, n):        
#            return end
#        else: return None

#def Re(start):
#    if np.equal(start, [0, 0, 1]).all():
#        return None
#    elif start[2] > 0:
#        end = np.add(start, [1, 1, -1])
#        return end
#    else: return None
#
#
#def C1(start, side):
#    if np.equal(start, [0, 0, 1]).all():
#        return None
#    elif nCr(start[side], 2) > 0:
#        move = np.zeros(3, int)
#        move[side] = -1
#        end = np.add(start, move)
#        return end
#    else: return None
#
#def C2(start):
#    if np.equal(start, [0, 0, 1]).all():
#        return None
#    elif start[0]*start[1] > 0:
#        end = np.add(start, [-1, -1, 1])
#        return end
#    else: return None
#    
#def C3(start, side):
#    if np.equal(start, [0, 0, 1]).all():
#        return None
#    elif start[side]*start[2] > 0:
#        move = np.zeros(3, int)
#        move[side] = -1
#        end = np.add(start, move)
#        return end
#    else: return None
#
#def C4(start):
#    if np.equal(start, [0, 0, 1]).all():
#        return None
#    elif nCr(start[2], 2) > 0:
#        end = np.add(start, [0, 0, -1])
#        return end
#    else: return None
        
def Re(start):
    if start[2] > 0:
        end = np.add(start, [1, 1, -1])
        return end
    else: return None


def C1(start, side):
    if nCr(start[side], 2) > 0:
        move = np.zeros(3, int)
        move[side] = -1
        end = np.add(start, move)
        return end
    else: return None

def C2(start):
    if start[0]*start[1] > 0:
        end = np.add(start, [-1, -1, 1])
        return end
    else: return None
    
def C3(start, side):
    if start[side]*start[2] > 0:
        move = np.zeros(3, int)
        move[side] = -1
        end = np.add(start, move)
        return end
    else: return None

def C4(start):
    if nCr(start[2], 2) > 0:
        end = np.add(start, [0, 0, -1])
        return end
    else: return None

def C1L(start):
    return C1(start, 0)
    
def C1R(start):
    return C1(start, 1)

def C3L(start):
    return C3(start, 0)
    
def C3R(start):
    return C3(start, 1)
    
def no_chira(trans):
    if trans[0] is "C":
        trans = (trans.replace("R", "")).replace("L", "")
    return trans

moves = {"Re": Re,
         "C1L": C1L,
         "C1R": C1R,
         "C2": C2,
         "C3L": C3L,
         "C3R": C3R,
         "C4": C4}

def make_transition_df(states):
    num_states = np.shape(states)[0]
    state_names = find_stateNames(states)
    transition_matrix = np.reshape(["0"]*(num_states ** 2), (num_states, num_states))
    transition_df = pd.DataFrame(transition_matrix, index = state_names, columns = state_names)
    for start in states:
        if not np.equal(start, [0, 0, 1]).all():
            for m in moves:
                end = moves[m](start)
                if end is not None:
                    if transition_df[vectoString(end)][vectoString(start)] == "0":
                        transition_df[vectoString(end)][vectoString(start)] = m
                    else: 
                        transition_df[vectoString(end)][vectoString(start)] = transition_df[vectoString(end)][vectoString(start)] + "or" + m
    return transition_df

#==============================================================================
# Graphviz
#==============================================================================

def strip(s1):
    s2 = ''
    if s1.find('C') != -1 or s1.find('Re') != -1 or s1.find('M') != -1:
        s2 = ''        
        if s1.find('C') != -1:
            s2 = s2 + 'C'
        if s1.find('Re') != -1:
            s2 = s2 + 'Re'
        if s1.find('M') != -1:
            s2 = s2 + 'M'
    elif s1 != '0' and s1 != '':
        s2 = 'x'
    else: s2 = '0'
    return s2

def removePopulationLabel(state): 
    vec = [0, 0, 0]
    if state[0] == '{':
        i_total = int(state[1])+int(state[6])
        j_total = int(state[2])+int(state[7])
        k_total = int(state[3])+int(state[8])
        vec = [i_total, j_total, k_total]
    elif state[0] in ['0', '1', '2', '3']:
        return [int(state[0]), int(state[1]), int(state[2])]
    else:
        print "Problem with state"
    return vec

def findStateColor(state):
    color = ''
    ## state_lab from state
    ## Just 2 states
    state_lab = vectoString(removePopulationLabel(state))
    if state_lab == '002':
        color = 'lightcoral'
    elif state_lab == '111':
        color = 'orange'
    elif state_lab == '220':
        color = 'lemonchiffon'
    elif state_lab == '101':
        color = 'darkolivegreen2'
    elif state_lab == '210': 
        color = 'mediumaquamarine'
    elif state_lab == '011':
        color = 'lightblue1'
    elif state_lab == '120':
        color = 'lightskyblue'
    elif state_lab == '110':
        color = 'darkorchid1'
    elif state_lab == '001':
        color = 'lightpink'
    ## additional 11 states for 3 states
    elif state_lab == '003':
        color = 'aquamarine'
    elif state_lab == '012':
        color = 'bisque'
    elif state_lab == '021':
        color = 'blue'
    elif state_lab == '102':
        color = 'blueviolet'
    elif state_lab == '112':
        color = 'brown1'
    elif state_lab == '121':
        color = 'chartreuse'
    elif state_lab == '130':
        color = 'cyan3'
    elif state_lab == '201':
        color = 'hotpink'
    elif state_lab == '211':
        color = 'magenta'
    elif state_lab == '221':
        color = 'lightgoldenrod'
    elif state_lab == '230':
        color = 'palegreen1'
    elif state_lab == '310':
        color = 'plum'
    elif state_lab == '320':
        color = 'salmon'
    elif state_lab == '330':
        color = 'yellow'
    else: 
        print state+" produced no color"        
        color = 'white'
    return color

def findTransColor(trans):
    color = ''
    if strip(trans) == 'C':
        color = 'navy'
    elif strip(trans) == 'Re':
        color = 'darkgreen'
    elif strip(trans) == 'M':
        color = 'red'
    else: color = 'black'
    return color

def make_graphviz(transition_df, filename):
    # dot -T png -O FILENAME.dot
    t_mat = transition_df.as_matrix()
    dot = open(filename, "w")
    dot.write("digraph markov_chain {\n")    
    dot.write("\t"+"rankdir = LR\n")
    dot.write("\t"+"node [shape = circle, style = filled, color = lightskyblue];\n")
    for i in range(np.shape(t_mat)[0]):
            for j in range(np.shape(t_mat)[1]):
                if t_mat[i][j] is not "0":
                    dot.write("\t" + transition_df.index[i] + " -> "+transition_df.index[j]+" [ label = "+strip(t_mat[i][j])+" ];\n")
    dot.write("}")
    dot.close()

def make_graphviz_colors(transition_df, filename):
    # dot -T png -O FILENAME.dot
    dot = open(filename, "w")
    dot.write("digraph markov_chain {\n")    
    dot.write("\t"+"rankdir = LR\n")
    for state in list(transition_df):
        dot.write("\t" + state + "[shape = circle, style = filled, color = "+findStateColor(state)+"];\n")
    t_mat = transition_df.as_matrix()
    for i in range(np.shape(t_mat)[0]):
            for j in range(np.shape(t_mat)[1]):
                if t_mat[i][j] is not "0":
                    dot.write("\t" + transition_df.index[i] + " -> "+transition_df.index[j]+" [ label = "+strip(t_mat[i][j])+", color = "+findTransColor(t_mat[i][j])+"  ];\n")
    dot.write("}")
    dot.close()

n2_df = make_transition_df(find_states(2))
n3_df = make_transition_df(find_states(3))
##
#make_graphviz_colors(n2_df, "n2_recomb_color.dot")
#make_graphviz_colors(n3_df, "n3_recomb_color.dot")
#
##n2_df
##
## Then in terminal run:
## dot -T png -O TEST_n3.dot
#
#n2_df
#n3_df.head()

#==============================================================================
# Adding migration
#==============================================================================

def twovec_tostring(vec1, vec2):
    return "{" + str(vec1[0]) + str(vec1[1]) + str(vec1[2])+"}{" + str(vec2[0]) + str(vec2[1]) + str(vec2[2])+"}"

def mig_states_size(n):
    return int(float(1)/120*(10*n + 329*n**2 + 275*n**3 + 90*n**4 + 15*n**5 + n**6))

def find_mig_states(n, states):
    mig_state_strings = [""]*mig_states_size(n)
    place = 0
    for state in states:
        i = int(state[0])
        j = int(state[1])
        k = int(state[2])
        for ia in range(i+1):
            for ja in range(j+1):
                for ka in range(k+1):
                    mig_state_strings[place] = "{"+str(ia)+str(ja)+str(ka)+"}"+"{"+str(i-ia)+str(j-ja)+str(k-ka)+"}"
                    place += 1
                    #mig_state_strings = filter(None, mig_state_strings)
                    #np.shape(mig_state_strings)
    return mig_state_strings
    
## So as not to double count, only look at lineage++
## In the code block it will look at pop A and B
## Therefore lineageA++ and lineageA-- will be considered

def get_vecA(state):
    return [int(state[1]), int(state[2]), int(state[3])]

def get_vecB(state):
    return [int(state[6]), int(state[7]), int(state[8])]

def is_mig_endState(vecA, vecB):
    is_end1 = np.equal(vecA, [0, 0, 1]).all() and np.equal(vecB, [0, 0, 0]).all()
    is_end2 = np.equal(vecB, [0, 0, 1]).all() and np.equal(vecA, [0, 0, 0]).all()
    return is_end1 or is_end2

def MigtoB(startvecA, startvecB, idx):
    if is_mig_endState(startvecA, startvecB):
        return None
    elif startvecA[idx] > 0:
        add_vec = np.zeros(3, int)
        add_vec[idx] = 1
        sub_vec = [x*(-1) for x in add_vec]
        endA = np.add(startvecA, sub_vec)
        endB = np.add(startvecB, add_vec)
        return "{"+vectoString(endA)+"}{"+vectoString(endB)+"}"
    else: return None

def MigtoA(startvecA, startvecB, idx):
    if is_mig_endState(startvecA, startvecB):
        return None
    elif startvecB[idx] > 0:
        add_vec = np.zeros(3, int)
        add_vec[idx] = 1
        sub_vec = [x*(-1) for x in add_vec]
        endB = np.add(startvecB, sub_vec)
        endA = np.add(startvecA, add_vec)
        return "{"+vectoString(endA)+"}{"+vectoString(endB)+"}"
    else: return None

def MitoB(start):
    return MigtoB(get_vecA(start), get_vecB(start), 0)    

def MitoA(start):
    return MigtoA(get_vecA(start), get_vecB(start), 0)    

def MjtoB(start):
    return MigtoB(get_vecA(start), get_vecB(start), 1)    

def MjtoA(start):
    return MigtoA(get_vecA(start), get_vecB(start), 1)    

def MktoB(start):
    return MigtoB(get_vecA(start), get_vecB(start), 2)    

def MktoA(start):
    return MigtoA(get_vecA(start), get_vecB(start), 2)    

migMoves = {"MitoB": MitoB, "MitoA": MitoA, "MjtoB": MjtoB, "MjtoA": MjtoA, "MktoB": MktoB, "MktoA": MktoA}

#==============================================================================
#TESTING 
#==============================================================================

#mig_states = find_mig_states(2, list(n2_df))
#
#num_states = np.shape(mig_states)[0]
#t_matrix = np.reshape(["0"]*(num_states ** 2), (num_states, num_states))
#t_df = pd.DataFrame(t_matrix, index = mig_states, columns = mig_states) 
#
#t_matrix = np.reshape(["0"]*(num_states ** 2), (num_states, num_states))
#t_df = pd.DataFrame(t_matrix, index = mig_states, columns = mig_states) 
#
#def test_run(start, params):
#    dic = make_dictionary(np.vectorize(float)(stringtoVec(start)+params))
#    for m in moves:
#        ## Possible coalescence/re in pop A
#        end_vec = moves[m](get_vecA(start))
#        if end_vec is not None:
#            end_string = "{"+vectoString(end_vec)+"}{"+vectoString(get_vecB(start))+"}"
#            if end_string in t_df.index and start in t_df.index:
#                print m+'A'
#                print process_prob_name(m+'A', dic)
#            else: print "start: " + start + " " + m +" to end: " + end_string
#        ## Possible coalescence/re in pop B        
#        end_vec = moves[m](get_vecB(start))
#        if end_vec is not None:
#            end_string = "{"+vectoString(get_vecA(start))+"}{"+vectoString(end_vec)+"}"
#            if end_string in t_df.index and start in t_df.index:
#                print m+'B'
#                print process_prob_name(m+'B', dic)
#            else: print "start: " + start + " " + m +" to end: " + end_string
#    for mig in migMoves:
#        end_string = migMoves[mig](start)
#        if end_string is not None:
#            if end_string in t_df.index and start in t_df.index:
#                print mig
#                print process_prob_name(mig, dic)
#            else: print "start: " + start + " " + m +" to end: " + end_string
#
#start = '{010}{101}'
#test_run(start, [1.5, 0.1])

# probs_mig[name]
#{000}{001} ok
#{001}{000} ok
#{000}{002} ok
#{001}{001} ok
#{002}{000} ok
#{000}{011} ok
#{001}{010} ok
#{010}{001} ok
#{011}{000} ok
#{000}{101} ok
#{001}{100} ok
#{100}{001} ok
#{101}{000} ok
#{000}{110} ok
#{010}{100} ok
#{100}{010} ok
#{110}{000} ok
#{000}{111} ok 
#{001}{110} ok
#{010}{101}
#{011}{100}
#{100}{011}
#{101}{010}
#{110}{001}
#{111}{000}
#{000}{120}
#{010}{110}
#{020}{100}
#{100}{020}
#{110}{010}
#{120}{000}
#{000}{210}
#{010}{200}
#{100}{110}
#{110}{100}
#{200}{010}
#{210}{000}
#{000}{220}
#{010}{210}
#{020}{200}
#{100}{120}
#{110}{110}
#{120}{100}
#{200}{020}
#{210}{010}
#{220}{000}


#==============================================================================
# TESTING
#==============================================================================

def make_migT_df(mig_states):
    num_states = np.shape(mig_states)[0]
    t_matrix = np.reshape(["0"]*(num_states ** 2), (num_states, num_states))
    t_df = pd.DataFrame(t_matrix, index = mig_states, columns = mig_states) 
    for start in t_df.index:  
        for m in moves:
            ## Possible coalescence/re in pop A
            end_vec = moves[m](get_vecA(start))
            if end_vec is not None:
                end_string = "{"+vectoString(end_vec)+"}{"+vectoString(get_vecB(start))+"}"
                if end_string in t_df.index and start in t_df.index:
                    if t_df[end_string][start] == "0":
                        t_df[end_string][start] = m+"A"
                    #elif  t_df[end_string][start] is not m: 
                        #t_df[end_string][start] = t_df[end_string][start] + "or" + m+"A"
                else: print "start: " + start + " " + m +" to end: " + end_string
            ## Possible coalescence/re in pop B        
            end_vec = moves[m](get_vecB(start))
            if end_vec is not None:
                end_string = "{"+vectoString(get_vecA(start))+"}{"+vectoString(end_vec)+"}"
                if end_string in t_df.index and start in t_df.index:
                    if t_df[end_string][start] == "0":
                        t_df[end_string][start] = m+"B"
                    #elif t_df[end_string][start] is not m: 
                        #t_df[end_string][start] = t_df[end_string][start] + "or" + m+"B"
                else: print "start: " + start + " " + m +" to end: " + end_string
        ## Migration
        for mig in migMoves:
            end_string = migMoves[mig](start)
            if end_string is not None:
                if end_string in t_df.index and start in t_df.index:
                    if t_df[end_string][start] == "0":
                        t_df[end_string][start] = mig
                    #elif  t_df[end_string][start] is not mig:
                        #t_df[end_string][start] = t_df[end_string][start] + "or" + mig
                else: print "start: " + start + " " + m +" to end: " + end_string
    return t_df

def gviz_name(state):
    return "A"+vectoString(get_vecA(state))+"B"+vectoString(get_vecB(state))    
    
def make_graphviz_mig(transition_df, filename):
    # dot -T png -O FILENAME.dot
    t_mat = transition_df.as_matrix()
    dot = open(filename, "w")
    dot.write("digraph markov_chain {\n")    
    dot.write("\t"+"rankdir = LR\n")
    dot.write("\t"+"node [shape = circle, style = filled, color = lightskyblue];\n")
    for i in range(np.shape(t_mat)[0]):
            for j in range(np.shape(t_mat)[1]):
                if t_mat[i][j] is not "0":
                    #dot.write("\t" + gviz_name(transition_df.index[i]) + " -> "+gviz_name(transition_df.index[j])+" [ label = "+t_mat[i][j]+" ];\n")
                    dot.write("\t" + gviz_name(transition_df.index[i]) + " -> "+gviz_name(transition_df.index[j])+";\n")
    dot.write("}")
    dot.close()

def make_graphviz_mig_color(transition_df, filename):
    # dot -T png -O FILENAME.dot
    t_mat = transition_df.as_matrix()
    dot = open(filename, "w")
    dot.write("digraph markov_chain {\n")    
    dot.write("\t"+"rankdir = LR\n")
    for state in list(transition_df):
        dot.write("\t"+gviz_name(state)+" [shape = circle, style = filled, color = "+findStateColor(state)+"];\n")    
    #dot.write("\t"+"node [shape = circle, style = filled, color = lightskyblue];\n")
    for i in range(np.shape(t_mat)[0]):
            for j in range(np.shape(t_mat)[1]):
                if t_mat[i][j] is not "0":
                    #dot.write("\t" + gviz_name(transition_df.index[i]) + " -> "+gviz_name(transition_df.index[j])+" [ label = "+t_mat[i][j]+" ];\n")
                    dot.write("\t" + gviz_name(transition_df.index[i]) + " -> "+gviz_name(transition_df.index[j])+" [ color = "+findTransColor(t_mat[i][j])+" ];\n")
    dot.write("}")
    dot.close()

n = 2
states = find_states(n)
state_strings = find_stateNames(states)
mig_states = find_mig_states(n, states)
n2_mig_df = make_migT_df(mig_states)

n = 3
states = find_states(n)
state_strings = find_stateNames(states)
mig_states = find_mig_states(n, states)
n3_mig_df = make_migT_df(mig_states)

#make_graphviz_mig_color(n2_mig_df, 'n2_mig_color.dot')
#make_graphviz_mig_color(n3_mig_df, 'n3_mig_color.dot')

#==============================================================================
# Different tree heights for n = 2
#==============================================================================
#count = 0
#for state in list(n3_mig_df):
#    if n3_mig_df.loc[state, '{000}{001}'] != '0' or n3_mig_df.loc[state, '{001}{000}'] != '0':
#        count += 1   
#        print count        
#        print state
                
# 002 equal
# 011 equal
# 101 unequal
# 110 unequal

## {000}{002}
## {002}{000}
## {000}{011}
## {011}{000}
## {000}{101}
## {101}{000}
## {000}{110}
## {110}{000}
## Same for n = 2 and n = 3, also the same 
      

## Turn diagram matrix into a transition matrix with real values
## Create a function that takes in R and M as parameters
## Need to convert letter matrix to number matrix
## Then also need to repalce C, R, and M with real values
## How do??

probs = {"Re": "(k*R)/(nCr(i+j+k, 2)+k*R)",
         "C1L": "(nCr(i, 2))/(nCr(i+j+k, 2)+k*R)",
         "C1R": "(nCr(j, 2))/(nCr(i+j+k, 2)+k*R)",
         "C2": "(i*j)/(nCr(i+j+k, 2)+k*R)",
         "C3L": "(i*k)/(nCr(i+j+k, 2)+k*R)",
         "C3R": "(j*k)/(nCr(i+j+k, 2)+k*R)",
         "C4": "(nCr(k, 2))/(nCr(i+j+k, 2)+k*R)"}

probs_mig={"ReA": "(kA*R)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C1LA": "(nCr(iA,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C1RA": "(nCr(jA,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C2A": "(iA*jA)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C3LA": "(iA*kA)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C3RA": "(jA*kA)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C4A": "(nCr(kA,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "ReB": "(kB*R)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C1LB": "(nCr(iB,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C1RB": "(nCr(jB,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C2B": "(iB*jB)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C3LB": "(iB*kB)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C3RB": "(jB*kB)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "C4B": "(nCr(kB,2))/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MitoB": "(iA*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MitoA": "(iB*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MjtoB": "(jA*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MjtoA": "(jB*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MktoB": "(kA*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)",
           "MktoA": "(kB*M)/(nCr(iA+jA+kA,2)+kA*R+nCr(iB+jB+kB,2)+kB*R+iA*M+iB*M+jA*M+jB*M+kA*M+kB*M)"}

def process_eq(eqString, str_to_real):
    while(eqString.find('nCr') != -1):
        ncr_start = eqString.find('nCr(')
        ncr_end = eqString.find(")")
        [f1, f2] = eqString[(ncr_start+len('nCr(')):ncr_end].split(",")
        f1_num = process_eq(f1, str_to_real)
        f2_num = process_eq(f2, str_to_real)
        ncr_out = nCr(f1_num, f2_num)
        eqString = eqString.replace(eqString[(ncr_start):(ncr_end+1)], str(ncr_out))
    eq_vec = eqString.split("+")
    Sum = 0
    for term in eq_vec:
        Prod = 1
        for factor in term.split("*"):
            if factor in str_to_real:
                Prod *= str_to_real[factor]
            elif isInt(factor): Prod *= float(factor)
            else: print "Problem with "+factor
        Sum += Prod
    return Sum
    
#==============================================================================
# Need to make process_eq_str deal with the numbers and also R, m not real
#def process_eq_str(eqString, dic):
#    while(eqString.find('nCr') != -1):
#        ncr_start = eqString.find('nCr(')
#        ncr_end = eqString.find(")")
#        [f1, f2] = eqString[(ncr_start+len('nCr(')):ncr_end].split(",")
#        f1_num = process_eq(f1, dic)
#        f2_num = process_eq(f2, dic)
#        ncr_out = nCr(f1_num, f2_num)
#        eqString = eqString.replace(eqString[(ncr_start):(ncr_end+1)], str(ncr_out))
#    eq_vec = eqString.split("+")
#    Sum = 0
#    for term in eq_vec:
#        Prod = 1
#        for factor in term.split("*"):
#            if factor in str_to_real:
#                Prod *= dic[factor]
#            elif isInt(factor): Prod *= float(factor)
#            else: Prod
#        Sum += Prod
#    return Sum

#==============================================================================

def make_dictionary_str(all_args):
    dictionary = None    
    if len(all_args) == 3:
        dictionary = {'i': all_args[0],
                      'j': all_args[1],
                      'k': all_args[2],
                      'R': 'R'}
    elif len(all_args) == 6:
        dictionary = {'iA': all_args[0],
                      'jA': all_args[1],
                      'kA': all_args[2],
                      'iB': all_args[3],
                      'jB': all_args[4],
                      'kB': all_args[5],
                      'R': 'R',
                      'M': 'M'}
    else: print "Arguments to make_dictionary wrong, returning None"
    return dictionary     

def make_dictionary(all_args):
    dictionary = None    
    if len(all_args) == 4:
        dictionary = {'i': all_args[0],
                      'j': all_args[1],
                      'k': all_args[2],
                      'R': all_args[3]}
    elif len(all_args) == 8:
        dictionary = {'iA': all_args[0],
                      'jA': all_args[1],
                      'kA': all_args[2],
                      'iB': all_args[3],
                      'jB': all_args[4],
                      'kB': all_args[5],
                      'R': all_args[6],
                      'M': all_args[7]}
    else: print "Arguments to make_dictionary wrong, returning None"
    return dictionary     

def process_prob_name(prob_name, dic1):
    if len(dic1) == 4:
        dic2 = probs
    elif len(dic1) == 8:
        dic2 = probs_mig
    else: 
        print "Problem in process prob with input dictionary"
        return None
    eqString = dic2[prob_name]
    #print eqString
    [string_numerator, string_denominator] = eqString.split("/")
    string_numerator = string_numerator[1:(len(string_numerator)-1)]
    #print string_numerator
    string_denominator = string_denominator[1:(len(string_denominator)-1)]
    #print string_denominator
    prob_real = float(process_eq(string_numerator, dic1))/float(process_eq(string_denominator, dic1))
    #print float(process_eq(string_numerator, dic1))
    #print float(process_eq(string_denominator, dic1))
    return prob_real


def df_to_tmat(input_df, params):
    df = input_df.copy()
    for state1, row in df.iterrows():
        dic = make_dictionary(np.vectorize(float)(stringtoVec(state1)+params))
        for state2 in list(df):
            trans = row[state2]
            if trans == '0': 
                df[state2][state1] = float(0)        
            elif trans != '0':
                df[state2][state1] = process_prob_name(trans, dic)
    return df.as_matrix()


n2_tmat = df_to_tmat(n2_df, [1.5])
for i in range(len(n2_tmat)):
    if not np.isclose(sum(n2_tmat[i, :]), 1.0):
        print "Problem with row "+str(i)+", sum is "+str(sum(n2_tmat[i, :]))

n2_mig_tmat = df_to_tmat(n2_mig_df, [1, 0.01])
for i in range(len(n2_tmat)):
    if not np.isclose(sum(n2_tmat[i, :]), 1.0):
        print "Problem with row "+str(i)+", sum is "+str(sum(n2_tmat[i, :]))

## default tolerance is 10^(-5)

#==============================================================================
# Running process
#==============================================================================

#eqht_states = ['002', '110', '{000}{002}', '{002}{000}', '{000}{110}', '{110}{000}']
#uneqht_states = ['101', '011', '{000}{101}', '{101}{000}', '{000}{011}', '{011}{000}']

#==============================================================================
# Tree heights unequal is enters 4,5,6,7,8
# ['110', '011', '210', '120', '110']
#==============================================================================

uneqht_nomig_states = ['101', '011', '210', '120', '110']
uneqht_mig_states = []

for state in list(n2_mig_df):
    if vectoString(removePopulationLabel(state)) in uneqht_nomig_states:
        uneqht_mig_states.append(state)
        
for state in uneqht_mig_states:
    print vectoString(removePopulationLabel(state))

len(uneqht_mig_states)
#uneqht_states = []
#uneqht_states.append(uneqht_nomig_states)
#uneqht_states.append(uneqht_mig_states)

#uneqht_num = []
#for state in uneqht_nomig_states:
#    if state in mig_states:
#        uneqht_num.append(mig_states.index(state))
#    else:
#        print "Problem with state"+state
#        break
#
#
#uneqht_num = np.zeros(len(uneqht_states))
#for i in range(len(uneqht_states)):
#    state = uneqht_states[i]
#    if state in mig_states: uneqht_num[i] = mig_states.index(state)
#    elif state in nomig_states: uneqht_num[i] = nomig_states.index(state)
#    else: 
#        print "Problem with " + state
#        break
#


def run_process_fullpath(tmat):
    if np.shape(tmat)[0] == 9: 
        start = list(n2_df).index('002')
        end = [list(n2_df).index('001')]
    elif np.shape(tmat)[0] == 46: 
        start = list(n2_mig_df).index('{001}{001}')
        end = [list(n2_mig_df).index('{000}{001}'), list(n2_mig_df).index('{001}{000}')] 
    steps = 0
    total_steps = 10000
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
            path[steps] = i
    return(path[0:(steps+1)].astype(int))
  
def run_process(tmat):
    if np.shape(tmat)[0] == 9: 
        start = list(n2_df).index('002')
        end = [list(n2_df).index('001')]
    elif np.shape(tmat)[0] == 46: 
        start = list(n2_mig_df).index('{001}{001}')
        end = [list(n2_mig_df).index('{000}{001}'), list(n2_mig_df).index('{001}{000}')] 
    steps = 0
    total_steps = 10000
    #path = np.zeros(total_steps)
    i = start
    prev = None
    #path[0] = i
    while i not in end and i not in uneqht_mig_states and steps < total_steps:
        row = tmat[i,:]
        if not np.isclose(sum(tmat[i, :]), 1.0):
            print "Not equal to 1 at " + str(i)
            break
        else:
            j = rv_discrete(values=(range(len(row)),row)).rvs(size=1)[0]
            prev = i
            i = j
            steps = steps + 1    
            #path[steps] = i
    if i in uneqht_states:
        return(i)
    elif i in end:
        return(prev)
    else: 
        print "Process ended when in state: "+str(i)
        return None

def equal_heights(num):
    if num in eqht_num: 
        return True
    elif num in uneqht_num:
        return False
    else: 
        print "Problem with index: "+str(num)        
        return None


# 002 equal
# 110 equal
# 101 unequal
# 011 unequal


#eqht_states = ['002', '110', '{000}{002}', '{002}{000}', '{000}{110}', '{110}{000}']
#uneqht_states = ['101', '011', '{000}{101}', '{101}{000}', '{000}{011}', '{011}{000}']
#
#eqht_num = np.zeros(len(eqht_states))
#for i in range(len(eqht_states)):
#    state = eqht_states[i]
#    if state in mig_states: eqht_num[i] = mig_states.index(state)
#    elif state in nomig_states: eqht_num[i] = nomig_states.index(state)
#    else: 
#        print "Problem with " + state
#        break
#
#uneqht_num = np.zeros(len(uneqht_states))
#for i in range(len(uneqht_states)):
#    state = uneqht_states[i]
#    if state in mig_states: uneqht_num[i] = mig_states.index(state)
#    elif state in nomig_states: uneqht_num[i] = nomig_states.index(state)
#    else: 
#        print "Problem with " + state
#        break
#
#eqht_num = eqht_num.astype(int)
#uneqht_num = uneqht_num.astype(int)
#
#eqht_num
#uneqht_num
#
#def equal_heights_fullpath(path, mig_tag):
#    ## m_tag is boolean flag for if migration included or not
#    states = []
#    pathlen = len(path)
#    if mig_tag: 
#        states = mig_states
#        end = [states.index('{001}{000}'), states.index('{000}{001}')]
#    else: 
#        states = nomig_states
#        end = [states.index('001')]
#    if path[len(path)-1] not in end:
#        print "Problem: last state was " + str(path[len(path)-1])
#        return None
#    elif states[path[pathlen-2]] in eqht_states: return True
#    elif states[path[pathlen-2]] in uneqht_states: return False
#    else:
#        print "Problem: penultimate state was "+str(path[len(path)-2])
#        return None
#
#def equal_heights(num):
#    if num in eqht_num: 
#        return True
#    elif num in uneqht_num:
#        return False
#    else: 
#        print "Problem with index: "+str(num)        
#        return None

#==============================================================================
# Run simulations
#==============================================================================

nomig_states = list(n2_df)
mig_states = list(n2_mig_df)

R_vals = [10**(-3), 10**(-2), 10**(-1), 10**(0), 10**(1), 10**(2)]
M_vals = [10**(1), 10**(0), 10**(-1), 10**(-2)]

sims_m1 = np.zeros(len(R_vals))

M = 1
R_index = 4
R = R_vals[R_index]
tmat = df_to_tmat(n2_mig_df, [R, M])

num_sims = 10000
sims = 0
time1 = time.time()
for i in range(num_sims):
    sims += equal_heights(run_process(tmat))
sims_m1[R_index] = float(sims)/num_sims
time.time()-time1
time1 = time.time()


#==============================================================================
### Running simulations with no migration... decreases then increases?? 
#==============================================================================

#R_index = 2
#R = R_vals[R_index]
#tmat = df_to_tmat(n2_df, [R])

#num_sims = 100
#sims = 0
#time1 = time.time()
#for i in range(num_sims):
#    state = run_process(tmat)
#    print "At "+str(i)+": "+nomig_states[state]+" is "+str(equal_heights(state))
#sims_noMig[R_index] = float(sims)/num_sims
#time.time()-time1
#time1 = time.time()
#
#n2_df
#
#num_sims = 10000
#sims = 0
#time1 = time.time()
#for i in range(num_sims):
#    sims += equal_heights(run_process(tmat))
#sims_noMig[R_index] = float(sims)/num_sims
#time.time()-time1
#time1 = time.time()

# sims_noMig[0] = 0.99909999999999999 #  8.853049039840698
# sims_noMig[1] = 0.99480000000000002  # 0.15776090224583944
# sims_noMig[2] = 0.94350000000000001 # 0.2178394357363383
# sims_noMig[3] = 0.78879999999999995 # 0.5808733503023783
# sims_noMig[4] = 0.91779999999999995 # 1.5120146989822387
# sims_noMig[5] = 0.98999999999999999 # 1.8248772025108337


#==============================================================================
### Why does probability of equal tree height decrease then increase? 
#==============================================================================

#n2_df['001']
#tmat[:,list(n2_df).index('001')]
#n2_mig_df[['{001}{000}', '{000}{001}']]
#M = 0.1
#tmat = df_to_tmat(n2_mig_df, [R, M])
#tmat[:, list(n2_mig_df).index('{001}{000}')]


#f = open('pequalht_M1.csv','w')
#f.write('R,M,pequal_ht\n')
#for R_index in range(len(sims_m1)):
#    f.write(str(R_vals[R_index])+str(M)+','+str(sims_m1[R_index])+'\n')
#f.close()

#sims_m1 = np.zeros(len(R_vals))
#
#M = 1
#R_index = 4
#R = R_vals[R_index]
#tmat = df_to_tmat(n2_mig_df, [R, M])
#
#num_sims = 10000
#sims = 0
#time1 = time.time()
#for i in range(num_sims):
#    sims += equal_heights(run_process(tmat))
#sims_m1[R_index] = float(sims)/num_sims
#time.time()-time1
#time1 = time.time()

## Changed process eq

# sims_m1[0] = 0.99729999999999996 # 0.875335419178009
# sims_m1[1] = 0.96830000000000005 # 1.0028655529022217
# sims_m1[2] = 0.77910000000000001 # 78.71453285217285
# sims_m1[3] = 0.55740000000000001 #
# sims_m1[4] = 0.90200000000000002 # 4.4523774822553
# sims_m1[5] = #


#f = open('pequalht_M0001R1.txt','w')
#f.write('M = '+str(M)+', R = '+str(R)+": "+str(sims_m0001[R_index]))
#f.close()


# sims_m0001[0] = 0.99580000000000002 # 0.593396250406901
# sims_m0001[1] = 0.95930000000000004 # 3.502032967408498
# sims_m0001[2] =  0.95930000000000004 # 29.849683133761086


#sims_m001 = np.zeros(len(R_vals))
#sims_m001[0] = 0.99580000000000002 # 0.3521139661471049
#sims_m001[1] = 0.96189999999999998 # 0.5994714657465617
#sims_m001[2] = 0.72250000000000003 # 2.9409512162208555
#sims_m001[3] = 0.49869999999999998 # 15.244563631216685
#sims_m001[4] = 0.90749999999999997 # 30.366686634222667
#sims_m001[5] = 0.9889 #
#
#sims_m01 = np.zeros(len(R_vals))
#sims_m01[0] = 0.996 # 0.5848225673039754
#sims_m01[1] = 0.96579999999999999 # 25.798851013183594
#sims_m01[2] = 0.73919999999999997 # 0.8481764674186707
#sims_m01[3] = 0.51119999999999999 # 2.961593035856883
#sims_m01[4] = 0.90590000000000004 # 4.991067500909169
#sims_m01[5] = 0.9889 # 5.631102331479391
#
#sims_m1 = np.zeros(len(R_vals))
#sims_m1[0] = 0.99770000000000003 # 5.94907337029775 
#sims_m1[1] = 0.97150000000000003 # 0.8901789665222168
#sims_m1[2] = 0.78310000000000002 # 1.1001760999361674
#sims_m1[3] = 0.56059999999999999 # 2.4536336143811543
#sims_m1[4] = 0.90259999999999996 # 4.148743200302124
#sims_m1[5] = 0.99029999999999996 # 263.89752197265625
#

#sims_m10 = np.zeros(len(R_vals))
#sims_m10[0] = 0.99790000000000001 # 5.760323683420817
#sims_m10[1] = 0.9718 # 6.03586608171463
#sims_m10[2] = 0.79690000000000005 # 7.044858415921529
#sims_m10[3] = 0.58899999999999997 # 788.2147619724274
#sims_m10[4] = 0.90129999999999999 # 21.941480016708375
#sims_m10[5] = 0.98980000000000001 # 21.922354833285013
#
#sims_df = pd.DataFrame([sims_m001, sims_m01, sims_m1, sims_m10], columns = list(R_vals), index=list(M_vals))
#
#sims_df.to_csv('sims_df.csv')


#f = open('pequalht_M01.csv','w')
#f.write('R,M,pequal_ht\n')
#for R_index in range(len(R_vals)):
#    f.write(str(R_vals[R_index])+','+str(M)+','+str(sum(sims)/num_sims)+'\n')

#f = open('pequalht_M10.csv','w')
#f.write('R,M,pequal_ht\n')
#for R_index in range(len(M_vals)):
#    for M_index in range(len(M_vals)):
#        tmat = df_to_tmat(n2_mig_df, [R_vals[R_index], M_vals[M_index]])
#        time1 = time.time()
#        num_sims = 10000
#        sims = np.zeros(num_sims)
#        time1 = time.time()
#        for i in range(num_sims):
#            sims[i] = equal_heights(run_process(tmat))
#        f.write(str(R_vals[R_index])+','+str(M_vals[M_index])+','+str(sum(sims)/num_sims)+'\n')
#        print 'Finished with R:'+str(R)+', M:'+str(M_vals[M_index])+' sims'
#        print time.time()-time1
#        time1 = time.time()
#f.close()


## In 'pequalht_M1.csv','w' all M = 1 values

## How long does each one take?

#for R_index in range(len(R_vals)):
#    for M_index in range(len(M_vals)):
#        print "R = "+ str(R_vals[R_index])+", M = "+str(M_vals[M_index])+": "       
#        time1 = time.time()
#        tmat = df_to_tmat(n2_mig_df, [R_vals[R_index], M_vals[M_index]])
#        time1 = time.time()
#        p = equal_heights(run_process(tmat))
#        time_pequal = (time.time()-time1)
#        print (time_pequal*10000)/60

# R = 0.001, M = 10: 
# 1.33848190308
# 356.944402217865/60
# 5.94907337029775

# R = 0.001, M = 1: 
# 1.10650062561

# R = 0.001, M = 0.1: 
# 0.356674194336

# R = 0.001, M = 0.01: 
# 0.572999318441

# R = 0.01, M = 10: 
# 22.7641661962

# R = 0.01, M = 1: 
# 1.74351533254

# R = 0.01, M = 0.1: 
# 0.263333320618

# R = 0.01, M = 0.01: 
# 0.269333521525

# R = 0.1, M = 10: 
# 6.26516342163

# R = 0.1, M = 1: 
# 0.505010286967

# R = 0.1, M = 0.1: 
# 0.528335571289

# R = 0.1, M = 0.01: 
# 3.60933939616

# R = 1, M = 10: 
# 44.7358290354

# R = 1, M = 1: 
# 6.61249955495

# R = 1, M = 0.1: 
# 4.22298908234

# R = 1, M = 0.01: 
# 2.24900245667

# R = 10, M = 10: 
# 13.8501723607

# R = 10, M = 1: 
# 2.33002503713

# R = 10, M = 0.1: 
# 3.37966283162

# R = 10, M = 0.01: 
# 53.3243020376

# R = 100, M = 10: 
# 10.8083486557

# R = 100, M = 1: 
# 6.02682431539

# R = 100, M = 0.1: 
# 0.977675120036

# R = 100, M = 0.01: 
# 10.0778341293

#time1 = time.time()
#tmat = df_to_tmat(n2_mig_df, [R_vals[R_index], M_vals[M_index]])
#time.time()-time1
#num_sims = 10000
#sims = np.zeros(num_sims)
#time1 = time.time()
#for i in range(num_sims):
#    sims[i] = equal_heights(run_process(tmat))
#    #f.write(str(R_vals[R_index])+','+str(M_vals[M_index])+','+str(sum(sims)/num_sims)+'\n')
#    #print 'Finished with R:'+str(R)+', M:'+str(M_vals[M_index])+' sims'
#    #print time.time()-time1
#time1 = time.time()
#



#==============================================================================
# Lumping
#==============================================================================

##

#def combine_rc(df, name1, name2):
#    ## new_col defines the way of combining two cols    
#    new_col = df.loc[:, name1] + df.loc[:, name2]
#    df[name1] = new_col
#    df = df.drop([name2], axis = 1)
#    df = df.rename(index=str, columns={name1: (name1+name2)})
#    ## new_row defines the way of combining two row    
#    new_row = df.loc[name1,:]+df.loc[name2,:]
#    df.loc[name1,:] = new_row
#    df = df.drop([name2])
#    df = df.rename({name1: (name1+name2)})
#    return df
#
#def combine_rows(row1, row2): 
#    return np.vectorize(string_add)(row1.apply(no_chira), row2.apply(no_chira))
#
#def take_one(s1, s2):
#    onething = ''
#    if s1 == '0' and s2 != '0':
#        onething = s2
#    elif s1 != '0' and s2 == '0':
#        onething = s1
#    elif no_chira(s1) == no_chira(s2):
#        onething = s1
#    else: print "Problem with " + s1 + " and " + s2
#    return onething
#
#def combine_cols(col1, col2):
#    return np.vectorize(take_one)(col1.apply(no_chira), col2.apply(no_chira))
#
#def lump(df, name1, name2):
#    ## new_col defines the way of combining two cols
#    #print "Original df: "
#    #print df
#    col1 = df.loc[:, name1]
#    col2 = df.loc[:, name2]        
#    df[name1] = combine_cols(col1, col2)
#    #print "Combined cols: "
#    #print df
#    df = df.drop([name2], axis = 1)
#    df = df.rename(index=str, columns={name1: (name1+name2)})
#    #print "Drop cols, rename: "
#    #print df    
#    ## new_row defines the way of combining two row
#    row1 = df.loc[name1, :]
#    row2 = df.loc[name2, :]   
#    df.loc[name1,:] = combine_rows(row1, row2)
#    #print "Combine rows: "
#    #print df
#    df = df.drop([name2])
#    df = df.rename({name1: (name1+name2)})
#    #print "Drop rows, rename: "
#    #print df
#    return df
#
#def lump_df(df):
#    state_names = list(df)
#    for state in state_names:
#        if vectoString(state) in list(df):
#            if state[0] != state[1]:
#                asym_state = state[1] + state[0] + state[2]
#                if asym_state not in state_names:
#                    print asym_state + " is not in states"
#                df = lump(df, vectoString(state), vectoString(asym_state))
#    return df
#
#n2_df = make_transition_df(find_states(2))
#n2_df_lumped = lump_df(n2_df)
#n2_df_lumped
#make_graphviz(n2_df_lumped, "lumped_n2.dot")
#
#n3_df = make_transition_df(find_states(3))
#n3_df_lumped = lump_df(n3_df)
#make_graphviz(n3_df_lumped, "lumped_n3.dot")
#
#
### Things to consider:
### Does lumping more than 2 things make it different?
### Lumping by populations instead of lineage number
#
#
#def lumpMig(df, name1, name2):
#    ## new_col defines the way of combining two cols
#    #print "Original df: "
#    #print df
#    col1 = df.loc[:, name1]
#    col2 = df.loc[:, name2]        
#    df[name1] = combine_cols(col1, col2)
#    #print "Combined cols: "
#    #print df
#    df = df.drop([name2], axis = 1)
#    df = df.rename(index=str, columns={name1: (name1+name2)})
#    #print "Drop cols, rename: "
#    #print df    
#    ## new_row defines the way of combining two row
#    row1 = df.loc[name1, :]
#    row2 = df.loc[name2, :]   
#    df.loc[name1,:] = combine_rows(row1, row2)
#    #print "Combine rows: "
#    #print df
#    df = df.drop([name2])
#    df = df.rename({name1: (name1+name2)})
#    #print "Drop rows, rename: "
#    #print df
#    return df


### Names are going to get too long
### Ned to deal with that

#short_names = {"s0": vectoString(states[0])}
#
#def lump_dfMig(df):
#    state_names = list(df)
#    for state in state_names:
#        if vectoString(state) in list(df):
#            if state[1] != state[2]:
#                asym_state = state[2] + state[1] + state[3]
#                if asym_state not in state_names:
#                    print asym_state + " is not in states"
#            if state[6] != state[7]:
#                asym_state = state[7] + state[6] + state[8]
#                if asym_state not in state_names:
#                    print asym_state + " is not in states"
#                df = lump(df, vectoString(state), vectoString(asym_state))
#            if state[1] != state[2]:
#                asym_state = state[2] + state[1] + state[3]
#                if asym_state not in state_names:
#                    print asym_state + " is not in states"
#            
#    return df
#
#
#
#
#n = 2
#states = find_states(n)
#state_strings = find_stateNames(states)
#
#np.shape(states)
#mig_states = find_mig_states(n, states)
#
#np.shape(mig_states)
#
#for state in mig_states:
#    print state
#
#t_df = make_migT_df(mig_states)
#t_df.head()
#
#t_df = t_df.applymap(strip)
#t_df.head()
#
#for idx in t_df.index:
#    print idx
#
#make_graphviz_mig(t_df, "n3_ReMig.dot")
