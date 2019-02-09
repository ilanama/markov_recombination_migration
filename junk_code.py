# -*- coding: utf-8 -*-
"""
Created on Thu Sep 13 14:52:11 2018

@author: ilanaarbisser
"""


#############################################
##### One way to make transition matrix #####
#############################################

for s1 in range(num_states):
    i1 = states[s1, 0]
    j1 = states[s1, 1]
    k1 = states[s1, 2]
    for s2 in range(num_states):
        i2 = states[s2, 0]
        j2 = states[s2, 1]
        k2 = states[s2, 2]
        transition_matrix[s1, s2] = move(i1, j1, k1, i2, j2 ,k2)

#############################################
### Another way to make transition matrix ###
#############################################

## Let's just start with explicit transitions
## I'm getting rid of c3_l/c3_r because they're the same move, but one requires 2 left/right and one requires 1 left/right and 1 double
## For n = 2, these are mutual exclusive
## For n > 2, these are not mutually exclusive. So I'll just determine the probability of either the 
transitions = {"Re": [1, 1, -1], "C1_ell": [-1, 0, 0], "C1_r": [0, -1, 0], "C2": [-1, -1, 1], "C3_ell": [-1, 0, 0], "C3_r": [0, -1, 0], "C4": [0, 0, -1]}

transition_matrix = np.reshape(["0"]*(num_states ** 2), (num_states, num_states))
transition_matrix[0, 0] = "1"
transition_matrix
transition_df = pd.DataFrame(transition_matrix, index = state_names, columns = state_names)

for s1 in states:    
    for t in transitions:
        move = transitions[t]
        s2 = np.add(s1, move)
        if isinbounds(s2):
            transition_df[vecstateString(s1)][vecstateString(s2)] = t 
transition_df = transition_df.T

transition_df

#######################################
########## More crap ##################
#######################################


## Types of transitions:
## 1 double --> 1 left, 1 right: i++, j++, k--; c
## 1 left, 1 right --> 1 double: i--, j--, k++; re
## 1 left, 1 double --> 1 double: i--, j, k; c
## 1 right, 1 double --> 1 double: i, j--, k; c
## 2 double --> 1 double: i, j, k--; c
## 2 left --> 1 left: i--, j, k; c
## 2 left --> 1 left: i, j--, k; c
## 1 double --> 1 double

## if i2, j2, k2 is in bounds
## i--, j--, k++; c

## i++, j++, k--; re
## i, j, k--; c

## for now this is enough info
## left can recombine with other left or double, which have two diff probabilities
## i--, j, k; c 
## i, j--, k; c

def move(i1, j1, k1, i2, j2, k2):
    move = "0"
    if k1 == 1 and i1 == 0 and j1 == 0 and k2 == 1 and i2 == 0 and j2 == 0:
        move = "1"
    if k2 == k1:
        if i2 == (i1 - 1) and j2 == j1:
            if isinbounds(i2, j2, k2):
                move = "c"
        elif i2 == (i1 - 1) and j2 == j1:
            if isinbounds(i2, j2, k2):            
                move = "c"
    elif k2 == (k1 - 1):
        if i2 == i1 and j2 == j1:
            if isinbounds(i2, j2, k2):            
                move = "c"
        elif i2 == (i1 + 1) and j2 == (j1 + 1):
            if isinbounds(i2, j2, k2):            
                move = "re"
    elif k2 == (k1 + 1):
        if i2 == (i2 - 1) and j2 == (j2 - 1):
            if isinbounds(i2, j2, k2):            
                move = "c"
    return move
