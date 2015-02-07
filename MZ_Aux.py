# -*- coding: utf-8 -*-
"""
---------------------------------------------
Library for auxiliary functions
---------------------------------------------

Log of changes:

    2015.02.06      Functions for motifs analysis

"""


import cPickle as pickle



"""  ------------------------------------------------------------
//  Fast calculation of the cubic root.
"""

def CubeRoot( X ):

    if X < 0.0:
        return -1.0 * pow(-X, 1.0 / 3.0)
        
    return pow(X, 1.0 / 3.0)



def setBit(int_type, offset):
    
    mask = 1 << offset
    return(int_type | mask)


def clearBit(int_type, offset):
    
    mask = ~(1 << offset)
    return(int_type & mask)



def save_object(obj, filename):
    with open(filename, 'wb') as output:
        pickle.dump(obj, output, pickle.HIGHEST_PROTOCOL)

def load_object(filename):
    with open(filename, 'rb') as input:
        obj = pickle.load(input)
    return obj

