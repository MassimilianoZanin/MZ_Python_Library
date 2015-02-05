# -*- coding: utf-8 -*-
"""
---------------------------------------------
Library for Time Series analysis
---------------------------------------------

Log of changes:

    2015.01.30 Added function for Granger Causality

"""


from statsmodels.tsa.stattools import grangercausalitytests
import numpy



def GrangerTest( X, Y, maxlag ):

    myVerbose = False
    
    fullTS = [Y, X]
    fullTS = numpy.asarray( fullTS ).T.tolist()
    
    try:
        res = grangercausalitytests( fullTS, maxlag, True, myVerbose )
    except:
        return (1.0, 0)
        pass
    
    minFValue = 1.0
    bestLag = 0
    for n in res.keys():
        temp = res[n][0]
        temp = temp['params_ftest']
        temp = temp[0]
        
        if temp < minFValue:
            minFValue = temp
            bestLag = n

    return (minFValue, bestLag)
        