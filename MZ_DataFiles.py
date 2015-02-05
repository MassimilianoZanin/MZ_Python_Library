# -*- coding: utf-8 -*-
"""
---------------------------------------------
Library for Data Files handling
---------------------------------------------

Log of changes:

    2015.01.30 Added function for loading CSV files to float

"""


import csv


def ReadCSVToFloat( fileName, myDelimiter ):

    AllData = []
    
    with open(fileName, 'rb') as csvfile:
        
        spamreader = csv.reader( csvfile, delimiter = myDelimiter, quotechar = '|' )
        
        for row in spamreader:
            tempRow = []
            
            for item in row:
                tempRow.append( float(item) )
                
            AllData.append(tempRow)
            
    return AllData
    
    
    
    