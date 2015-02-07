# -*- coding: utf-8 -*-
"""
---------------------------------------------
Library for Complex Networks analysis
---------------------------------------------

Log of changes:

    2015.02.06      Fast motifs analysis

"""


import MZ_Aux
import numpy




"""  ------------------------------------------------------------
//  Pre-identifies the nodes that are not connected.
//  Used to improve calculation performance on very sparse networks.
"""

def GetConnectedness( AM, numNodes ):
    
    IsConnected = [False] * numNodes
    
    for n1 in xrange(0, numNodes):
        for n2 in xrange(0, numNodes):
            if AM[n1][n2] > 0 or AM[n2][n1] > 0:
                IsConnected[n1] = True
                break
            
    return IsConnected



def _DecryptMotifs( MotifsArray ):
    
    Motifs = [0] * 13
    
    
    # Motif 1
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Motifs[ 0 ] = Motifs[ 0 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5)
    Motifs[ 0 ] = Motifs[ 0 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 0 ] = Motifs[ 0 ] + MotifsArray[ Offset ]


    # Motif 2
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 6)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 5)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 7)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 7)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 1 ] = Motifs[ 1 ] + MotifsArray[ Offset ]


    # Motif 3
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 6)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 2 ] = Motifs[ 2 ] + MotifsArray[ Offset ]


    # Motif 4
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5)
    Motifs[ 3 ] = Motifs[ 3 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 7)
    Motifs[ 3 ] = Motifs[ 3 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6)
    Motifs[ 3 ] = Motifs[ 3 ] + MotifsArray[ Offset ]


    # Motif 5
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 7)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 4 ] = Motifs[ 4 ] + MotifsArray[ Offset ]


    # Motif 6
    
    Offset = MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 5 ] = Motifs[ 5 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 5 ] = Motifs[ 5 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5)
    Motifs[ 5 ] = Motifs[ 5 ] + MotifsArray[ Offset ]


    # Motif 7
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 7)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 6 ] = Motifs[ 6 ] + MotifsArray[ Offset ]


    # Motif 8
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6)
    Motifs[ 7 ] = Motifs[ 7 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 7 ] = Motifs[ 7 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 7 ] = Motifs[ 7 ] + MotifsArray[ Offset ]


    # Motif 9
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 8 ] = Motifs[ 8 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 7)
    Motifs[ 8 ] = Motifs[ 8 ] + MotifsArray[ Offset ]


    # Motif 10
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 7)
    Motifs[ 9 ] = Motifs[ 9 ] + MotifsArray[ Offset ]


    # Motif 11
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 10 ] = Motifs[ 10 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 10 ] = Motifs[ 10 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 10 ] = Motifs[ 10 ] + MotifsArray[ Offset ]


    # Motif 12
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Offset += MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Offset += MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 3)
    Offset += MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Offset += MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 7)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Offset += MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]
    
    Offset = MZ_Aux.setBit(0, 2) +  MZ_Aux.setBit(0, 3)
    Offset += MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 11 ] = Motifs[ 11 ] + MotifsArray[ Offset ]


    # Motif 13
    
    Offset = MZ_Aux.setBit(0, 1) +  MZ_Aux.setBit(0, 2)
    Offset += MZ_Aux.setBit(0, 3) +  MZ_Aux.setBit(0, 5) +  MZ_Aux.setBit(0, 6) +  MZ_Aux.setBit(0, 7)
    Motifs[ 12 ] = Motifs[ 12 ] + MotifsArray[ Offset ]

    return Motifs





def GetMotifs( AM, numNodes ):
    
    MotifsArray = numpy.empty( (512, 1) )
    MotifsArray.fill( 0 )

    IsConnected = GetConnectedness( AM, numNodes )
    
    for node1 in xrange(0, numNodes):
        
        if IsConnected[node1] == False:
            continue
        
        for node2 in xrange(node1 + 1, numNodes):
            if IsConnected[node2] == False:
                continue
            
            XReg = 0
            
            if AM[node1][node2] > 0:
                XReg = MZ_Aux.setBit(XReg, 1)
            if AM[node2][node1] > 0:
                XReg = MZ_Aux.setBit(XReg, 3)

            for node3 in xrange(node2 + 1, numNodes):
                
                if IsConnected[node3] == False:
                    continue
                
                Reg = XReg
                if AM[node1][node3] > 0:
                    XReg = MZ_Aux.setBit(XReg, 2)
                if AM[node2][node3] > 0:
                    XReg = MZ_Aux.setBit(XReg, 5)
                if AM[node3][node1] > 0:
                    XReg = MZ_Aux.setBit(XReg, 6)
                if AM[node3][node2] > 0:
                    XReg = MZ_Aux.setBit(XReg, 7)
                
                MotifsArray[ Reg ] = MotifsArray[ Reg ] + 1

    Motifs = _DecryptMotifs( MotifsArray )
    return Motifs




