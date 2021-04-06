#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 20:54:07 2020

@author: katherinefilpolopez
"""

import re
from math import log, exp
import MultipleSeqAlignment as MS

seqA1 = 'RFdataseq1.txt'
tfile = open(seqA1, 'r')
myText = tfile.read()

seqA2 ='RFdataseq2.txt'
tfile2 = open(seqA2, 'r')
myText2 = tfile2.read()

seqA3 ='RFdataseq3.txt'
tfile3 = open(seqA3, 'r')
myText3 = tfile3.read()

seqA4 = 'RFdataseq4.txt'
tfile4 = open(seqA4, 'r')
myText4 = tfile4.read()

seqA4b = 'RFdataseq5V2.txt'
tfile4b = open(seqA4b, 'r')
myText4b = tfile4b.read()

seqA5 ='RFdataseq5.txt'
tfile5 = open(seqA5, 'r')
myText5 = tfile5.read()

seqA5b ='RFdataseq5V2.txt'
tfile5b = open(seqA5b, 'r')
myText5b = tfile5b.read()

seqA6 ='RFdataseq6.txt'
tfile6 = open(seqA6, 'r')
myText6 = tfile6.read()

seqA7 ='RFdataseq7.txt'
tfile7 = open(seqA7, 'r')
myText7 = tfile7.read()

#print(MS.pairAlignScore(myText5b, myText, MS.DNA_2))
#print(MS.pairAlignScore(myText5b, myText2, MS.DNA_2))
#print(MS.pairAlignScore(myText5b, myText3, MS.DNA_2))
#print(MS.pairAlignScore(myText5b, myText4, MS.DNA_2))
#
#print(MS.pairAlignScore(myText5b, myText6, MS.DNA_2))
#print(MS.pairAlignScore(myText5b, myText7, MS.DNA_2))


def getDistanceMatrix(seqs, simMatrix):

  n = len(seqs)
  #NbyN
  matrix = [[0.0] * n for x in range(n)]
  
  maxScores = [calcSeqSimilarity(x, x, simMatrix) for x in seqs]

  for i in range(n-1):
    seqA = seqs[i]
  
    for j in range(i+1,n):
      seqB = seqs[j]
      
      score, alignA, alignB = sequenceAlign(seqA, seqB, simMatrix)
      maxScore = max(maxScores[i],maxScores[j])
      dist = maxScore - score
      
      matrix[i][j] = dist
      matrix[j][i] = dist

  return matrix

#ourdistMatrix = getDistanceMatrix(myText, BLOSUM62)
#print(ourdistMatrix)

#Tree construction
def getDistToJunction(distMatrix, i, j):
  
  n = len(distMatrix)
  row = distMatrix[i]
  column = distMatrix[j]
 
  dist = distMatrix[i][j] + (sum(row)-sum(column))/(n-2)
  dist *= 0.5

  return dist

#get joinPair
def getJoinPair(distMatrix):
    n = len(distMatrix) 
    minQ = None
    joinPair = None
    
    for i in range(n-1):
        # 
        sumRow = sum(distMatrix[i]) 
        
        #
        for j  in range(i+1,n):
            #
            sumNext = sum(distMatrix[j])
            
            #
            dist = distMatrix[i][j]
            q = (n-2)*dist -sumRow - sumNext
            #print(i,j,q)
            
            if (minQ is None) or (q < minQ):
                minQ = q
                joinPair = [i,j]

    return joinPair
    
#print(getJoinPair(ourdistMatrix))

def neighbourJoinTree(distMatrix):
 
  joinOrder = []
  n = len(distMatrix)
  tree = list(range(n))  
  
  while n > 2:

    x, y = getJoinPair(distMatrix)

    node = (tree[x], tree[y])
    joinOrder.append(node)
    tree.append(node)

    del tree[y]
    del tree[x]

    distX = getDistToJunction(distMatrix, x, y)
    distY = getDistToJunction(distMatrix, y, x)
  
    distMatrix.append([0] * (n+1))

    for i in range(n):
      if i not in (x,y):

        dist = (distMatrix[x][i]-distX) + (distMatrix[y][i]-distY)
        dist *= 0.5
  
        distMatrix[i].append(dist)
        distMatrix[n][i] = dist

    del distMatrix[y]
    del distMatrix[x]
  
    for row in distMatrix:
      del row[y]
      del row[x]

    n -= 1

  tree = tuple(tree)
  joinOrder.append(tree)
  
  return tree, joinOrder

distMatrix = [[0, -2157, -1784, -1918, -1667,-1627, -2138],
               [-2157, 0, -1844, -1868, -1854, -1768, -2072],
               [-1784, -1844, 0, -1638, -1814, -1551, -1977],
               [-1918, -1868, -1638, 0, -1625, -1549, -1865],
               [-1667, -1854, -1814, -1625, 0, -1470, -1848],
               [-1627, -1768, -1551, -1549, -1470, 0, -1728],
               [-2138, -2072, -1977, -1865, -1848, -1728, 0]]

distMatrix2 = [[0, -2157, -1784, -1558, -1558,-1627, -2138],
               [-2157, 0, -1844, -1590, -1590, -1768, -2072],
               [-1784, -1844, 0, -1341, -1341, -1551, -1977],
               [-1918, -1868, -1638, 0, -1394, -1549, -1865],
               [-1667, -1854, -1814, -1369, 0, -1470, -1848],
               [-1627, -1768, -1551, -1304, -1304, 0, -1728],
               [-2138, -2072, -1977, -1560, -1560, -1728, 0]]

tree, treeJoinOrder = neighbourJoinTree(distMatrix)
print(tree)

tree2, treeJoinOrder2 = neighbourJoinTree(distMatrix2)
print(tree2)
       
seqs = [myText,
        myText2,
        myText3,
        myText4, 
        myText5,
        myText6,
        myText7]


#print('\nConsensus paired alignment')
 
#MS.consensusMultipleAlign(seqs, 0.25, MS.BLOSUM62)