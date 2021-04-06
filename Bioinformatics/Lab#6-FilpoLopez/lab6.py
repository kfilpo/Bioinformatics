# -*- coding: utf-8 -*-
"""
NeighborTree
Created on Sat Feb 10 14:47:42 2018

@author: Kerri Norton
"""

from math import log, exp

#from MultipleSeqAlignment import profile, profileAlign
from Dictionaries import STANDARD_GENETIC_CODE as SGC, DNA_1
from sequenceAlignment import DNA_2, BLOSUM62, sequenceAlign, calcSeqSimilarity


#calculate distance matrix 
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

#ourdistMatrix = getDistanceMatrix(seqs, BLOSUM62)
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


def getJoinPairforUPGMA(element, distMatrix):
    n = len(distMatrix) 
    minQ = None
    joinPair = None
    
    for i in range(n-1):
        for j  in range(i+1,n):
        
            dist = distMatrix[i][j]
            dist2 = distMatrix[j][i]
            q = (element *dist) + (element *dist2)
            print(i,j,q)
            
            if (minQ is None) or (q < minQ):
                minQ = q
                joinPair = [i,j]

    return joinPair

def UPGMATree(distMatrix):
 
  joinOrder = []
  n = len(distMatrix)
  tree = list(range(n))  
  
  while n > 2:

    x, y = getJoinPairforUPGMA(1, distMatrix)

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


distMatrix2 = [[0, 20, 60, 100, 90],
               [20, 0, 50, 90, 80],
               [60, 50, 0, 40, 50],
               [100, 90, 40, 0, 30],
               [90, 80, 50, 30, 0]]

tree, treeJoinOrder = neighbourJoinTree(distMatrix2)
#print(getJoinPair(distMatrix2))

#print(tree) # Result : (((7, (0, 1)), (4, 5)), ((2, 3), (6, 8)))

distMatrix2test = [[0, 20, 60, 100, 90],
               [20, 0, 50, 90, 80],
               [60, 50, 0, 40, 50],
               [100, 90, 40, 0, 30],
               [90, 80, 50, 30, 0]]

tree2, treeJoinOrder2 = UPGMATree(distMatrix2test)

print(tree2)