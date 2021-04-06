#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 21:20:41 2020

@author: katherinefilpolopez
"""
from matplotlib import pyplot

from random import randint, sample

from numpy import array, cov, diag, dot, linalg, ones
from numpy import outer, random, sqrt, vstack, zeros
import lab5 as MUT
import MultipleSeqAlignment as MS

listData = []

seq1 = 'GCACGTATTGATTGGCCTGTACCTA'
listData.append(seq1)
seq11 = 'GCACGTATTGATTGGCCTGTACCTA'
listData.append(seq11)
seq21 = 'GCACGTATTGATTGGCCTGTACCTA'
listData.append(seq21)
test = MUT.Mutation(seq1)
test2 = MUT.Mutation(seq11)
test3 = MUT.Mutation(seq21)

seq2 = test.frameshiftDeletion()
listData.append(seq2)
seq3 = test.frameshiftInsertion()
listData.append(seq3)
seq4 = test.pointMutation()
listData.append(seq4)
seq5 = test.frameshiftInsertion()
listData.append(seq5)
seq6 = test.pointMutation()
listData.append(seq6)
seq7 = test.frameshiftDeletion()
listData.append(seq7)
seq8 = test.frameshiftInsertion()
listData.append(seq8)
seq9 = test.pointMutation()
listData.append(seq9)
seq10 = test.frameshiftDeletion()
listData.append(seq10)

seq12 = test.frameshiftDeletion()
listData.append(seq12)
seq13 = test.frameshiftInsertion()
listData.append(seq13)
seq14 = test.pointMutation()
listData.append(seq14)
seq15 = test.frameshiftInsertion()
listData.append(seq15)
seq16 = test.pointMutation()
listData.append(seq16)
seq17 = test.frameshiftDeletion()
listData.append(seq17)
seq18 = test.frameshiftInsertion()
listData.append(seq18)
seq19 = test.pointMutation()
listData.append(seq19)
seq20 = test.frameshiftDeletion()
listData.append(seq20)

seq22 = test.frameshiftDeletion()
seq23 = test.frameshiftInsertion()
seq24 = test.pointMutation()
seq25 = test.frameshiftInsertion()
seq26 = test.pointMutation()
seq27 = test.frameshiftDeletion()
seq28 = test.frameshiftInsertion()
seq29 = test.pointMutation()
seq30 = test.frameshiftDeletion()


#MS.consensusMultipleAlign(listData, 0.25, MS.BLOSUM62))
testData = []

i = 0
for i in range(len(listData)-1):
    testData.append(MS.pairAlignScore(listData[i], listData[i+1], MS.DNA_2))
    i = i +1
    
print(testData)


def euclideanDist(vectorA, vectorB):
  
  diff = vectorA-vectorB
  
  return sqrt(dot(diff,diff))


def findNeighbours(data, distFunc, threshold):
  
  neighbourDict = {}
  
  n = len(data)
  for i in range(n):
    neighbourDict[i] = []

  for i in range(0,n-1):
    for j in range(i+1,n):
      dist = distFunc(data[i], data[j])
      
      if dist < threshold:
        neighbourDict[i].append(j)
        neighbourDict[j].append(i)

  return neighbourDict

def kMeans(data, k, centers=None):
  
  if centers is None:
    centers = array( sample(list(data), k) )  # list() not needed in Python 2

  change = 1.0
  prev = []

  while change > 1e-8:

    clusters = [[] for x in range(k)]
    for vector in data:
      diffs = centers - vector
      dists = (diffs * diffs).sum(axis=1)
      closest = dists.argmin()
      clusters[closest].append(vector)
     
    change = 0
    for i, cluster in enumerate(clusters):
      cluster = array(cluster)
      center = cluster.sum(axis=0)/len(cluster)
      diff = center - centers[i]
      change += (diff * diff).sum()
      centers[i] = center
    
  return centers, clusters


def kMeansSpread(data, k):

  n = len(data)
  index = randint(0, n-1)
  indices = set([index])
  
  influence = zeros(n)
  while len(indices) < k:
    diff = data - data[index]
    sumSq = (diff * diff).sum(axis=1) + 1.0
    influence += 1.0 / sumSq
    index = influence.argmin()
    
    while index in indices:
      index = randint(0, n-1)
    
    indices.add(index)    
  
  centers = vstack([data[i] for i in indices])
    
  return kMeans(data, k, centers)



if __name__ == '__main__':
    
  print("\nK-means clustering\n")

#  testDataA = random.random((1000,2)) # No clumps

  centers, clusters = kMeans(testData, 3)


  #testDataB1 = random.normal(0.0, 2.0, (100,2))
  #testDataB2 = random.normal(7.0, 2.0, (100,2))
  #testDataB = vstack([testDataB1, testDataB2]) # Two clumps

  centers, clusters = kMeans(testData, 2)


  colors = ['#FF0000','#00FF00','#0000FF',
            '#FFFF00','#00FFFF','#FF00FF']

  for i, cluster in enumerate(clusters):
     x, y = zip(*cluster)
     color = colors[i % len(colors)]
     pyplot.scatter(x, y, c=color, marker='o')

  x, y = zip(*centers)
  pyplot.scatter(x, y, s=40, c='black', marker='o')
  pyplot.show()

