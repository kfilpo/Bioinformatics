# -*- coding: utf-8 -*-
"""
Created on Mon Apr 23 19:58:42 2018

@author: Kerri Norton
"""

from matplotlib import pyplot

from random import randint, sample

from numpy import array, cov, diag, dot, linalg, ones
from numpy import outer, random, sqrt, vstack, zeros


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

  testDataA = random.random((1000,2)) # No clumps

  centers, clusters = kMeans(testDataA, 3)


  testDataB1 = random.normal(0.0, 2.0, (100,2))
  testDataB2 = random.normal(7.0, 2.0, (100,2))
  testDataB = vstack([testDataB1, testDataB2]) # Two clumps

  centers, clusters = kMeans(testDataB, 2)


  colors = ['#FF0000','#00FF00','#0000FF',
            '#FFFF00','#00FFFF','#FF00FF']

  for i, cluster in enumerate(clusters):
     x, y = zip(*cluster)
     color = colors[i % len(colors)]
     pyplot.scatter(x, y, c=color, marker='o')

  x, y = zip(*centers)
  pyplot.scatter(x, y, s=40, c='black', marker='o')
  pyplot.show()
