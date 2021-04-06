# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 14:53:23 2018

@author: Kerri-Ann Norton
Markov and Hidden Markov Model
"""

from matplotlib import pyplot
from numpy import array, empty, identity, dot, ones, zeros, log
from scipy.special import comb 
from scipy.stats import binom, geom, poisson


def binomialProbability(n, k, p):

  return comb(n, k) * p**k * (1-p) ** (n-k)


def getNextGenPop (currentPop, randVar):
  
  progeny = randVar.rvs(size=currentPop)
  nextPop = progeny.sum()
  return nextPop


def viterbi(obs, pStart, pTrans, pEmit):

  nStates = len(pStart)
  states = range(nStates)

  scores = empty(nStates)
  scoresPrev = pStart + pEmit[:,obs[0]]
  pathsPrev = dict([(i, [i]) for i in states])

  for val in obs[1:]:
    paths = {}

    for i in states:

      options = scoresPrev + pTrans[:,i] + pEmit[i,val]
      bestState = options.argmax()
      scores[i] = options.max()   
      paths[i] = pathsPrev[bestState] + [i]

    pathsPrev = paths
    scoresPrev = array(scores)

  endState = scores.argmax()
  logProb = scores.max()

  return logProb, paths[endState]


def forwardBackward(obs, pStart, pTrans, pEmit):

  n = len(obs)
  nStates = len(pStart)
  I = identity(nStates)

  fwd = empty([n+1,nStates])
  fwd[0] = pStart

  for i, val in enumerate(obs):
    fProb = dot(pEmit[:,val]*I, dot(pTrans, fwd[i])) 
    fwd[i+1] = fProb / fProb.sum()

  bwd = ones(nStates)
  smooth = empty([n+1,nStates])
  smooth[-1] = fwd[-1]

  for i in range(n-1, -1, -1):
    bwd = dot(pTrans, dot(pEmit[:,obs[i]]*I, bwd))
    bwd /= bwd.sum()
    prob = fwd[i] * bwd
    smooth[i] = prob / prob.sum()

  return smooth



if __name__ == '__main__':
    
     # Simple population Markov Chain

  p = 0.0025
  geomRandomVar = geom(p)
  lengths = array(range(1, 1000))
  probs = geomRandomVar.pmf(lengths)

  pyplot.plot(lengths, probs)
  pyplot.show()

  from scipy.stats import poisson
  rate = 1.02
  poissRandVar = poisson(rate)

  pop = 25
  for i in range(100):
    pop = getNextGenPop(pop, poissRandVar)
    print("Generation:%3d Population:%7d" % (i, pop))




  # HMM test example - protein burried/exposed states

  expTypes = ['-','*',]
  aaTypes = ['A','C','D','E','F','G','H','I','K','L',
             'M','N','P','Q','R','S','T','V','W','Y']

  # Initialisation

  nExp = len(expTypes)		# Number of exposure categories
  nAmino = len(aaTypes)		# Number of amino acid types
  nStates = nExp * nAmino	# Number of HMM states

  pStart = zeros(nStates, float)		 # Starting probabilities
  pTrans = zeros((nStates, nStates) , float) # Transition probabilities
  pEmit  = zeros((nStates, nAmino) , float)	 # Emission probabilities

  indexDict = {}
  stateDict = {}
  index = 0
  for exposure in expTypes:
    for aminoAcid in aaTypes:
      stateKey = (exposure, aminoAcid)
      indexDict[stateKey] = index
      stateDict[index] = stateKey
      index += 1
  
  # Estimate transition probabilities from database counts
  
  fileName = 'PdbSeqExposureCategories.txt'
  fileObj = open(fileName,'r')

  line1 = fileObj.readline()
  line2 = fileObj.readline()

  while line1 and line2:
    sequence = line1.strip()
    exposure = line2.strip()

    n = len(sequence)  
    for i in range(n-2):
      aa1, aa2 = sequence[i:i+2]
      exp1, exp2 = exposure[i:i+2]
      stateKey1 = (exp1, aa1)
      stateKey2 = (exp2, aa2)
      index1 = indexDict.get(stateKey1)
      index2 = indexDict.get(stateKey2)

      if index1 is None or index2 is None:
        continue

      pStart[index1] += 1.0
      pTrans[index1, index2] += 1.0

    pStart[index2] += 1
    line1 = fileObj.readline()
    line2 = fileObj.readline()

  pStart /= pStart.sum()
  for i in range(nStates):
    pTrans[i] /= pTrans[i].sum()
  
  # Setup trivial emission probabilities
  
  for exposure in expTypes:
    for aminoIndex, aminoAcid in enumerate(aaTypes):
      stateIndex = indexDict[(exposure, aminoAcid)]
      pEmit[stateIndex, aminoIndex] = 1.0
  
  # HMM test data
  
  seq = "MYGKIIFVLLLSEIVSISASSTTGVAMHTSTSSSVTKSYISSQTNDTHKRDTYAATPRAH"\
        "EVSEISVRTVYPPEEETGERVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK"\
        "SPSDVKPLPSPDTDVPLSSVEIENPETSDQ"
        
  obs = [aaTypes.index(aa) for aa in seq]

  adj = 1e-99
  logStart = log(pStart+adj)
  logTrans = log(pTrans+adj)
  logEmit  = log(pEmit+adj)
  
  # Best route - Viterbi
  
  logProbScore, path = viterbi(obs, logStart, logTrans, logEmit)                            

  bestExpCodes = ''.join([stateDict[i][0] for i in path])

  print(seq)
  print(bestExpCodes)
  
  # Positional state probabilitied - Forward-backward
  
  smooth = forwardBackward(obs, pStart, pTrans, pEmit)

  buriedList = []
  exposeList = []

  for values in smooth:
    exposeList.append(sum(values[:20]))
    buriedList.append(sum(values[20:]))

  xAxisValues = list(range(len(exposeList))) # Sequence positions 

  pyplot.plot(xAxisValues, exposeList, c='#A0A0A0')
  pyplot.plot(xAxisValues, buriedList, c='#000000')
  pyplot.show()
