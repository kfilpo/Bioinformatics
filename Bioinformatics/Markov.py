# -*- coding: utf-8 -*-
"""
Created on Tue Apr  3 13:06:33 2018

@author: Kerri-Ann Norton
"""

from matplotlib import pyplot
from numpy import array, empty, identity, dot, ones, zeros, log
from scipy.misc import comb 
from scipy.stats import binom, geom, poisson


def getNextGenPop (currentPop, randVar):
  
  progeny = randVar.rvs(size=currentPop)
  print("Progeny: ", progeny)
  nextPop = progeny.sum()
  print("NewPop: ", nextPop)
  return nextPop


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
  for i in range(50):
    pop = getNextGenPop(pop, poissRandVar)
    print("Generation:%3d Population:%7d" % (i, pop))