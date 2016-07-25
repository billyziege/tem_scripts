import sys
import os
import argparse
import numpy as np
from scipy import stats
from coordinate_vector_3d import Cartesian3DVector, Cylindrical3DVector, Spherical3DVector

class MyDist(stats.rv_continuous):
  """
  A generic class to provide common features to my distributions.
  """

  def __init__(self,mean=0,scale=1,upper=None,lower=None):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    self.mean = mean
    self.scale = scale
    self.upper = upper
    self.lower = lower
    stats.rv_continuous.__init__(self,momtype=0,a=self.lower,b=self.upper,name="mydist")

  def getSample(self,n):
    """
    Get a sample of n values returned as a list.
    """
    sample = self.rvs(size=n)
    sample = [s + self.mean for s in sample] #Moves distribution
    sample = [s * self.scale for s in sample] #Scales
    return sample

class GaussianGen(MyDist):
  """
  A class to provide convenient access to gaussian distributions.
  """

  def _pdf(self, x):
    """The probability density function."""
    return np.exp(-(x)**2 / 2. ) / np.sqrt(2.0 * np.pi)


class CircularUniformGen(MyDist):
  """
  A class to provide convenient access to circular uniform distribution.
  """

  def __init__(self,mean=0,scale=1,upper=None):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    MyDist(self,mean,scale,upper,0)
     
  def _pdf(self, x):
    """The probability density function."""
    return 1./np.pi

class CircularElipticalGen(MyDist):
  """
  A class to provide convenient access to circular eliptical distribution.
  """

  def __init__(self,mean=0,scale=1,upper=None):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    MyDist(self,mean,scale,upper,0)
     
  def _pdf(self, x):
    """The probability density function."""
    return 3./(2. * np.pi)*np.sqrt(1.-x**2)


class PolarGaussianGen():

  def __init__(self,scale=1):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    self.r_dist = MyDist(0,scale,1,0)
   
  def getSample(self,n):
    """
    Get a sample of n duples returned as a list of duples.
    Using the Box-Mueller theorem.
    """
    x_and_y_sample = []
    for i in range(n):
      x_and_y_sample.append(box_mueller(self.r_dist.scale))
    return x_and_y_sample

  def generateCartesian3DVector(self,n):
    """
    Returns a list of n Cartesian3DVectors with x and y coordinates distributed 
    according to the polar distribution and z set to 0.
    """
    output = []
    x_and_y_coords = self.getSample(n)
    for x_and_y_coord in x_and_y_coords:
      r = Cartesian3DVector(x_and_y_coord[0],x_and_y_coord[1],0)
      output.append(r)
    return output
    
class AnnularGaussianGen(PolarGaussianGen):

  def __init__(self,scale=1,upper=1,lower=0):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    self.r_dist = MyDist(0,scale,upper,lower)
   
  def getSample(self,n):
    """
    Get a sample of n duples returned as a list of duples.
    Using the Box-Mueller theorem.
    """
    x_and_y_sample = []
    for i in range(n):
      x_and_y_sample.append(box_mueller_ring(self.r_dist.scale,self.r_dist.lower,self.r_dist.upper))
    return x_and_y_sample

class PolarUniformGen():

  def __init__(self,scale=1):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    self.r_dist = CircularUniformGen(0,scale,1)
   
  def getSample(self,n):
    """
    Get a sample of n duples returned as a list of duples.
    """
    r_sample = self.r_dist.getSample(n)
    phi_sample = stats.uniform.rvs(size=n)
    phi_sample = [2 * np.pi * s for s in phi_sample] #Scales
    return zip(r_sample,phi_sample)

  def generateCartesian3DVector(self,n):
    """
    Returns a list of n Cartesian3DVectors with x and y coordinates distributed 
    according to the polar distribution and z set to 0.
    """
    output = []
    polar_coords = self.getSample(n)
    for polar_coord in polar_coords:
      r = Cylindrical3DVector(polar_coord[0],polar_coord[1],0)
      output.append(r.convertToCartesian())
    return output
    
class PolarElipticalGen(PolarUniformGen):

  def __init__(self,scale=1):
    """
    Wraps the important parameters we want to vary to the initialization function.
    """
    self.r_dist = CircularElipticalGen(0,scale,1,0)

def box_mueller(scale):
  """
  Gets a random r and phi value distributed normally.
  This is the so called polar form and is faster for
  getting an entire distribution.
  """
  wrand=0.0
  while(wrand==0 or wrand>=1):
    urand, vrand = stats.uniform.rvs(size=2)
    urand=urand*2.0-1.0
    vrand=vrand*2.0-1.0
    wrand=urand**2.0+vrand**2.0
  wrand=np.sqrt(-2.0*np.log(wrand)/wrand)
  x=urand*wrand*scale
  y=vrand*wrand*scale
  return (x,y)

def box_mueller_ring(scale,min_r,max_r):
  """
  Gets a ring of random r and phi value from
  a normal distribution.
  This is the so called basic form and 
  does not require rejection at the expense
  of calculating trig functions.  However, since
  only a thin ring is desired, the rejection issue
  of the polar form is bad. 
  """
  pi = 3.14159
  scaled_min_r = min_r/scale
  scaled_max_r = max_r/scale
  max_u1 = np.exp(-0.5*scaled_min_r**2)
  min_u1 = np.exp(-0.5*scaled_max_r**2)
  diff_u1 = max_u1 - min_u1
  u1 = stats.uniform.rvs(scale=diff_u1)+min_u1
  u2 = stats.uniform.rvs()
  r = np.sqrt(-2*np.log(u1))
  theta = 2*pi*u2
  #sys.stderr.write(str(diff_u1)+":"+str(scaled_min_r)+"-"+str(scaled_max_r)+":"+str(r)+":"+str(scale)+"\n")
  return (r*np.cos(theta)*scale,r*np.sin(theta)*scale)
  

def probability_from_normal_distribution(steps,width=3):
  """
  Breaks up the normal distribution from -width * standard
  deviation to width * standard deviation into normalized
  probabilities at each step. 
  """
  total_probability = 0 #for normalization. 
  probability = []
  stepsize = 2.*width/steps
  for i in range(steps):
    current_probability = stats.norm.cdf(stepsize*(i+1) - width) - stats.norm.cdf(stepsize*i - width)
    probability.append(current_probability)
    total_probability += current_probability
  probability = [p/total_probability for p in probability]
  return probability


