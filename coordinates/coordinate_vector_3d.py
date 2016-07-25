import sys
import os
import math
import numpy as np
import argparse
from coordinate_vector import CoordinateException, CoordinateVector

class Cartesian3DVector(CoordinateVector):
  """
  The standard x, y, z vector.
  """

  def __init__(self,x=0,y=0,z=0):
    """
    Initializes the vector by default to 0,0,0
    """
    order = ["x","y","z"]
    values = {"x":float(x), "y":float(y), "z":float(z)}
    CoordinateVector.__init__(self,order,values)

  def convertToCylindrical(self):
    """
    Returns a cylindrical coordinate system vector with the cooridnates
    rho, phi, z that correspond to the same x, y, and z values.
    """
    rho = math.sqrt(self.x**2 + self.y**2)
    if rho == 0:
      return Cylindrical3DVector(0,0,self.z)
    if self.y >= 0:
      phi = acos(self.x/rho)
    else:
      phi = 2*np.pi - math.acos(self.x/rho)
    return Cylindrical3DVector(rho,phi,self.z)

  def convertToSpherical(self):
    """
    Returns a spherical coordinate system vector with the cooridnates
    r, phi, and theta that correspond to the same x, y, and z values.
    """
    rhosq = self.x**2 + self.y**2
    rho = math.sqrt(rhosq)
    r = math.sqrt(rhosq + self.z**2)
    if r == 0:
      return Spherical3DVector()
    theta = math.acos(self.z/r)
    if self.y >= 0:
      phi = acos(self.x/rho)
    else:
      phi = 2*np.pi - math.acos(self.x/rho)
    return Spherical3DVector(r,phi,theta)

  def __cross_product__(self,other):
    """
    Defines the cross product in 3 dimensions.
    """
    if isinstance(other,CoordinateVector):
      if len(other) != 3:
        raise CoordinateException("Error: Attempting to do the cross product on a vector that is not 3 dimensions.")
      values = {}
      values["x"] = self.y*other.z - self.z*other.y
      values["y"] = self.z*other.x - self.x*other.z
      values["z"] = self.x*other.y - self.y*other.x
      return CoordinateVector(self.order,values)
    raise CoordinateException("Using the cross product for an unknown class.")

class Cylindrical3DVector(CoordinateVector):
  """
  The standard rho, phi, z vector.
  """

  def __init__(self,rho=0,phi=0,z=0):
    """
    Initializes the vector by default to 0,0,0
    """
    order = ["rho","phi","z"]
    values = {"rho":float(rho), "phi":float(phi), "z":float(z)}
    CoordinateVector.__init__(self,order,values)

  def convertToCartesian(self):
    """
    Returns a cartesian coordinate system vector with the cooridnates
    x, y, and z that correspond to the same rho, phi, and z values.
    """
    x = self.rho * math.cos(self.phi)
    y = self.rho * math.sin(self.phi)
    return Cartesian3DVector(x,y,self.z)

  def convertToSpherical(self):
    """
    Returns a spherical coordinate system vector with the cooridnates
    r, phi, and theta that correspond to the same rho, phi, and z values.
    """
    r = math.sqrt(self.rho**2 + self.z**2)
    if r == 0:
      return Spherical3DVector(0,self.phi,0)
    theta = math.acos(self.z/r)
    return Spherical3DVector(r,self.phi,theta)

class Spherical3DVector(CoordinateVector):
  """
  The standard r, phi, and theta vector.
  """

  def __init__(self,r=0,phi=0,theta=0):
    """
    Initializes the vector by default to 0,0,0
    """
    order = ["r","phi","theta"]
    values = {"r":float(r), "phi":float(phi), "theta":float(theta)}
    CoordinateVector.__init__(self,order,values)

  def convertToCartesian(self):
    """
    Returns a cartesian coordinate system vector with the cooridnates
    x, y, and z that correspond to the same r, phi, and theta values.
    """
    x = self.r * math.cos(self.phi)*math.sin(self.theta)
    y = self.r * math.sin(self.phi)*math.sin(self.theta)
    z = self.r*math.cos(self.theta)
    return Cartesian3DVector(x,y,z)

  def convertToCylindrical(self):
    """
    Returns a cylindrical coordinate system vector with the cooridnates
    rho, phi, and z that correspond to the same r, phi, and theta values.
    """
    z = self.r*math.cos(self.theta)
    rho = math.sqrt(self.r**2 -z**2)
    return Cylindrical3DVector(rho,self.phi,z) 
