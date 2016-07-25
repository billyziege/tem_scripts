import sys
import os
import math
import argparse

class CoordinateException(Exception):
  pass

class CoordinateVector():
  """
  Base class to handle methods applicable to all coordinate systems.
  """

  def __init__(self,order=[],values={}):
    """
    Order is the name of the coordinates and the order in which they will be referred, and
    values are the components along each coordinate for the vector.
    """
    if len(values.keys()) != len(order):
      raise CoordinateException("Expect all coordinates in a coordinate vector to have a value")
    self.setOrder(order)
    self.setVector(values)
      
  def setOrder(self,order):
    """
    Sets the order to the provided order list.
    """
    self.order = order

  def getOrder(self):
    """
    Gets the order of the vector.
    """
    return self.order

  def setVector(self,values):
    """
    Sets the coordinate as an attribute to the provided value.
    """
    for coordinate in self.order:
      setattr(self,coordinate,values[coordinate])

  def getVector(self):
    """
    Get the values of the vector in the set order.
    """
    vector = []
    for coordinate in self.order:
      vector.append(getattr(self,coordinate))
    return vector

  def __str__(self):
    """
    Returns the coordinate vector in the specified order separated by the given delimiter.
    """
    return " ".join([str(x) for x in self.getVector()])

  def __abs__(self):
    """
    Returns the L2 norm of the vector.
    """
    return math.sqrt(self * self)

  def __len__(self):
    """
    Returns the number of coordinates in the vector.
    """
    return int(len(self.order))

  def __add__(self,other):
    """
    Defines the meaning of addition.
      1. Other is scalar, add scalar to each component returning vector.
      2. Other is vector, add by components returning vector.
    """
    if isinstance(other,float) or isinstance(other,int):
      values = {}
      vector = self.getVector()
      for i in range(len(self)):
        values[self.order[i]] = other+vector[i]
      instance = CoordinateVector(self.order,values)
      instance.__class__ = self.__class__
      return instance
    if isinstance(other,CoordinateVector):
      if len(other) != len(self):
        raise CoordinateException("Error: Attempting to add two vectors of different lengths.")
      values = {}
      vector = self.getVector()
      other_vector = other.getVector()
      for i in range(len(self)):
        values[self.order[i]] = other_vector[i]+vector[i]
      instance = CoordinateVector(self.order,values)
      instance.__class__ = self.__class__
      return instance
    raise CoordinateException("Using addition for an unknown class.")
    
  def __iadd__(self,other):
    """
    Defines the meaning of +=.
    """
    return self.__add__(other)  

  def __sub__(self,other):
    """
    Defines the meaning of subtraction.
      1. Other is scalar, subtract scalar from each component returning vector.
      2. Other is vector, subtract by components returning vector.
    """
    return self.__add__(-1*other)

  def __isub__(self,other):
    """
    Defines the meaning of -=.
    """
    return self.__sub__(other)  

  def __rmul__(self,other):
    """
    Defines the meaning of multiplication.
      1.  Other is scaler, scales vector returning vector.
      2.  Other is vector, dot product returning scalar.
    """
    if isinstance(other,float) or isinstance(other,int):
      values = {}
      vector = self.getVector()
      for i in range(len(self)):
        values[self.order[i]] = other*vector[i]
      instance = CoordinateVector(self.order,values)
      instance.__class__ = self.__class__
      return instance
    if isinstance(other,CoordinateVector):
      if len(other) != len(self):
        raise CoordinateException("Error: Attempting to multiply two vectors of different lengths.")
      sqsum = 0
      vector = self.getVector()
      other_vector = other.getVector()
      for i in range(len(self)):
        sqsum += other_vector[i]*vector[i]
      return sqsum
    raise CoordinateException("Using multiply for an unknown class.")

  def __div__(self,other):
    """
    Defines the meaning of division.  Only works when other is a scalar.
      1.  Other is scaler, scales vector returning vector.
    """
    if isinstance(other,float) or isinstance(other,int):
      values = {}
      vector = self.getVector()
      for i in range(len(self)):
        values[self.order[i]] = vector[i]/other
      instance = CoordinateVector(self.order,values)
      instance.__class__ = self.__class__
      return instance
    raise CoordinateException("Using division for an unknown class.")
      
  def copy(self):
    """
    Returns a duplicate object with the same values but different reference locations.
    """
    return 1*self
      

