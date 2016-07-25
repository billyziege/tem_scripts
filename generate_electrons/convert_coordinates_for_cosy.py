import sys
import os
import math
import argparse
import copy
import numpy as np
try:
  import cPickle as pickle
except:
  import pickle
from scipy import constants
from phase_volume_6d import ParticlePhaseCoordinates
from coordinate_vector_3d import Cartesian3DVector

class TimeBin():
  """
  A class to provide access to lists of timepoints that will
  enter the simulaton at the same time.
  """
  
  def __init__(self,time_points):
    self.time_points = time_points
    #time_points.sort(key=float)
    self.bin_time = 0
    self.number_points = len(time_points)
    self.dt = 0

  def __iter__(self):
    """
    Allows the ability to iterate over the time points in
    the time bin without referencing the time_points attribute.
    """
    for time_point in self.time_points:
      yield time_point

  def __eq__(self,other):
    """
    Two time bins are equal if they have the same bin time.
    """
    if not isinstance(other,TimeBin):
      return False
    return self.bin_time == other.bin_time

  def __add__(self,other):
    """
    Two time bins can be added only if they have the same bin time...
    then there time points and number_points are added.
    """
    if not isinstance(other,TimeBin):
      raise Exception("Trying to add a non-time bin object to a time bin.")
    if not self == other:
      raise Exception("To add time bins, their bin times need to be equal.")
    return TimeBin(self.time_points + other.time_points)

  def __iadd__(self,other):
    """
    Allows the use of +=.
    """
    return self.__add__(other)
     
class TimeBins():
  """
  A class to handle all of the time bins together.
  """

  def __init__(self,time_points=[],number_bins=0,max_time=0):
    self.time_bins = []
    self.max_time = max_time
    if number_bins > 0:
      self.assignBins(time_points,number_bins)
      self = self.resolveBins()
      self.calcDt()
    
  def __len__(self):
    """
    Returns the number of bins.
    """
    return len(self.time_bins)

  def __iter__(self):
    """
    Provides access to the time bins without reference to the 
    time_bins attribute.
    """
    for time_bin in self.time_bins:
      yield time_bin

  def addBin(self,time_points):
    """
    Adds the list of time points as a time bin to
    the time bins attribute.
    """
    if isinstance(time_points,TimeBin):
      self.time_bins.append(time_points)
    else:
      self.time_bins.append(TimeBin(time_points))

  def assignBins(self,time_points,number_bins,even_time=True):
    """
    Splits the time_sample list into a list of lists where
    each list has len(time_sample)/time_steps elements in it or
    the time steps are the same.
    """
    time_points = sorted(time_points,key=float)
    number_time_points = len(time_points)
    max_time_step = float(self.max_time)/float(number_bins)
    if even_time:
      for bin_number in range(number_bins):
        current_bin = []
        bin_time = max_time_step*(bin_number+1)
        try:
          while float(time_points[0]) <= float(bin_time):
            current_time_point = time_points.pop(0)
            current_bin.append(current_time_point)
        except:
          pass
        self.addBin(current_bin)
        self.time_bins[-1].bin_time = bin_time
    else:
      average_number_time_points = float(number_time_points)/float(number_bins) #This need not be an integer
      last_number_time_points = 0
      last_time = 0
      while last_number_time_points < number_time_points:
        current_time_point = time_points.pop(0)
        current_bin = []
        current_bin.append(current_time_point)
        #while len(current_bin) < np.floor(average_number_time_points) and current_time_point-last_time < max_time_step and last_number_time_points+len(current_bin) < number_time_points:  
        while len(current_bin) < np.floor(average_number_time_points) and last_number_time_points+len(current_bin) < number_time_points:  
          current_time_point = time_points.pop(0)
          current_bin.append(current_time_point)
        self.addBin(current_bin)
        last_number_time_points += len(current_bin)
        last_time = current_time_point

  def resolveBins(self):
    """
    Goes through and determines if neighboring time bins have the same
    bin time.  If so, they are joined.
    """
    same_count = 0
    resolved_bins = TimeBins(max_time=self.max_time)
    for i in range(len(self)-1):
      while same_count > 0:
        same_count -= 1
        continue
      current_time_bin = self.time_bins[i]
      count = 1
      next_time_bin = self.time_bins[i+count]
      while i<len(self) and next_time_bin == current_time_bin:
        current_time_bin += next_time_bin
        same_count += 1
        count += 1
        next_time_bin = self.time_bins[i+count]
      resolved_bins.addBin(current_time_bin)
    return resolved_bins

  def calcDt(self):
    """
    Determines the dt between the current time bin and the next.  If
    it is the last time, the max time is used as "next".
    """
    if len(self.time_bins) == 0:
      return
    current_time_bin = self.time_bins[0]
    for i in range(1,len(self.time_bins)):
      next_time_bin = self.time_bins[i]
      current_time_bin.dt = next_time_bin.bin_time - current_time_bin.bin_time
      current_time_bin = next_time_bin
    current_time_bin.dt = self.max_time - current_time_bin.bin_time
    
         
def connect_phase_and_write_output(time_bins,output_file,coordinates_table,mass,turn_off_advance=False):
  """
  Takes a list of lists of time points, treats them as bins
  with the time of the minimum time point, connects the phase coordinates, and 
  writes the output to the output_file.
  """
  with open(output_file, 'w') as output:
    for time_bin in time_bins:
      dt = time_bin.dt
      front_output = [str(time_bin.bin_time),str(dt),str(time_bin.number_points)]
      for row in coordinates_table:
        if row["t"] not in time_bin:
          continue
        x = Cartesian3DVector(row["x"],row["y"],row["z"]) 
        setattr(x,"z",1e-20) #Make sure the particle is not sitting at 0.
        p = Cartesian3DVector(row["px"],row["py"],row["pz"]) 
        particle = ParticlePhaseCoordinates(mass,x=x,p=p)
        if not turn_off_advance:
          advanced_position = particle.advancePosition(time_bin.bin_time - row["t"])
          particle.setPosition(advanced_position)
        output.write(" ".join(front_output) + " " + str(particle) + "\n")

def read_coordinates_file(coordinates_file):
  """
  Reads the time and particle coordinates from a pickled dict and returns a dict. 
  """
  with open(coordinates_file, 'r') as f:
    data = pickle.load(f)
  return data

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generates electron macroparticle phase coordinates in the format needed for COSY.')
  parser.add_argument('coordinates_file', type=str, help='Path to the coordinates file.')
  parser.add_argument('time_steps', type=int, help='The number of timesteps desired (estimate).')
  parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is IniConditions.txt in the current working directory.', default="IniConditions.txt")
  parser.add_argument('-n','--number_of_electrons_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons simulated by each macroparticle in the COSY simulation. Default is 100.', default=100)
  parser.add_argument('--turn_off_advance', dest="turn_off_advance", action="store_true", help='Turn off the advancing of the electrons from the cathode.', default=False)
  args = parser.parse_args()

  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_particle = args.number_of_electrons_per_macroparticle*mass_of_electron
  coordinates_table = read_coordinates_file(args.coordinates_file)
  time_points = [row["t"] for row in coordinates_table]
  time_points = sorted(time_points,key=float)
  time_bins = TimeBins(time_points,args.time_steps,time_points[-1])
  connect_phase_and_write_output(time_bins,args.output_file,coordinates_table,mass_of_particle,args.turn_off_advance)
  


