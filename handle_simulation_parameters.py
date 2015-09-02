import sys
import os
import csv
import argparse
import math

class ParameterSet():
 
  def __init__(self,shape,n_particles,e_per_particle,applied_field,timesteps,lazer_duration=50):
    self.shape = shape
    if self.shape == 'G':
      self.shape_code = 0
    if self.shape == 'U':
      self.shape_code = 1
    if self.shape == 'E':
      self.shape_code = 2
    self.n_particles = n_particles
    self.e_per_particle = e_per_particle
    self.electron_model = str(n_particles) + "x" + str(e_per_particle)
    self.applied_field = applied_field
    self.lazer_duration = lazer_duration
    self.sigma_t = lazer_duration/( 2 * math.sqrt( 2 * math.log(2) ) )*pow(10,-15)
    self.timesteps = timesteps
    self.main_output_path = None
    self.stepnames = {}
    
  def __getitem__(self,key):
    return getattr(self,key)

  def generateMainOutputDirectoryName(self):
    """
    Converts the primary attributes to
    a directory name.
    """
    filepath = []
    filepath.append(str(self.shape))
    filepath.append(str(self.electron_model))
    filepath.append(str(self.applied_field + "MV-per-m"))
    filepath.append(str(self.timesteps + "_steps"))
    return "_".join(filepath)

  def createMainOutputPath(self,working_directory):
    if self.main_output_path is None:
      self.main_output_path = os.path.join(working_directory,self.generateMainOutputDirectoryName())
    if not os.path.exists(self.main_output_path) or not os.path.isdir(self.main_output_path):
      os.mkdir(self.main_output_path)

  def stepIsTracked(self,stepname,order = None):
    for key, value in self.stepnames.iteritems():
      if value == stepname:
        if order is None:
          return True
        else:
          if int(key) == int(order):
            return True
          else:
            raise Exception('Step ' + stepname + ' is already tracked with order ' + key + '. You are trying to add it with order ' + order + '.')
    return False

  def readStepnamesFromDir(self):
    if self.main_output_path is None:
      raise Exception('Main output path needs to exist before a step may be added.')
    for pathname in os.listdir(self.main_output_path):
      if os.path.isdir(os.path.join(self.main_output_path,pathname)):
        pieces = pathname.split('-')
        if len(pieces) == 2:
          self.stepnames[str(pieces[0])] = pieces[1] #0 is the number, 1 is the stepname 

  def addStepname(self,stepname):
    numbers = self.stepnames.keys()
    i = 1
    while True:
      if str(i) in numbers:
        i += 1
      else:
        break
    self.stepnames[str(i)] = stepname

  def returnStepnamePath(self,stepname):
    self.readStepnamesFromDir()
    if not self.stepIsTracked(stepname):
      self.addStepname(stepname)
    for key, value in self.stepnames.iteritems():
      if value == stepname:
        return os.path.join(self.main_output_path,str(key) + "-" + value)

  def createStepnamePath(self,stepname):
    stepname_path = self.returnStepnamePath(stepname)
    if not os.path.exists(stepname_path):
      os.mkdir(stepname_path)

  def returnInitialConditionsFilename(self):
    filename_pieces = ['IniCondition']
    filename_pieces.append(self.shape)
    filename_pieces.append(self.n_particles)
    filename_pieces.append(self.e_per_particle)
    self.initial_conditions_file = "-".join(filename_pieces) + ".txt"
    return self.initial_conditions_file

  def returnPEMainFilename(self):
    if self.shape == 'G':
      self.pulse_file = 'PE-Main_Gaussian.fox'
    if self.shape == 'U':
      self.pulse_file = 'PE-Main_Uniform.fox'
    if self.shape == 'E':
      self.pulse_file = 'PE-Main_Elliptical.fox'
    return self.pulse_file

  def returnHoleFilename(self):
    if self.shape == 'G':
      self.hole_file =  'PosHoleField-Gaussian.txt'
    if self.shape == 'U':
      self.hole_file = 'PosHoleField-Uniform.txt'
    if self.shape == 'E':
      self.hole_file = 'PosHoleField-Elliptical.txt'
    return self.hole_file

  def extractTotalParticlesFromIniConditionFile(self,filepath):
    with open(filepath) as f:
      for line in f:
        if "total particles generated:" in line:
          pieces = line.split("total particles generated:")
          self.total_particles_generated = pieces[1].rstrip().strip(" ")
    return

class ParameterSetList():
  
  def __init__(self):
    self.list = []

  def addParameterSet(self,shape,n_particles,e_per_particle,applied_field,timesteps,lazer_duration=50):
    self.list.append(ParameterSet(shape,n_particles,e_per_particle,applied_field,timesteps,lazer_duration))

  def injectParametersFromCsv(self,filepath,lazer_duration=50):
    """
    Reads the csv file containing all of the
    parameters we want to run and puts them into the
    object.
    """
    with open(filepath) as csvfile:
      parameters = csv.DictReader(csvfile)
      for row in parameters:
        shape = row["Shape"]
        n_particles = row["NParticles"]
        e_per_particle = row["ElectronsPerParticle"]
        applied_field = row["Fa"]
        timesteps = row["Timesteps"]
        self.addParameterSet(shape,n_particles,e_per_particle,applied_field,timesteps,lazer_duration)


if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package')
  args = parser.parse_args()

