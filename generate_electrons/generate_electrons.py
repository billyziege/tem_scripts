import sys
import os
import math
import argparse
import copy
import numpy as np
from scipy import constants
from scipy import stats
from phase_volume_6d import Phase6DVolume
from coordinate_vector import CoordinateVector
from coordinate_vector_3d import Cartesian3DVector, Cylindrical3DVector, Spherical3DVector
import distributions

class MyDict(dict):

  def __missing__(self,key):
    return None

class ProximityGrid():
  """
  A class to quickly handle checking to make sure a new 
  particle is not too close to other particles.
  """

  def __init__(self,stepsize,order=["x","y","z"]):
    """
    Initializes the grid and the stepsize.
    """
    self.grid = MyDict() #Not perfect since grid[key] may have value, but okay.
    self.order = order #order of the dimensions.
    self.stepsize = stepsize

  def getVectorFromGridPoint(self,grid_point_vector):
    """
    Returns the vector associated with the grid point if any.
    """
    current_grid_level = self.grid
    for coordinate in self.order:
      current_grid_level = current_grid_level[str(getattr(grid_point_vector,coordinate))]
      if current_grid_level is None:
        break
    return current_grid_level

  def getGridPointNeighborhood(self,grid_point_vector):
    """
    Returns the vectors associated with the grid points within 1 step of the povided grid point.
    """
    neighborhood_grid_point_vectors = []
    neighborhood_grid_point_vectors.append(grid_point_vector)
    neighborhood_vectors =[]
    current_grid_level = self.grid
    for coordinate in self.order:
      for vector in neighborhood_grid_point_vectors:
        new_neighborhood = []
        current_coordinate_value = int(getattr(vector,coordinate))
        current_coordinate_value_plus_one = current_coordinate_value + 1
        current_coordinate_value_minus_one = current_coordinate_value + 1
        vector_plus_one = vector.copy()
        vector_minus_one = vector.copy()
        setattr(vector_plus_one,coordinate,current_coordinate_value_plus_one)
        setattr(vector_minus_one,coordinate,current_coordinate_value_minus_one)
        true_vector_plus_one = self.getVectorFromGridPoint(vector_plus_one)
        if true_vector_plus_one is not None:
          new_neighborhood.append(vector_plus_one)
          neighborhood_vectors.append(true_vector_plus_one)
        true_vector_minus_one = self.getVectorFromGridPoint(vector_minus_one)
        if true_vector_minus_one is not None:
          new_neighborhood.append(vector_minus_one)
          neighborhood_vectors.append(true_vector_minus_one)
      neighborhood_grid_point_vectors.extend(new_neighborhood)
    return neighborhood_vectors

  def checkVector(self,input_vector):
    """
    Returns true if only no other vector in the grid is within the grid stepsize.
    """
    grid_point_vector = self.convertVectorToGridPointVector(input_vector)
    vector = self.getVectorFromGridPoint(grid_point_vector)
    if vector is not None:
      return False
    for vector in self.getGridPointNeighborhood(grid_point_vector):
      if abs(input_vector-vector) < self.stepsize:
        return False
    return True
    

  def addGridPoint(self,input_vector):
    """
    Adds the provided cartesian 3d vector to the nearest grid point.
    """
    grid_point_vector = self.convertVectorToGridPointVector(input_vector)
    value = MyDict()
    previous_coordinate = self.order[-1]
    coordinate_value = str(getattr(grid_point_vector,previous_coordinate))
    value[coordinate_value] = input_vector
    for coordinate in reversed(self.order):
      if coordinate == previous_coordinate:
        continue
      coordinate_value = str(getattr(grid_point_vector,coordinate))
      previous_coordinate_value = str(getattr(grid_point_vector,previous_coordinate))
      value[coordinate_value] = copy.deepcopy(value)
      if previous_coordinate_value != coordinate_value:
        del value[previous_coordinate_value]
      previous_coordinate = coordinate
    self.grid.update(value)
    
  def convertVectorToGridPointVector(self,input_vector):
    """
    Finds the closest grid point to each of the components of the 
    input coordinate vector.
    """
    values = {}
    if set(input_vector.order) != set(self.order):
      raise Exception("Trying to convert a vector ont to a grid, but the dimensions are different.")
    for coordinate in self.order:
      values[coordinate] = int(np.round(getattr(input_vector,coordinate)/self.stepsize + 0.5))
    output_vector = CoordinateVector(input_vector.order,values)
    return output_vector
    

class TimeBin():
  """
  A class to provide access to lists of timepoints that will
  enter the simulaton at the same time.
  """
  
  def __init__(self,time_points):
    self.time_points = time_points
    #time_points.sort(key=float)
    self.bin_time = time_points[-1]
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

  def assignBins(self,time_points,number_bins):
    """
    Splits the time_sample list into a list of lists where
    each list has len(time_sample)/time_steps elements in it.
    """
    time_points = sorted(time_points,key=float)
    number_time_points = len(time_points)
    average_number_time_points = float(number_time_points)/float(number_bins) #This need not be an integer
    max_time_step = float(self.max_time)/float(number_bins)
    last_number_time_points = 0
    last_time = 0
    while last_number_time_points < number_time_points:
      current_time_point = time_points[0]
      del time_points[0]
      current_bin = []
      current_bin.append(current_time_point)
      #while len(current_bin) < np.floor(average_number_time_points) and current_time_point-last_time < max_time_step and last_number_time_points+len(current_bin) < number_time_points:  
      while len(current_bin) < np.floor(average_number_time_points) and last_number_time_points+len(current_bin) < number_time_points:  
        current_time_point = time_points[0]
        del time_points[0]
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
    
         
def generate_phase_and_write_output(time_bins,output_file,polar_distribution,proximity_threshold,turn_off_advance=False,time_shift=0,**kwargs):
  """
  Takes a list of lists of time points, treats them as bins
  with the time of the minimum time point, generate macroparticles
  using a photoemission process a the time point, and distributes the particles according to
  the provided distribution_parameters.   Writes the output to the output_file.
  """
  last_time = 0
  with open(output_file, 'w') as output:
    for time_bin in time_bins:
      dt = time_bin.dt
      front_output = [str(time_bin.bin_time),str(dt),str(time_bin.number_points)]
      proximity_grid = ProximityGrid(proximity_threshold)
      phase_volume = simulate_photoemission_processes(time_bin.number_points,**kwargs)
      for i in range(time_bin.number_points):
        particle = phase_volume.particles[i]
        random_position = polar_distribution.generateCartesian3DVector(1).pop()
        #sys.stderr.write("generated a random position.")
        particle.setPosition(random_position)
        if turn_off_advance:
          while proximity_grid.checkVector(random_position) is False: #Rechoose the x,y coordinates if the point fails the proximity check
            print "Proximity alert."
            random_position = polar_distribution.generateCartesian3DVector(1).pop()
            particle.setPosition(random_position)
          proximity_grid.addGridPoint(random_position)  
          particle.setPosition(random_position)
        else:
          advanced_position = particle.advancePosition(time_bin.bin_time - time_bin.time_points[i] + time_shift)
          while proximity_grid.checkVector(advanced_position) is False: #Rechoose the x,y coordinates if the point fails the proximity check
            print "Proximity alert."
            random_position = polar_distribution.generateCartesian3DVector(1).pop()
            particle.setPosition(random_position)
            advanced_position = particle.advancePosition(time_bin.bin_time - time_bin.time_points[i])
          proximity_grid.addGridPoint(advanced_position)  
          particle.setPosition(advanced_position)
        output.write(" ".join(front_output) + " " + str(particle) + "\n")
      last_time=time_bin.bin_time

def generate_jenni_phase_and_write_output(probability,dt,total_number_particles,output_file,polar_distribution,proximity_threshold,**kwargs):
  """
  Recreates Jenni's algorithm for electron generation.
  """
  current_time=0
  with open(output_file, 'w') as output:
    for p in range(len(probability)):
      current_time += dt
      number_of_particles = int(np.round(total_number_particles*probability[p]))
      front_output = [str(current_time),str(dt),str(number_of_particles)]
      proximity_grid = ProximityGrid(proximity_threshold)
      phase_volume = simulate_photoemission_processes(number_of_particles,jennis_dist=True,**kwargs)
      for i in range(number_of_particles):
        particle = phase_volume.particles[i]
        random_position = polar_distribution.generateCartesian3DVector(1).pop()
        particle.setPosition(random_position)
        while proximity_grid.checkVector(random_position) is False: #Rechoose the x,y coordinates if the point fails the proximity check
          print "Proximity alert."
          random_position = polar_distribution.generateCartesian3DVector(1).pop()
          particle.setPosition(random_position)
        proximity_grid.addGridPoint(random_position)  
        particle.setPosition(random_position)
        output.write(" ".join(front_output) + " " + str(particle) + "\n")

def simulate_photoemission_processes(n,e_fermi=5.000E-6,e_work=4.45E-6,e_photon=4.65E-6,number_of_electrons_per_macroparticle=100,jennis_dist=False,uniform_momentum=False,**kwargs):
  """
  Simulates the photoemission process n times resulting in n particles with the simulated momentum at the (0,0,0) position.
  All variables beginning with e_ are energies in MeV units.
  """
  phase_volume = Phase6DVolume()
  mass_of_electron = constants.physical_constants["electron mass energy equivalent in MeV"][0]
  mass_of_particle = number_of_electrons_per_macroparticle*mass_of_electron
  for i in range(n):
    #e_emitted_electron=(stats.uniform.rvs())*(e_photon-e_work)+e_fermi+e_work
    if uniform_momentum:
      velocity = Cartesian3DVector(z=5.2e-4)
    else:
      e_emitted_electron=np.sqrt(stats.uniform.rvs())*(e_photon-e_work)+e_fermi+e_work
      if jennis_dist:
        e_emitted_electron=stats.uniform.rvs()*(e_photon-e_work)+e_fermi+e_work
      while e_emitted_electron == e_fermi + e_work: #Prevent problems since we need to have some energy beyond E_fermi and E_work.
        e_emitted_electron=np.sqrt(stats.uniform.rvs())*(e_photon-e_work)+e_fermi+e_work
        if jennis_dist:
          e_emitted_electron=stats.uniform.rvs()*(e_photon-e_work)+e_fermi+e_work
      theta_max=math.acos(math.sqrt((e_fermi+e_work)/e_emitted_electron))
      theta = math.acos(1 - (1 - math.sqrt((e_fermi+e_work)/e_emitted_electron))*stats.uniform.rvs())
      if jennis_dist:
        theta = theta_max*stats.uniform.rvs()
      while theta == theta_max: #Prevent problems since we are on [0,theta_max)
        theta = math.acos(1 - (1 - math.sqrt((e_fermi+e_work)/e_emitted_electron))*stats.uniform.rvs())
        if jennis_dist:
          theta = theta_max*stats.uniform.rvs()
      phi = stats.uniform.rvs(scale=2*np.pi)
      while phi == 2*np.pi: #Prevent double counting since we are on [0,2*pi)
        phi = stats.uniform.rvs(scale=2*np.pi)
      r = math.sqrt(2*e_emitted_electron/mass_of_electron) #maginitude of velocity
      velocity = Spherical3DVector(r=r,phi=phi,theta=theta).convertToCartesian()
      velocity.z = math.sqrt(velocity.z*velocity.z - 2*(e_fermi+e_work)/mass_of_electron)
    phase_volume.addParticle(mass_of_particle,v=velocity)
  return phase_volume
  
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generates electron macroparticle phase coordinates in the format needed for COSY.')
  parser.add_argument('number_of_macroparticles', type=int, help='The number of macroparticles desired.')
  parser.add_argument('time_steps', type=int, help='The number of timesteps desired (estimate).')
  parser.add_argument('-n','--number_of_electrons_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons simulated by each macroparticle in the COSY simulation. Default is 100.', default=100)
  parser.add_argument('-t','--std-t', dest="std_t", type=float, help='The standard deviation of the gaussian of electron emission in time. Default is about 21 fs.', default=50./( 2. * math.sqrt( 2. * math.log(2.) ) )*pow(10,-15))
  parser.add_argument('-r','--std-r', dest="std_r", type=float, help='The standard deviation of the gaussian of electron emission in radial components. Default is about 96 um.', default=math.sqrt(81*115)*1e-6)
  parser.add_argument('-p','--emission_profile', dest="profile", type=str, help='Describes the profile of the electrons in the radial direction.  Options are Gaussian (G), Elliptical (E), Uniform (U), or Annular (A).  Default is Gaussian.', default="Gaussian")
  parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is IniConditions.txt in the current working directory.', default="IniConditions.txt")
  parser.add_argument('--time-upper-bound', dest="t_upper_bound", type=float, help='The multiple of the std_t above the mean that will be permitted from the distribution.  Default is 3.', default=3.)
  parser.add_argument('--time-lower-bound', dest="t_lower_bound", type=float, help='The multiple of the std_t below the mean that will be permitted from the distribution.  Default is -3.', default=-3.)
  parser.add_argument('--time-shift', dest="time_shift", type=float, help='The amount of time the particle advances before injection into the simulation.  Default is zero.', default=0.)
  parser.add_argument('--fermi_energy', dest="e_fermi", type=float, help='The fermi energy of the electron source. Default is 5.000E-6.', default=5.000E-6)
  parser.add_argument('--work_potential', dest="e_work", type=float, help='The work potential of the electron source. Default is 4.450E-6.', default=4.450E-6)
  parser.add_argument('--photon_energy', dest="e_photon", type=float, help='The energy of the lazer\'s photon. Default is 4.650E-6.', default=4.650E-6)
  parser.add_argument('--proximity_threshold', dest="proximity_threshold", type=float, help='The minimum allowed distance between generated particles.', default=1E-10)
  parser.add_argument('--turn_off_advance', dest="turn_off_advance", action="store_true", help='Turn off the advancing of the electrons from the cathode.', default=False)
  parser.add_argument('--jennis_dist', dest="jennis_dist", action="store_true", help='Use Jennis distribution instead of my corrected.', default=False)
  parser.add_argument('--uniform_momentum',dest='uniform_momentum',action='store_true',help='Enforces that all particles have the same, average momentum when generated.  Default is off.', default=False)
  parser.add_argument('--min_r',dest='min_r',type=float,help='Used with the annular distribution.  Specifies the inner radius of the Annulus.  Default is 0.', default=0)
  parser.add_argument('--max_r',dest='max_r',type=float,help='Uased with the annular distribution.  Specifies the outer radius of the Annulus.  Default is 1.', default=1)
  args = parser.parse_args()

  photoemission_parameters = {}
  photoemission_parameters["number_of_electrons_per_macroparticle"] = args.number_of_electrons_per_macroparticle
  photoemission_parameters["e_fermi"] = args.e_fermi
  photoemission_parameters["e_work"] = args.e_work
  photoemission_parameters["e_photon"] = args.e_photon

  #args.std_r = args.std_r/np.sqrt(args.number_of_electrons_per_macroparticle)   

  if args.profile == 'G' or args.profile == 'Gaussian':  
    polar_distribution = distributions.PolarGaussianGen(scale=args.std_r)
  elif args.profile == 'E' or args.profile == 'Elliptical':  
    polar_distribution = distributions.PolarEllipticalGen(scale=math.sqrt(5)*args.std_r)
  elif args.profile == 'U' or args.profile == 'Uniform':  
    polar_distribution = distributions.PolarUniformGen(scale=2*args.std_r)
  elif args.profile == 'A' or args.profile == 'Annular':  
    polar_distribution = distributions.AnnularGaussianGen(scale=2*args.std_r,upper=args.max_r,lower=args.min_r)
  else:
    raise Exception("An unsupported profile has been provided.")

  #sys.stderr.write("Passed distribution intialization") 
  max_time = (args.t_upper_bound-args.t_lower_bound)*args.std_t
  if args.jennis_dist:
    probability = distributions.probability_from_normal_distribution(args.time_steps)
    generate_jenni_phase_and_write_output(probability,6.*args.std_t/args.time_steps,args.number_of_macroparticles,args.output_file,polar_distribution,photoemission_parameters,args.proximity_threshold)
  else:
    time_dist = distributions.GaussianGen(-1*args.t_lower_bound,args.std_t,args.t_upper_bound,args.t_lower_bound)
    time_sample = time_dist.getSample(args.number_of_macroparticles)
    time_sample.sort(key=float)
    time_bins = TimeBins(time_sample,args.time_steps,max_time=max_time)
    generate_phase_and_write_output(time_bins,args.output_file,polar_distribution,args.proximity_threshold,jennis_dist=args.jennis_dist,turn_off_advance=args.turn_off_advance,uniform_momentum=args.uniform_momentum,time_shift=args.time_shift,**photoemission_parameters)
  


