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

def uniform_line(dx,z,y,r_max,middle_point=True):
  """
  Return coordinates that have equidistance to their nearest neighbors in 1D.  
  All are lieing on a line of y=y and z=z.  
  """
  points = []
  if middle_point:
    x = 0
    if x**2 + y**2 < r_max**2:
      points.append(Cartesian3DVector(x=x,y=y,z=z))
  else:
    x = dx/2.
    if x**2 + y**2 < r_max**2:
      points.append(Cartesian3DVector(x=x,y=y,z=z))
      points.append(Cartesian3DVector(x=-x,y=y,z=z))
  x += dx
  while x**2 + y**2 < r_max**2:
    points.append(Cartesian3DVector(x=x,y=y,z=z))
    points.append(Cartesian3DVector(x=-x,y=y,z=z))
    x += dx
  return points
  
def uniform_plane(d,z,r_max,y_offset = 0):
  """
  Returns coordinates for a plane where nearest neighbors are equidistant.
  All points have coordinate z=z.
  """
  points = []
  #Do all the lines with a middle point.
  points.extend(uniform_line(d,z,y_offset,r_max,True))
  y_upper = y_offset + d
  y_lower = y_offset - d
  while y_upper < r_max or abs(y_lower) < r_max:
    points.extend(uniform_line(d,z,y_upper,r_max,True))
    points.extend(uniform_line(d,z,y_lower,r_max,True))
    y_upper += d
    y_lower -= d
  #Now do all the lines with no middle point.
  y_upper = y_offset + d/2.
  y_lower = y_offset - d/2.
  while y_upper < r_max or abs(y_lower) < r_max:
    points.extend(uniform_line(d,z,y_upper,r_max,False))
    points.extend(uniform_line(d,z,y_lower,r_max,False))
    y_upper += d
    y_lower -= d
  return points

def uniform_cylinder(d,r,l,z_min):
  """
  Returns coordinates for a cylinder where nearest neighbors are equidistant (d).
  Lenght of the cylinder is l and radius is r.
  """
  points = []
  z = z_min
  y_offset = 0.0
  while z < z_min + l:
    points.extend(uniform_plane(d,z,r,y_offset))
    y_offset = (y_offset + d/4.) % d
    z += d/2.
  return points

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Generates electron macroparticle phase coordinates in the format needed for COSY.')
  parser.add_argument('number_of_macroparticles', type=int, help='The number of macroparticles desired.')
  parser.add_argument("chirp", type=float, help='The longitudinal chirp.')
  parser.add_argument("chirpx", type=float, help='The transverse chirp.')
  #parser.add_argument('-n','--number_of_electrons_per_macroparticle', dest="number_of_electrons_per_macroparticle", type=int, help='The number of electrons simulated by each macroparticle in the COSY simulation. Default is 100.', default=100)
  parser.add_argument('-r','--radius', dest="radius", type=float, help='The standard deviation of the radius of the pulse.  Units in mm.  Default is 1 mm.', default=1.)
  #parser.add_argument('-p','--emission_profile', dest="profile", type=str, help='Describes the profile of the electrons in the radial direction.  Options are Gaussian (G), Elliptical (E), Uniform (U), or Annular (A).  Default is Gaussian.', default="Gaussian")
  parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is IniConditions.txt in the current working directory.', default="IniConditions.txt")
  parser.add_argument('--proximity_threshold', dest="proximity_threshold", type=float, help='The minimum allowed distance between generated particles.', default=20E-9)
  parser.add_argument('-e','--emittance', dest="ez", type=float, help='The longitudinal emittance in units of mm * mrad.  Default is 0.1.', default=0.1)
  parser.add_argument('--transverse-emittance', dest="er", type=float, help='The transverse emittance in units of mm * mrad.  Default is 0.1.', default=0.1)
  parser.add_argument('--dt', dest="dt", type=float, help='The longitudinal standard deviation of the pulse.  Units are in ps.  Default is 2.5 ps.', default=2.5)
  parser.add_argument('--vz', dest='vz', type=float, help='The speed of the pulse in m/s.  Default corresponds to a 100 keV electron pulse.', default=1.64e8)
  parser.add_argument('--dvz', dest='dvz', type=float, help='The standard deviation of the speed of the pulse in m/s.  Specify instead of emittance in z directon.', default=None)

  args = parser.parse_args()

  speed_light = constants.physical_constants["speed of light in vacuum"][0]
  beta = args.vz/speed_light
  gamma = np.sqrt(1/(1 - beta**2))

  polar_distribution = distributions.PolarUniformGen(scale=2*args.radius*1e-3)

  dz = args.vz*args.dt*1e-12#Standard deviation of z in units of m
  z_max = dz #Assuming uniform distribution
  z_min = -dz
  #dalphaz = args.ez/dz*1e-3 #The 1e-3 is due to the ez having units mm while dz is in m
  #dvz = dalphaz*speed_light/gamma*1e-3
  if args.dvz is None:
    dvz = speed_light*args.ez/(np.sqrt(6)*dz)*1e-6 #1e-3 due to mm*mrad
  else:
    dvz = args.dvz
  dvx = speed_light*args.er/(args.radius)*1e-3 #Both er and radius are in mm, so no conversion there.  The 1e-3 is due to mrad.
  #dvx = dalphax*speed_light/gamma
  #This next bit does random particle generation.
  for i in range(0,args.number_of_macroparticles):
    random_position = polar_distribution.generateCartesian3DVector(1).pop()
    z = stats.uniform.rvs(size=1,loc=-np.sqrt(3)*dz,scale=2*np.sqrt(3)*dz)
    random_position.z = z[0]

    vxmin = random_position.x*args.chirpx - dvx*np.sqrt(3)
    vx = stats.uniform.rvs(size=1,loc=vxmin,scale=2*np.sqrt(3)*dvx)

    vymin = random_position.y*args.chirpx - dvx*np.sqrt(3)
    vy = stats.uniform.rvs(size=1,loc=vymin,scale=2*np.sqrt(3)*dvx)

    vzmin = args.vz + z*args.chirp - dvz*np.sqrt(3)
    vz = stats.uniform.rvs(size=1,loc=vzmin,scale=dvz*np.sqrt(3)*2)
    random_velocity = Cartesian3DVector(vx,vy,vz)

    print str(random_position) + " " + str(random_velocity)

  #This next bit does evenly space generation
  #d = 2*np.cbrt(2*np.pi*(z_max-z_min)*(args.radius*1e-3)**2/args.number_of_macroparticles)
  #points = uniform_cylinder(d,2*args.radius*1e-3,z_max-z_min,z_min)
  #for point in points:
  #  vxmin = point.x*args.chirpx - dvx*np.sqrt(3)
  #  vx = stats.uniform.rvs(size=1,loc=vxmin,scale=2*np.sqrt(3)*dvx)

  #  vymin = point.y*args.chirpx - dvx*np.sqrt(3)
  #  vy = stats.uniform.rvs(size=1,loc=vymin,scale=2*np.sqrt(3)*dvx)

  #  vzmin = args.vz + point.z*args.chirp - dvz*np.sqrt(3)
  #  vz = stats.uniform.rvs(size=1,loc=vzmin,scale=dvz*np.sqrt(3)*2)
  #  random_velocity = Cartesian3DVector(vx,vy,vz)

  #  print str(point) + " " + str(random_velocity)

