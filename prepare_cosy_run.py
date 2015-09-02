import sys
import os
import shutil
import csv
import argparse
from template import fill_template
from handle_simulation_parameters import ParameterSetList

parser = argparse.ArgumentParser(description='Builds the directory structure and the cosy run directory according to the provided parameter file. Assumes the the initial conditions have already been generated.')
parser.add_argument('filepath', type=str, help='The csv file with the specifications for the simulations desired to be run.')
parser.add_argument('-w','--working-directory', dest="working_dir", type=str, help='The path where the output directories should originate from.',default="/mnt/home/zerbe/working_dir")
parser.add_argument('-r','--run_template', dest="run_template", type=str, help='The template with the qsub script.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/run_120ps.template")
parser.add_argument('-s','--setup_template', dest="setup_template", type=str, help='The template of the setup parameters.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/PE-setup_120ps.template")
parser.add_argument('-f','--positive_hole_files_dir', dest="pos_dir", type=str, help='The directory where the files are stored describing the positive hole geometry.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/Files")
parser.add_argument('-b','--bin_file', dest="bin_file", type=str, help='The path to the cosy binary.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/cosy.bin")
args = parser.parse_args()

parameter_set_list = ParameterSetList()
parameter_set_list.injectParametersFromCsv(args.filepath)
for parameter_set in parameter_set_list.list:
  parameter_set.createMainOutputPath(args.working_dir)
  parameter_set.createStepnamePath("Cosy")
  initial_conditions_file = os.path.join(parameter_set.returnStepnamePath("Generate_initial_conditions"),parameter_set.returnInitialConditionsFilename())
  pulse_file = os.path.join(args.pos_dir,parameter_set.returnPEMainFilename())
  hole_file = os.path.join(args.pos_dir,parameter_set.returnHoleFilename())
  parameter_set.extractTotalParticlesFromIniConditionFile(initial_conditions_file)
  os.symlink(initial_conditions_file, os.path.join(parameter_set.returnStepnamePath("Cosy"),parameter_set.initial_conditions_file))
  os.symlink(pulse_file, os.path.join(parameter_set.returnStepnamePath("Cosy"),parameter_set.pulse_file))
  os.symlink(hole_file, os.path.join(parameter_set.returnStepnamePath("Cosy"),parameter_set.hole_file))
  os.symlink(args.bin_file, os.path.join(parameter_set.returnStepnamePath("Cosy"),"cosy.bin"))
  with open(os.path.join(parameter_set.returnStepnamePath("Cosy"),"PE-setup.txt"),'w') as output:
    template_string = fill_template(args.setup_template,parameter_set)
    output.write(template_string)
  with open(os.path.join(parameter_set.returnStepnamePath("Cosy"),"cosy_run.sh"),'w') as output:
    template_string = fill_template(args.run_template,parameter_set)
    output.write(template_string)
