import sys
import os
import csv
import argparse
import subprocess
import select
from template import fill_template
from handle_simulation_parameters import ParameterSetList

parser = argparse.ArgumentParser(description='Builds the directory structure and the initial input directory according to the provided parameter file.')
parser.add_argument('filepath', type=str, help='The csv file with the specifications for the simulations desired to be run.')
parser.add_argument('-w','--working-directory', dest="working_dir", type=str, help='The path where the output directories should originate from.',default="/mnt/home/zerbe/working_dir")
parser.add_argument('-t','--template_file', dest="template_file", type=str, help='The input template for the generating initial conditions script.',default="/mnt/home/zerbe/src/tem_simulations/1-Generate_init_file/input.template")
parser.add_argument('-g','--generate_script', dest="generate_script", type=str, help='The path for generating initial conditions script.',default="/mnt/home/zerbe/src/tem_simulations/1-Generate_init_file/a.out")
parser.add_argument('-p','--probability_script', dest="probability_script", type=str, help='The path for generating the electrons script.',default="/mnt/home/zerbe/src/tem_simulations/1-Generate_init_file/generate_probability.R")
parser.add_argument('-d','--delta_t', dest="delta_t", type=str, help='The duration of the laser pule in fs.  Default is 50 fs.',default=50)
args = parser.parse_args()

parameter_set_list = ParameterSetList()
parameter_set_list.injectParametersFromCsv(args.filepath,lazer_duration=args.delta_t)
for parameter_set in parameter_set_list.list:
  parameter_set.createMainOutputPath(args.working_dir)
  parameter_set.createStepnamePath("Generate_initial_conditions")
  command =[]
  command.append(args.probability_script)
  command.append(str(parameter_set.timesteps))
  command.append(os.path.join(parameter_set.returnStepnamePath("Generate_initial_conditions"),"probability.txt"))
  subprocess.call(command)
  with open(os.path.join(parameter_set.returnStepnamePath("Generate_initial_conditions"),"input.txt"),'w') as output:
    template_string = fill_template(args.template_file,parameter_set)
    output.write(template_string)
  os.chdir(parameter_set.returnStepnamePath("Generate_initial_conditions"))
  command2 = []
  command2.append(args.generate_script)
  subprocess.call(command2)
