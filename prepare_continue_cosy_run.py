import sys
import os
import shutil
import csv
import argparse
from template import fill_template
from handle_simulation_parameters import ParameterSetList
from cosy_output import CosyOutput
from math import floor

def make_continuing_initial_condition_file(cosy_output_dir,step_to_continue_from,output_file):
  """
  Identifies the files in the cosy_output_dir that correspond to the
  step to continue from and concats these files into the output file.
  """
  file_list = []
  for file in os.listdir(cosy_output_dir):
    if os.path.isfile(os.path.join(cosy_output_dir,file)) and file.startswith(str(step_to_continue_from) + "-x-"):
      file_list.append(os.path.join(cosy_output_dir,file))
  with open(output_file,'w') as out_file_handle:
    for file in file_list:
      with open(file,'r') as in_file_handle:
        for line in in_file_handle:
          out_file_handle.write(line)
  if os.path.getsize(output_file) == 0:
    raise Exception("The initial conditions file is empty: " + output_file)

def edit_PE_setup_file(infile,outfile,parameter_set):
  """
  Uses the original PE setup file to write the new PE
  setup file.
  """
  new_regimes = process_regimes(infile,parameter_set.step_to_continue_from)
  with open(outfile,'w') as out_file_handle:
    count = 0 #line count for new file which is different from original
    with open(infile,'r') as in_file_handle:
      for line in in_file_handle:
        if count == 0:
          out_file_handle.write('1\n')
          out_file_handle.write(str(parameter_set.time_to_continue_from)+'\n')
          out_file_handle.write(str(os.path.basename(parameter_set.initial_conditions_file))+'\n')
          out_file_handle.write(str(parameter_set.total_particles_generated)+'\n')
        elif count < 12 and count > 2:
          out_file_handle.write(line)
        count += 1
    out_file_handle.write(str(len(new_regimes))+'\n')
    for new_regime in new_regimes:
      out_file_handle.write(new_regime["dt"]+'\n')
      out_file_handle.write(str(new_regime["steps"])+'\n')
  return
        
        

def process_regimes(infile,steps_done):
  """
  Reads the regimes (dt and steps) from the original PE setup file and
  calculates what the new ones should be from the step_to_continue_from.
  """
  new_regimes = []
  with open(infile) as f:
    lines = f.readlines()
    emmission_steps = int(lines[11])
    if int(steps_done) < emmission_steps:
      raise Exception("The run with the PE-setup.txt file, " + infile + ", needs to be rerun because the emmision process did not finish.")
    steps_done = int(steps_done) - emmission_steps
    number_regimes = int(lines[12])
    for n in range(number_regimes):
      dt_for_regime = lines[13 + 2*n].rstrip()
      steps_for_regime = int(lines[14 + 2*n])
      if steps_for_regime < int(steps_done):
        steps_done = int(steps_done) - steps_for_regime
        continue
      steps_for_regime = steps_for_regime - int(steps_done)
      new_regime = {"dt": dt_for_regime,
                    "steps": steps_for_regime}
      new_regimes.append(new_regime)
      steps_done = 0
  return new_regimes
      
    


parser = argparse.ArgumentParser(description='Builds the directory structure and the continue cosy run directory according to the provided parameter file.  Assumes that the regime where the electrons have been emmitted is finished and that at least the next regime has begun.')
parser.add_argument('filepath', type=str, help='The csv file with the specifications for the simulations desired to be continued.')
parser.add_argument('-w','--working-directory', dest="working_dir", type=str, help='The path where the output directories should originate from.',default="/mnt/home/zerbe/working_dir")
parser.add_argument('-r','--run_template', dest="run_template", type=str, help='The template with the qsub script.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/run_120ps.template")
parser.add_argument('-f','--positive_hole_files_dir', dest="pos_dir", type=str, help='The directory where the files are stored describing the positive hole geometry.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/Files")
parser.add_argument('-b','--bin_file', dest="bin_file", type=str, help='The path to the cosy binary.',default="/mnt/home/zerbe/src/tem_simulations/2-Cosy/cosy.bin")
args = parser.parse_args()

parameter_set_list = ParameterSetList()
parameter_set_list.injectParametersFromCsv(args.filepath)
for parameter_set in parameter_set_list.list:
  parameter_set.createMainOutputPath(args.working_dir)
  parameter_set.createStepnamePath("Cosy")

  #Idenfity step number to begin from.
  screen_file = os.path.join(parameter_set.returnStepnamePath("Cosy"),"screen.txt")
  cosy_output = CosyOutput()
  cosy_output.injectScreenFile(screen_file)
  last_step = cosy_output.returnTableValue("step number")
  parameter_set.step_to_continue_from = int((floor(float(last_step)/float(10))-1)*10)
  
  #Create output dir
  continue_dir = os.path.join(parameter_set.returnStepnamePath("Cosy"),"continue_from_" + str(parameter_set.step_to_continue_from))
  if not os.path.exists(continue_dir) or not os.path.isdir(continue_dir):
    os.mkdir(continue_dir)

  #Grab the correct step's files for new initial conditions file.
  parameter_set.initial_conditions_file = os.path.join(continue_dir,"InitCondition_from_" + str(parameter_set.step_to_continue_from) + ".txt")
  make_continuing_initial_condition_file(parameter_set.returnStepnamePath("Cosy"),parameter_set.step_to_continue_from, parameter_set.initial_conditions_file)

  #Edit existing PE-setup.txt file.
  relevant_row_numbers = cosy_output.returnRowNumbersWithKeyValuePair("step number",str(parameter_set.step_to_continue_from))
  parameter_set.time_to_continue_from = cosy_output.returnTableValue("time",row_number=relevant_row_numbers[0])
  parameter_set.total_particles_generated = cosy_output.returnTableValue("number macroparticles",row_number=relevant_row_numbers[0])
  pe_setup_in = os.path.join(parameter_set.returnStepnamePath("Cosy"),"PE-setup.txt")
  pe_setup_out = os.path.join(continue_dir,"PE-setup.txt")
  edit_PE_setup_file(pe_setup_in,pe_setup_out,parameter_set)

  #Link the rest of the files.
  pulse_file = os.path.join(args.pos_dir,parameter_set.returnPEMainFilename())
  hole_file = os.path.join(args.pos_dir,parameter_set.returnHoleFilename())
  initial_conditions_file = os.path.join(parameter_set.returnStepnamePath("Generate_initial_conditions"),parameter_set.returnInitialConditionsFilename())
  os.symlink(pulse_file, os.path.join(continue_dir,parameter_set.pulse_file))
  os.symlink(hole_file, os.path.join(continue_dir,parameter_set.hole_file))
  os.symlink(args.bin_file, os.path.join(continue_dir,"cosy.bin"))
  with open(os.path.join(continue_dir,"cosy_run.sh"),'w') as output:
    template_string = fill_template(args.run_template,parameter_set)
    output.write(template_string)
