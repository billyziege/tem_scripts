import argparse
from glob import glob
try:
  import cPickle as pickle
except:
  import pickle

if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Combines the input file(s) into the output file.')
  parser.add_argument('coordinate_files', nargs='+', type=str, help='Paths (supports *) to the input files to be combined.')
  parser.add_argument('-o','--output_file', dest="output_file", type=str, help='The path to where the output will be written.  Defualt is IniConditions.txt in the current working directory.', default="IniConditions.txt")
  args = parser.parse_args()

  all_files = [f for files in args.coordinate_files for f in glob(files)]

  data = []
  for filename in all_files:
    data.extend(pickle.load(open(filename,'r')))

  with open(args.output_file,'w') as f:
    pickle.dump(data,f,protocol=pickle.HIGHEST_PROTOCOL)
