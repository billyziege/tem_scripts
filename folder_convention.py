import sys
import os
import argparse

def parse_folder_for_metadata(path):
  """
  My naming convention stores metadata in the folder.  This extracts it.
  """
  output_dict = {}
  folder = os.path.basename(path)
  pieces = folder.split("_")
  output_dict["shape"] = pieces[0]
  output_dict["n_particles"], output_dict["e_per_particle"] = pieces[1].split("x")
  output_dict["total_electrons"] = str(int(output_dict["n_particles"])*int(output_dict["e_per_particle"]))
  output_dict["applied_field"] = pieces[2].replace("MV-per-m","")
  output_dict["steps"] = pieces[3]
  return output_dict



if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Test functions in this package')
  parser.add_argument('-f','--folder', dest="folder", type=str, help='The path(s) to the folder needing to be parsed.', default=None)
  args = parser.parse_args()

  if args.folder is not None:
    print str(parse_folder_for_metadata(args.folder))
