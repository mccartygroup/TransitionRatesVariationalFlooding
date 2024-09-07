import numpy as np
import argparse 

# calculate ECDF
def calculate_ecdf(times):
    counts = np.sort(times)
    norm = times.size
    ecdf = np.arange(1, norm + 1) / norm
    ecdf_data=np.column_stack((counts, ecdf))
    return ecdf_data

parser = argparse.ArgumentParser(description="Plot ecdf of some data")

parser.add_argument('--file','-f', required=True,help='input data file')
parser.add_argument('--out','-o', required=True,help='name of output file')
args = parser.parse_args()

fpt_filename = args.file 
out_filename = args.out

print("fpt file:",fpt_filename) 
fpt_file = open(fpt_filename,'r')

times = np.loadtxt(fpt_file,usecols=(0,))

ecdf_data = calculate_ecdf(times)

np.savetxt(out_filename, ecdf_data)
