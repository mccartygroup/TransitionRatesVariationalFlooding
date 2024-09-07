import numpy as np
from scipy.optimize import curve_fit
import argparse

def calculate_ecdf(times):
    counts = np.sort(times)
    norm = times.size
    ecdf = np.arange(1, norm + 1) / norm
    ecdf_data = np.column_stack((counts, ecdf))
    return ecdf_data

def model(x,p1,p2):
    b = 0.01
    return 1 - np.exp( -(pow(1.0 + b*x, p1) -1.0) * p2)

parser = argparse.ArgumentParser(description="Plot ecdf of some data")

parser.add_argument('--file','-f', required=True,help='input data file')
#parser.add_argument('--out','-o', required=True,help='name of output file')
args = parser.parse_args()

fpt_filename = args.file
#out_filename = args.out

print("fpt file:",fpt_filename)
fpt_file = open(fpt_filename,'r')

times = np.loadtxt(fpt_file,usecols=(0,))
n = len(times)
n_bootstraps = 50
i = n_bootstraps
kvalues = []
while i > 0:
    t = np.random.choice(times, size=n, replace=True)
    ecdf = calculate_ecdf(t)
    x_data = ecdf[:,0]
    y_data = ecdf[:,1]
    gamma0 = 0.6
    a = 12.02716
    b = 0.01

    p0 = np.zeros(2)
    p0[0] = gamma0*a + 1.0
    p0[1] = 0.0000001

    fit_params, _ = curve_fit(model,x_data,y_data,p0=p0)
    print(fit_params)

    print('gamma',(fit_params[0]-1.0)/a)
    k=fit_params[1]*fit_params[0]*b
    print('k',k)

    tau = 1.0/k
    print('tau',tau,'ps')
    tau_s = tau*1.e-12
    print('tau',tau_s,'s')
    k_o_s = 1.0/tau_s
    print('k_o',k_o_s,'s-1')

    kvalues.append(k_o_s)
    i = i - 1

kvalues = np.array(kvalues)
#print(np.mean(kvalues))
#print(np.std(kvalues))
lower_percentile = np.percentile(kvalues,30)
upper_percentile = np.percentile(kvalues,70)

print('lower 30th percentile',lower_percentile)
print('upper 70th percentile',upper_percentile)
