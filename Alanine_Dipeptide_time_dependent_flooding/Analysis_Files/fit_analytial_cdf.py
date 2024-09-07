import numpy as np
from scipy.optimize import curve_fit
import argparse

def model(x,A,B):
    return 1.0 - np.exp(A * (1 - np.exp(B * x )))

parser = argparse.ArgumentParser(description="Plot ecdf of some data")

parser.add_argument('--file','-f', required=True,help='input data file')
parser.add_argument('--out','-o', required=True,help='name of output file')
args = parser.parse_args()

ecdf_filename = args.file
out_filename = args.out

print("ecdf:",ecdf_filename)
ecdf_file = open(ecdf_filename,'r')

data = np.loadtxt(ecdf_file)

x_data = data[:,0]
y_data = data[:,1]
x_data_mean = np.mean(x_data)
print(x_data_mean)

x_data_norm = x_data / x_data_mean

initial_guess = (1,1)
fit_params, _ = curve_fit(model,x_data_norm,y_data,p0=initial_guess)

print(fit_params)

A = fit_params[0]
B_norm = fit_params[1]

B = B_norm / x_data_mean
k = A*B
print('k',k)
print('beta*gamma*a',B)

tcdf = model(x_data,A,B)
fout = open(out_filename,'w')

for i in range(len(x_data)):
        fout.write(str(x_data[i])+' '+str(tcdf[i])+'\n')

tau = 1.0/k
print('tau (ps)',tau)
print('tau (us)',tau/10**6)

beta = input("what is beta?")
avalue = input("what is a?")

timestep = 0.002

gamma = (B/float(beta)/float(avalue))*timestep
print('gamma',gamma)

fout.close()
