import numpy as np
from scipy.optimize import curve_fit
import argparse

def model(x,p1,p2):
    b = 0.01
    return 1 - np.exp( -(pow(1.0 + b*x, p1) -1.0) * p2)
    #return 1.0 - np.exp(A * (1 - np.exp(B * x )))

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
#x_data_mean = np.mean(x_data)
#print(x_data_mean)

#x_data_norm = x_data / x_data_mean

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


#print('beta*gamma*a',B)

#B = B/x_data_mean 

some_data = np.arange(np.min(x_data),np.max(x_data),0.1)
tcdf = model(some_data,fit_params[0],fit_params[1])
fout = open(out_filename,'w')

for i in range(len(some_data)):
        fout.write(str(some_data[i])+' '+str(tcdf[i])+'\n')

fout.close()


