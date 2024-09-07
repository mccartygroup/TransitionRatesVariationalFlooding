from math import *
import re

# try importing filterfalse tools
try:
    # Python 2
    from itertools import ifilterfalse as filterfalse
except ImportError:
    # Python 3
    from itertools import filterfalse

########## Edit the following definitions  ################
###########################################################
Temperature=260   # Temperature in Kelvins
dtmd=0.2           # timestep in COLVAR file (fs)
CV1col=2            # Column of comb1 CV (index starting from 1)
vesbiascol=4       # Column of VES bias (index starting from 1)
fname='COLVAR'     # name of COLVAR file
fout='COLVAR-RW'   # name of output COLVAR

################################################################
CV1col=CV1col-1              # python index starts from 0
vesbiascol=vesbiascol-1    # python index starts from 0

# for skipping headers
def iscomment(s):
    return s.startswith('#')

# Conversion from energy in units of kBT
kT=(2.479/298.0)*Temperature
beta=1.0/kT
print('kBT is',kT)
print('beta is',beta)

conversion=1e-6   # conversion factor ps ---> us

f=open(fname,'r')
fout=open(fout,'w')

timesum=0.0
previous_cv1 = None
for line in filterfalse(iscomment, f):
    line=line.strip()
    columns=line.split()
    cv1_str = columns[CV1col]
    bias_str = columns[vesbiascol]
    
    if len(columns) >= CV1col+1 and len(columns) >= vesbiascol+1:
        # Check if the columns are not empty
        if re.match(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$', cv1_str) and re.match(r'^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$', bias_str):
            cv1 = float(cv1_str)
            bias = float(bias_str)
            
            if previous_cv1 is not None:
                exp_arg = beta * previous_bias
                if exp_arg > 700:  # Approximate threshold to avoid overflow
                    timesum += dtmd
                else:
                    timesum = timesum + exp(exp_arg) * dtmd
                fout.write(str(timesum * conversion) + ' ' + str(previous_cv1) + '\n')
            
            previous_cv1 = cv1
            previous_bias = bias
        else:
            print(f"Skipping line: {line}, incomplete or non-numeric columns")
    else:
        print(f"Skipping line: {line}, incomplete line")

if previous_cv1 is not None:
    exp_arg = beta * previous_bias
    if exp_arg > 700:  # Approximate threshold to avoid overflow
        timesum += dtmd
    else:
        timesum = timesum + exp(exp_arg) * dtmd
    fout.write(str(timesum * conversion) + ' ' + str(previous_cv1) + '\n')

print('Unbiased first passage time (us)', timesum * conversion)

f.close()
fout.close()
