import numpy as np
from scipy.stats import ks_2samp
from scipy.optimize import curve_fit
from statsmodels.distributions.empirical_distribution import ECDF 
import matplotlib.pyplot as plt

f = open('rescaled-crossing-times.dat', 'r')

# define theoretical CDF
def func(x, tau):
    return 1 - np.exp(-x / tau)

x = []
count = 0
for line in f:
    line = line.strip()
    columns = line.split()
    x.append(float(columns[0]))
    count = count + 1

x = np.array(x)

# for numerical stability we divide by the mean
mu = x.mean()
x = x / mu

# now compute empirical CDF
ecdf = ECDF(x)

# plot ECDF
x1 = np.linspace(min(x), max(x))
y1 = ecdf(x1)
plt.step(x1 * mu, y1, 'k-', lw=3.)

# fit to theoretical CDF to obtain tau
popt, pcov = curve_fit(func, x1, y1)
tau = popt[0]

mean_data = mu
best_fit_tau = tau * mu

yfit = func(x1, tau)
# plot fit
plt.plot(x1 * mu, yfit, 'b-', lw=3.)

# for p-value
# now generate some random data with the same exponential distribution
np.random.seed(12345678)
x2 = np.random.exponential(1 / tau, 1000)

st, p = ks_2samp(x2, x)
p_value = p

plt.xscale('log')
plt.xlabel('time [s]')
plt.ylabel('Cumulative Probability')

ecdftcdf = "ECDF_TCDF_plot_rescaled.png"
plt.savefig(ecdftcdf)

# Save results to a file
output_filename = "rescaled-output.txt"
with open(output_filename, "w") as output_file:
    output_file.write(f"Mean of data: {mean_data}\n")
    output_file.write(f"Best fit tau: {best_fit_tau}\n")
    output_file.write(f"P-value: {p_value}\n")

print(f"Results saved to {output_filename}")
