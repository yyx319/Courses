from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
    

#Simplistic Poisson noise - assume that in each time bin, we have a probability
#of detection of 1%
prob_det = 0.01
n_samples = 1000
n_trials = 1000

#Create our sample array.
samples = np.zeros( (n_samples, n_trials) )
n_photons = n_samples*prob_det
plot_max = int(n_photons + 5*np.sqrt(n_photons))
x = np.arange(plot_max, dtype='int')

random_photons = (np.random.random( size=samples.shape ) < prob_det).astype(int)
#Just count up the photons in 100 bins (10 seconds)
manual_poisson = np.sum(random_photons, 0)
plt.clf()
plt.hist(manual_poisson, label='Sample-by-sample', bins=x-.5)
plt.xlabel('Number of Photons')
plt.ylabel('Frequency')
plt.title('Click to continue...')
plt.legend()
plt.pause(.001)
plt.ginput(1)

#Compare this to the numpy distribution
numpy_poisson = np.random.poisson(np.ones(n_trials)*n_photons)
plt.hist(numpy_poisson, label='np.random.poisson', bins=x-.5)

#Finally... compare to the analytic function for a poisson distribution
x_fact = [np.math.factorial(xx) for xx in x]
plt.plot(x, n_trials*n_photons**x*np.exp(-n_photons)/x_fact, label='Analytic Function')
plt.title('Poisson distribution for n={0:5.1f}'.format(n_photons))
plt.legend()
