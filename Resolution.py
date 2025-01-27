import numpy as np

Me = 0.51099895069
NUM_BINS = 2500 # Resolution of plots
E_MAX = 4.5 # MeV
E = np.linspace(0, E_MAX, NUM_BINS)
dE = E_MAX / NUM_BINS
