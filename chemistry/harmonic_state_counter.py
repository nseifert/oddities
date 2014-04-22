"""
HARMONIC (VIBRATIONAL) STATE COUNTER 
FOR POLYATOMIC MOLECULES.

Implements the state counting algorithm
from:

M. J. H. Kemper, J. M. F. van Dijk, H. M. Buck, 
Chem. Phys. Lett., 53 (1), 121.


Note:
State counter doesn't count for ground state,
so calculations of the partition function require
the addition of +1 to get the right number. 
"""

__author__ = 'Nathan Seifert'
__license__ = "MIT"
__version__ = "1.0"


# Import statements
import numpy as np
import scipy.stats as spstats
import scipy.misc as misc
from matplotlib import pyplot as pp

# Constant definitions
h = 6.62606957e-34
c = 299792458
R = 8.3144621
kB = 1.3806488e-23
kT_cm = 0.695039 # Boltzmann's constant in cm-1




def state_count(E_min, E_max, modes):

	"""
	state_count(E_min, E_max, modes):

	Given an array or list of harmonic mode frequencies (modes) and 
	energy bounds E_min and E_max (in units of the harmonic mode frequencies),
	state_count() will calculate the number of possible harmonic vibrational
	states and return an array of all possible states with their harmonic 
	energies. The returned array has shape (M,N+1) where M is the number 
	of available modes between the energy bounds and N is the number of 
	harmonic modes available in the target molecule. The extra column, the last 
	column of the returned arary, contains the energies of each possible mode. 
	"""

	min_mode = np.amin(modes)
	num_columns = int(np.ceil(E_max/min_mode))
	num_modes = len(modes)
	#print "NUMBER OF COLUMNS:", num_columns

	# Initialize E difference matrix
	E_mat = np.zeros((num_modes,num_columns))

	for i in range(num_modes):
	    for j in range(0, num_columns):
	            E_mat[i,j] = modes[i]



	i = 0 # This will be our internal counter for adding quanta. 

	total_states = 0 # This will be our state counter, not used except for debug

	back_counter = 1 # This corresponds to "j" in step (2) of the pseudocode from the paper.

	quanta = np.zeros(num_modes) # The quanta-counting array
	states = [] # Results of each quanta combination with E_min < total energy < E_max will be placed here.

	while back_counter >= 0: # Program terminates when "j" is zero or less, so this is our primary cutoff
	    E_tot = 0

	    for m in range(0, num_modes): # Step (1) from paper
	        if quanta[m] > 0:
	            for n in range(int(quanta[m])):
	                 E_tot = E_tot + E_mat[m,n]

	    while i < num_modes:
	        while E_tot <= E_max:
	            quanta[i] = quanta[i] + 1
	            E_tot = E_tot + E_mat[i,quanta[i]-1]

	            while E_tot < E_min: # "Sum is too low"
	                quanta[i] = quanta[i] + 1
	                E_tot = E_tot + E_mat[i,quanta[i]-1]

	            if E_tot <= E_max:
	                total_states = total_states + 1

	                temp = []
	                for m in range(0,num_modes):
	                    temp.append(quanta[m])
	                temp.append(E_tot)
	                states.append(temp)
	            #print temp

	        E_tot = E_tot - E_mat[i,quanta[i]-1]
	        quanta[i] = quanta[i] - 1
	        i = i + 1

	    back_counter = i - 2


	    while back_counter >= 0 and quanta[back_counter] == 0:
	        back_counter = back_counter - 1

	    if back_counter >= 0:

	        quanta[back_counter] = quanta[back_counter] - 1
	        for m in range(back_counter+1,num_modes):
	            quanta[m] = 0

	        i = back_counter + 1

	return np.array(states)


def calc_q_vs_beta(T_min, T_max, N, energy_vals, direct=1, harmonic_modes = None):
	"""
	A general helper function, takes in a range of temperature
	values, the complete set of energies for the available states,
	and returns a 2D array of ln (q_vib) vs beta (1/kT) values.

	Other arguments:
	N = number of temperature steps to calculate
	direct = if energy_vals is a set of normal mode state energies,
	direct = 1. Else direct = 0 if you want to use harmonic mode 
	frequencies for your q_vib calculation. If direct = 0, 
	then harmonic_modes must be supplied (energy_vals is not used).

	One caveat: State counting is not done in this function, so 
	energy_vals must contain states with an appropriate range of
	thermal availability for the Boltzmann distribution. (e.g.
	energy_vals contains energies well larger than kT_max)
	"""


	temp_range = np.linspace(T_min,T_max,N)
	beta_range = np.zeros(N)
	for i in range(0,N):
	    beta_range[i] = (temp_range[i]*kT_cm)**(-1.0)

	q_vs_beta = np.zeros((N,2))
	for i in range(0, N):
	    if direct == 1:
	        part_func = qvib_direct(energy_vals,temp_range[i])
	    if direct == 0 and harmonic_modes:
	        part_func = qvib_est(vib_temps,temp_range[i])
	    q_vs_beta[i,1] = np.log(part_func)
	    q_vs_beta[i,0] = beta_range[i]
	return q_vs_beta


def qvib_direct(energy_vals,temp):
	"""
	Returns the partition function from an
	input array energy_vals and temperature temp.
	Energy_vals is expected to be an array of energies
	derived from all accessable harmonic states.
	"""
	return np.sum(np.exp(-1.0*energy_vals/(kT_cm*temp)))+1


def qvib_est(modes,temp):
	"""
	Returns the estimated partition function
	from an input list of harmonic mode frequencies. 
	"""
	vib_temps = modes/kT_cm
	vib_part_funcs = (1-np.exp(-vib_temps/temp))**(-1)
	return reduce(lambda x, y: x*y, vib_part_funcs)



if __name__ == '__main__':

	""" 
	Below is some test code for NH3 
	to illustrate the exploitation of the above 
	functions. 
	"""

	# Normal mode frequencies for NH3
	NH3_modes = [950.0,1628.0,1628.0,3414.0,3414.0,3337.0]


	T = 1000.0 # Temperature of experiment
	kT_temp = T * kT_cm

	# states at E_max = 5 kT
	states_5 = state_count(0.0,5.0*kT_temp,NH3_modes)
	print qvib_direct(states_5[:,-1], T)
	# states at E_max = 10 kT
	states_10 = state_count(0.0,10.0*kT_temp,NH3_modes)
	print qvib_direct(states_10[:,-1], T)

	print 'States at 5 kT, T = %d K: %s' %(T,states_5.shape)
	print 'States at 10 kT, T = %d K: %s' %(T,states_10.shape)


	# Now we'll calculate the average vibrational energy
	# of NH3 over a variety of temperatures.
	# <E> = kT^2 (d(ln q)/dT)

	T_min = 100.0 #K
	T_max = 1000.0 
	N = 100
	temp_range = np.linspace(T_min,T_max,N) # Useful for plotting 

	NH3_states = state_count(0.0, kT_cm*10.0*T_max, NH3_modes)
	NH3_q_v_b = calc_q_vs_beta(T_min, T_max, 100, NH3_states[:,-1])

	pp.plot(temp_range,NH3_q_v_b[:,1])
	pp.show()

	# Now we'll calculate the average energy
	dx = temp_range[1] - temp_range[0]
	avg_energy = np.diff(NH3_q_v_b[:,1])/dx 

	for i in range(0,N-1):
		avg_energy[i] = avg_energy[i] * kT_cm * temp_range[i]**2

	# Now we will calculate C_v = d<E>/dT
	avg_cv = np.diff(avg_energy)/dx

	# Now we can compare this to the C_p data available from NIST Webbook
	# We have to subtract 3kT to remove translation and rotation
	# and another factor of kT for Cp <--> Cv conversion

	def nh3_cp_expt(T):
	    A = 19.99563 * 0.083593 # Conversion from J/molK --> 1/cmK
	    B = 49.77119 * 0.083593
	    C = -15.37599 * 0.083593
	    D = 1.921168 * 0.083593
	    E = 0.189174 * 0.083593
	    return A + B*(T/1000.0) + C*(T/1000.0)**2 + D*(T/1000.0)**3 + E/((T/1000.0)**2) - 4*kT_cm

	cv_expt = nh3_cp_expt(temp_range)
	pp.plot(temp_range[:N-2],avg_cv,'r') # Plots calculated C_v
	pp.plot(temp_range,cv_expt,'b') # Plots expt. fit C_v data
	pp.show()


	# Plots the difference between expt and calculated C_V
	differences = [avg_cv[i] - cv_expt[i] for i in range(0,len(avg_cv))]
	pp.plot(temp_range[:N-2],differences)
	pp.show()

