"""
IODINE FRANCK-CONDON FACTOR CALCULATOR
AND VIBRONIC SPECTRUM PREDICTOR

Uses a Morse potential formalism to calculate 
the Franck-Condon factors for the X-->B 
vibronic transition for iodine. 
Experimental data taken from NIST.

Requires Robert Johansson's "wavefunction"
library, found here:
https://github.com/jrjohansson/wavefunction
"""

__author__ = "Nathan Seifert"
__license__ = "MIT"
__version__ = "1.0"

# PART -1:
# Import statements
# =================
from wavefunction import *
from wavefunction.wavefunction1d import *
import numpy as np
from matplotlib import pyplot as pp

# PART 0:
# Initialization Steps
# ====================

# Constant definitions
h = 6.62606957e-34 # J s
h_ = h / (2 * np.pi) # hbar
c = 299792458.0  # m/s
mass_I = 126.90447 # amu
mu = (mass_I) ** 2 / (2 * mass_I) * 1.66053892E-27  # converts to Kg
N_a = 6.0221413E23 # mol^-1
cm_kJ_mol = 0.012  # Conversion factor from cm-1 to kJ/mol
kb = 0.695039 # cm-1
T = 298.15 # K

# Returns a Morse potential
def morse(x, args):
    re = args['re']
    b = args['beta']
    de = args['de']

    return de * (1 - np.exp(-b * (x - re))) ** 2

# Returns a Boltzmann factor
boltzmann = lambda e, t: np.exp(-e / (kb * t))

# Returns a Gaussian
gauss = lambda x, A, mu, sigma: A * np.exp(-(x - mu) ** 2 / (2. * sigma ** 2))

# Calculates the anharmonicity constant for the Morse potential
# Inputs must be in cm-1, returns constant in units of 1/Angstrom
beta = lambda de, ve: ve * np.pi * np.sqrt(2 * mu * c / (h * de)) * 1E-9


# Experimental potential parameters for I2, X (1sigmaG) state (in kJ/mol)
X_params = {'de': 18713.0 * cm_kJ_mol, 're': 2.6663, 'beta': beta(18713.0,214.502), 'te': 0.0}
print "BETA (1/A) for X state", X_params['beta']

# Experimental potential parameters for I2, B (3piU) state (in kJ/mol)
B_params = {'de': 5168.7 * cm_kJ_mol, 're': 3.0247, 'beta': beta(5168.7,125.967),'te': 15769.01}
print "BETA (1/A) for B state", B_params['beta']


# Define grid for solving Schrodinger equation
x_min = X_params['re'] - 1.0
x_max = B_params['re'] + 2.0
N = 750
x_space = np.linspace(x_min, x_max, N)
print 'GRID BOUNDS FOR PROBLEM --- X_MIN:', x_min, " X_MAX:", x_max

# Set up Morse potentials for plotting
u_X = morse(x_space, X_params)
u_B = morse(x_space, B_params)


# Part 1:
# Solving Schrodinger Equation
# =================================

# Kinetic energy prefactor for Hamiltonian
k = -h_ ** 2 / (2 * mu) * (1.0E10) ** 2 * N_a * 0.001
print "KINETIC ENERGY PREFACTOR", k, "kJ*As^2/mol"

# Helper function for solving Schrodinger equation
def solve_it(N, prefactor, func, grid, args):
    potential = assemble_u_potential(N,func,grid,args)
    ke_term = assemble_K(N,prefactor,grid[0],grid[-1])
    pot_term = assemble_V(N,potential,grid[0],grid[-1])

    return solve_eigenproblem(ke_term + pot_term) # returns tuple (evals,evecs)


# Solve Schrodinger Equation for both states
evals_X, evecs_X = solve_it(N, k, morse, x_space, X_params)
energies_X = evals_X.real * 84.0  # converts to cm-1
print "EIGENVALUES FOR FIRST FIVE VIBRATIONAL STATES OF X STATE:", energies_X[0:5]

evals_B, evecs_B = solve_it(N, k, morse, x_space, B_params)
energies_B = evals_B.real * 84.0
print "EIGENVALUES FOR FIRST FIVE VIBRATIONAL STATES OF B STATE:", energies_B[0:5]


# Plot first few wavefunctions
pp.plot(x_space, u_X)
for i in range(5):
    Y = evals_X[i] + 10 * np.real(evecs_X[i])

    wvfn_mask = np.where(Y > u_X - 20.0)
    evals_mask = np.where(Y > u_X)

    pp.plot(x_space[evals_mask], evals_X[i] * np.ones(np.shape(x_space))[evals_mask], 'k--')
    pp.plot(x_space[wvfn_mask], Y[wvfn_mask])

# Pull out vertical energy scale for B state for safe keeping (we'll need it later too)
te_B = B_params['te']
te_B_kJ = te_B * cm_kJ_mol

pp.plot(x_space,u_B+te_B_kJ)
for i in range(5):
    Y = evals_B[i] + 10 * np.real(evecs_B[i]) + te_B_kJ

    wvfn_mask = np.where(Y > u_B + te_B_kJ - 10.0)
    evals_mask = np.where(Y > u_B + te_B_kJ)

    pp.plot(x_space[evals_mask], evals_B[i] * np.ones(np.shape(x_space))[evals_mask] + te_B_kJ, 'k--')
    pp.plot(x_space[wvfn_mask], Y[wvfn_mask])

pp.xlim(2.2, 3.6)
pp.ylim(0, 300)
pp.show()



# PART 2:
# Calculate Franck-Condon Factors
# ================================


# Calculate our Boltzmann weightings and only account for X vibrational states with at least 1% population
weights = np.zeros(10)
for i in range(0, 10):  # 10 is arbitrary, but we shouldn't need this many at 298K
    weights[i] = (boltzmann(energies_X[i] - energies_X[0], T))
print 'Boltzmann weights for the first 10 vibrational states of X:', weights

# Set X & B vibrational state cutoffs for Franck-Condon factors...
X_cutoff = len(np.array(weights)[weights>0.01]) # States with proportional abundance > 1%
B_cutoff = 50 # Arbitrary
print X_cutoff

# Calculate Franck-Condon Factors
overlaps = np.zeros((B_cutoff,X_cutoff)).astype(np.complex)
for i in range(0, X_cutoff):
    for j in range(0, B_cutoff):
        ket = wavefunction_normalize(x_space, evecs_X[i])
        bra = wavefunction_normalize(x_space, evecs_B[j])

        overlaps[j, i] = inner_product(x_space, bra, ket) ** 2

# Print maximum factors
def print_max_factors(factor_matrix):
    print '\n==========='
    for i in range(np.shape(factor_matrix)[1]):
        print 'MAXIMUM FRANCK-CONDON FACTOR FOR X(v = %d) --> B(v = %d): %0.3f' \
              % (i,np.argmax(factor_matrix[:,i].real),np.amax(factor_matrix[:,i].real))
    print '===========\n'
print_max_factors(overlaps)


# Calculate transition frequencies. This array is a 1D array of 2D (freq,intensity) arrays. The first
# axis specifies the X vibrational state we're beginning at, and the second/third axes specify the frequency to a
# given B vibrational state and its intensity.
transitions = np.zeros((X_cutoff, B_cutoff, 2))
for i in range(0, X_cutoff):
    for j in range(0,B_cutoff):
        transitions[i, j, 0] = 1.E7 / ((energies_B[j] + te_B) - energies_X[i]) # converts to nm
        transitions[i, j, 1] = np.real(overlaps[j,i]) * weights[i]


# Part 3:
# Calculate predicted vibronic spectrum
# ==============================================

# Defines our frequency space we're plotting the spectrum on
NUM_POINTS = 5000
spec_x = np.linspace(np.amin(transitions[0,0,0]) - 300, np.amax(transitions[-1,0,0]) + 300, num=NUM_POINTS)

# Creates single spectra for each transition in transitions. Has the same arrangement as transitions, with
# three axes: (lower state, upper state, spectrum point)
spectrum = np.zeros((X_cutoff, B_cutoff,NUM_POINTS))
linewidth = 0.25 / 2.3548  # 0.25 nm FWHM
for i in range(X_cutoff):
    for j in range(B_cutoff):
        spectrum[i,j,:] = gauss(spec_x,transitions[i,j,1],transitions[i,j,0],linewidth)


# Combine all the spectra in spectrum into a final spectrum
final_spec = np.zeros((NUM_POINTS))
for i in range(X_cutoff):
    for j in range(B_cutoff):
        final_spec += spectrum[i, j, :]

# Plot composite spectrum
pp.xlabel(r'$wavelength \/(nm)$')
pp.ylabel(r'$Absorbance \/(arb units)$')
pp.xlim(450, 1000)
pp.title(
    r'$Vibronic\/ spectrum\/ for \/ the\/ X\/ ^1\Sigma^{+}_{g} \/ \rightarrow \/ B\/ ^3\Pi_{0+u} \/ transition \/ of \/ I_2$')
pp.plot(spec_x, final_spec, 'b')
pp.show()




