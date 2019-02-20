### The same basic algorithm of optical.py, but using Allen's formulation.
### Also implement interpolation and don't use command line arguments for input files.

import numpy as np
np.set_printoptions(threshold=np.nan)
import sys, os, re, shutil
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as col
import matplotlib.cm as cm
from scipy.interpolate import interp1d
import matplotlib.gridspec as gridspec

# Read about gridspec
#https://matplotlib.org/users/gridspec
# One important think to nore is that the DOS, realpart of self en and the imag part all are on the same energy grid
# This saves a lot of headache later.
# We need to use matplotlib.gridspec to put multiple plots with control over the spacing between them
# But we also need an inset figure.
# Read the following page for using fig.add_axes along with subplot. I am assuming if it works with subplot it
# will work with gridspec
# https://stackoverflow.com/questions/43326680/what-are-the-differences-between-add-axes-and-add-subplot
fig=plt.figure(figsize=(8, 8), dpi= 300)
gs = gridspec.GridSpec(2, 1)
gs.update(hspace=0.1)
ax1 = plt.subplot(gs[0,0])
ax2 = plt.subplot(gs[1,0], sharex=ax1)  # https://stackoverflow.com/questions/22511550/gridspec-with-shared-axes-in-python
plt.setp(ax1.get_xticklabels(), visible=False)
#f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
#ax1 = fig.add_subplot(111)
left1, bottom1, width1, height1 = [0.65, 0.7, 0.2, 0.2]
left2, bottom2, width2, height2 = [0.65, 0.3, 0.2, 0.2]
ax3 = fig.add_axes([left1, bottom1, width1, height1])
ax4 = fig.add_axes([left2, bottom2, width2, height2])
mpl.rc('font', family='sans-serif')
mpl.rc('font', serif='Helvetica Neue')
mpl.rc('text', usetex='false')

handles, labels = ax1.get_legend_handles_labels()

def interpolate(input_file, grid):
    data = []
    comments = []
    fd = open(input_file, 'r')
    for line in fd:
        if line.startswith('#'):
            comments.append(line)
        else:
            data.append([float(item) for item in line.strip().split()])
    fd.close()
    input_data = np.array(data)  # convert the list into a numpy array
    # interpolate the input data array
    x = input_data[:,0]
    y = input_data[:,1]

    # this is the central line of the code
    f = interp1d(x, y)
    #Get the minimum and maximum values of the input domain
    Emin = np.amin(x)
    Emax = np.amax(x)

    # The following is to cut the interpolation grid from goinng beyond the input domain which is not permissible
    # If the interpolation grid is already within the input domain it won't change it.
    grid_new = []
    for item in grid:
        if (item > Emin) and (item < Emax):
            grid_new.append(item)

    xnew = np.array(grid_new)
    array_fxnew = f(xnew)   # get the interpolated array
    return array_fxnew


def FD0(E, mu, T, k):
    return 1/(np.exp((E-mu)/(k*T)) + 1)

# List the location of input files by the path and the filenames
base_path="/Users/soham/Projects/dft/H3S/qe/epw_test_mode_resolved/report"
filename1 = "realpart.elself.interpolated.denser.dat"
filename2 = "linewidth.elself.interpolated.zero_near_zero.denser.dat"
filename3 = "Yundi_Kuan_H3S.dos1ev"
# Join the base path to filename to create the path for opening the files
path_to_realpart = os.path.join(base_path, filename1)
path_to_linewidth = os.path.join(base_path, filename2)
path_to_dos = os.path.join(base_path, filename3)

#####################################################################################################################
## Constants and grid spacing here
#####################################################################################################################
# Choose a spacing for w' grid. This is for the self-energy
# Note: Choose a spacing and endpoints such that w'+w also maps to the w'_array for every w.
# Otherwise one will have to interpolate the self-energy on the w'+w grid for every value of w. That will be unwise.
delta_omegap = 0.0001
# Choose a omega grid - this is the grid for the photon energy
delta_omega = 0.005
list_tau_inv = [0.025, 0.050, 0.075]  # noninteracting inverse scattering time in eV
num_plots = len(list_tau_inv)
array_tau_inv = np.array(list_tau_inv)
v = 0.25*(10**6)  # the Fermi velocity
e = 1.6*10**(-19)  # Charge of electron
Ef = 0.63  # DOS at Fermi level in states/eV
a = 5.6*5.292*(10**(-11))
epsilon = 8.854*10**(-12)  # Permittivity of free space
munot = 4*np.pi*10**(-7)  # vacuum permeability
hbar = 6.58*(10**(-13))  # reduced Plank's constant
T = 200 #K
# Boltzmann's constant
k = 8.6173303*10**(-5)  # in eV.K-1
mu = 0.0  # chemical potential
N0 = Ef/(e*0.5*a**3)
conduct_coeff = N0*(0.33*v**2)*(e**2)*hbar
w_p = np.sqrt(N0*(0.33*v**2)*(e**2)/epsilon)*(hbar/1000)  # plasma frequency in eV

#######################################################################################################################
colormap = plt.cm.plasma
plt.gca().set_prop_cycle(plt.cycler('color', colormap(np.linspace(0, 0.85, num_plots))))

print ("Conductance coefficient:", conduct_coeff, "Plasma Frequency:", w_p)

# Define your own input arrays - can be different (though still within the bounds) than those in the files
# This gives us the flexibility to check for intergation convergence
# We accomplish this by interpolating, soon to follow

#Note: As mentioned above, choose a spacing and endpoints such that w'+w also maps to the w'_array for every w.
# Otherwise one will have to interpolate the self-energy on the w'+w grid for every value of w. that will be unwise.
array_omegap = np.arange(-5, 3, delta_omegap) # The energy grid for the self-energy.
array_omega_a = np.arange(0.002, 0.02, delta_omega) # The energy grid for the incoming infrared photon
array_omega_b = np.arange(0.02, 0.80, delta_omega) # The energy grid for the incoming infrared photon
array_omega = np.concatenate((array_omega_a, array_omega_b))
array_M = interpolate(path_to_realpart, array_omegap)/1000.0
array_Gamma = interpolate(path_to_linewidth, array_omegap)/1000.0
# The following scheme shifts the energy grid so that M = 0 at omegap = 0
asign = np.sign(array_M)  # Get the signs of the M array
signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)  # check where the sign changes - those are the places where M crosses the omega1 = 0 line
signchange[0] = 0  # neglect the first location
sign_index = np.where(signchange==1)  # find the indices where the sign changes
shift_omegap = np.amin(np.abs(array_omegap[sign_index]))  # pick the smallest absolute value of omega1 where the sign changes - thats the origin
array_omegap = array_omegap - shift_omegap   # Shift the E array


# Some gymnastics so that we do not go beyond the bounds and try to read a nonexistant value for the self-energy at w'+w
# We truncate the upper energy end of the energy array without changing the self-energy array
# So that we have a valid self-energy for w'+w.
cutoff = array_omegap[-1] - array_omega[-1] # so that cutoff + array_omega[-1] gives a value where the self-energy exists
array_omegap_full = array_omegap
cutoff_location = np.searchsorted(array_omegap, cutoff) - 1
array_omegap = array_omegap[:cutoff_location]  #  So that w'+w does not give an energy where the self energy is not defined.




# Now we do the 3-dimensional integral
# we shall create a 2-D numpy array each for each of the Spectral functions A(E, omega1) and G(E, omega1 - omega)
# I use "w" and "omega" interchangably in the comments

for tau_inv in list_tau_inv:
    list_conduct = []
    list_conduct_drude = []
    # we now create a 1-D complex array of self energy from the real and imaginary part
    array_selfen = array_M - 1j*(array_Gamma + tau_inv) # Note the negative sign in th imaginary part - done to make the imaginary part always negative
    array_selfen_drude = 0.0*array_M  - 1j*np.full_like(array_Gamma, tau_inv) # zero real-part, constant imag part
    # find out the indices of the array array_omegap_full that correspond to array_omegap values
    array_omegap_indices = np.where(np.isin(array_omegap_full, array_omegap))
    # Note that the new omegap grid is truncated on the upper end - so we need to create a self-energy array for that truncated grid
    # we use numpy indexing: https://stackoverflow.com/questions/25201438/python-how-to-get-values-of-an-array-at-certain-index-positions
    array_selfen_omegap = array_selfen[array_omegap_indices]
    array_selfen_drude_omegap = array_selfen_drude[array_omegap_indices]
    for omega in array_omega: # the w array for the infrared photon
        # create the omegap+omega array:
        array_omegap_omega = array_omegap + omega
        # find out the indices of the array array_omegap_full that correspond to array_omegap_omega values
        # The following rounding of digits is needed because we shall be checking if array_wp_w is in array_omegap_full and they have different sigfigs (idk why)
        array_omegap_omega_indices = np.where(np.isin(np.around(array_omegap_full, decimals=7), np.around(array_omegap_omega, decimals=7)))
        # Create the self-energy array  on the  omegap+omega array
        array_selfen_omegap_omega = array_selfen[array_omegap_omega_indices]
        array_selfen_drude_omegap_omega = array_selfen_drude[array_omegap_omega_indices]
        # Here we calculate the 1-D array in eqn 10 of Allen 2015 PRB.
        Integrand = (1/omega)*(FD0(array_omegap, mu, T, k) - FD0(array_omegap_omega, mu, T, k))/(
                omega - array_selfen_omegap_omega + np.conj(array_selfen_omegap))
        drudeIntegrand = (1/omega)*(FD0(array_omegap, mu, T, k) - FD0(array_omegap_omega, mu, T, k))/(
                omega - array_selfen_drude_omegap_omega + np.conj(array_selfen_drude_omegap))

        conduct = (1j)*conduct_coeff*np.trapz(Integrand, x=array_omegap)
        conduct_drude = (1j)*conduct_coeff*np.trapz(drudeIntegrand, x=array_omegap)

        list_conduct.append(conduct)
        list_conduct_drude.append(conduct_drude)

    # Convert conductivity lists into numpy arrays
    array_conduct = np.array(list_conduct)
    array_conduct_drude = np.array(list_conduct_drude)
    # Get epsilon and Reflectivity
    array_epsilon = 1+ (1j)*array_conduct/(array_omega)
    array_epsilon_drude = 1+ (1j)*array_conduct_drude/(array_omega)
    array_N = np.sqrt(munot*array_epsilon)
    array_N_drude = np.sqrt(munot*array_epsilon_drude)
    array_R = (np.absolute((1-array_N)/(1+array_N)))**2
    array_R_drude = (np.absolute((1-array_N_drude)/(1+array_N_drude)))**2

    tau_inv = tau_inv*1000  # convert in meV
    print (" Min, max reflectivity:", np.amin(array_R), np.amax(array_R))

    #title = "optical_conductivity_Allen"
    #title = "optical_conductivity_Drude"
    #title = "optical_conductivity_Allen_Drude"
    #title = "Epsilon_Allen"
    title = "Reflectivity_Allen_temp"

    ax1.tick_params(axis='x', labelsize=16)
    ax1.tick_params(axis='y', labelsize=16)
    ax1.set_xlim(0,800)
    #ax1.set_ylim(0,100000)
    #ax1.set_xlabel("$\omega$ (meV)",fontsize=20)
    ax1.set_xlabel("$\omega$ (meV)",fontsize=20)
    ax1.set_ylabel("$4\\pi\sigma^{'}/\omega_{p}^{2}, 4\\pi\sigma^{'}/\omega_{p}^{2}(10^{8}\\times S/m)$",fontsize=20)
    ax2.set_ylabel("$4\\pi\sigma^{''}/\omega_{p}^{2}, 4\\pi\sigma^{''}/\omega_{p}^{2}(10^{8}\\times S/m)$",fontsize=20)
    ax3.set_xlabel("$\omega$ (meV)",fontsize=14)
    ax3.set_ylabel("$(\sigma^{'}_{Drude} - \sigma^{'})/\omega_{p}^{2} (10^{8}\\times S/m)$",fontsize=14)
    ax4.set_xlabel("$\omega$ (meV)",fontsize=14)
    ax4.set_ylabel("$(\sigma^{''}_{Drude} - \sigma^{''})/\omega_{p}^{2} (10^{8}\\times S/m)$",fontsize=14)
    #ax1.set_ylabel("$\sigma$",fontsize=20)
    #ax1.set_ylabel("$\epsilon^{''}$",fontsize=20)
    #ax1.set_ylabel("Reflectivity",fontsize=20)
    #ax1.plot(array_omega, np.real(array_epsilon),linewidth=2,label="Real epsilon")
    #ax1.plot(array_omega, np.imag(array_epsilon),linewidth=2,label="Imag epsilon")
    ax1.plot(array_omega*1000, 4*np.pi*np.real(array_conduct)/w_p**2,linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    ax1.scatter(array_omega*1000, 4*np.pi*np.real(array_conduct_drude)/w_p**2,s=7)#, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    ax2.plot(array_omega*1000, 4*np.pi*np.imag(array_conduct)/w_p**2,linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    ax2.scatter(array_omega*1000, 4*np.pi*np.imag(array_conduct_drude)/w_p**2,s=7)#, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    ax3.plot(array_omega*1000, 4*np.pi*(np.real(array_conduct_drude) - np.real(array_conduct))/w_p**2,linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    ax4.plot(array_omega*1000, 4*np.pi*(np.imag(array_conduct_drude) - np.imag(array_conduct))/w_p**2,linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    #ax1.scatter(array_omega*1000, np.imag(array_conduct)/(10**8),s=7)
    #ax1.plot(array_omega*1000, np.real(array_conduct_drude)/(10**8),linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    #ax1.scatter(array_omega*1000, np.imag(array_conduct_drude)/(10**8),s=7)
    #ax1.plot(array_omega*1000, -(np.real(array_conduct)-np.real(array_conduct_drude))/(10**8),linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    #ax1.scatter(array_omega*1000, -(np.imag(array_conduct) -np.imag(array_conduct_drude))/(10**8),s=7, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    #ax1.plot((1/tau_inv)*array_omega, np.real(array_conduct_drude)/(10**8),linewidth=3, label="$\\tau^{-1}$ = %4.3f eV" %tau_inv)
    #ax1.scatter((1/tau_inv)*array_omega, np.imag(array_conduct_drude)/(10**8),s=7)
    #ax1.scatter(array_omega, np.imag(array_conduct),s=7, color="red",label="Imag, eq 10")
    #ax1.plot(array_omega*1000, np.real(array_N),linewidth=3, label="$\\tau^{-1}$ = %4.1f meV" %tau_inv)
    #ax1.scatter(array_omega*1000, np.imag(array_N),s=7)
    #ax1.plot(array_omega*1000, array_R,linewidth=3, label="$\\tau^{-1}$ = %4.1f meV" %tau_inv)
    #ax1.plot(array_omega*1000, array_R_drude,linewidth=3, label="$\\ Drude tau^{-1}$ = %4.1f meV" %tau_inv)
    #ax1.text(-1.4,-0.0,'$\lambda=$%3.1f,$\Omega=$%3.1f,T=$\Omega$/%2.0f,\n NDOS=%d,\n N$\omega$=%d.'%(lambda_val,Omega,OmegabyT,NDOS,Nomega),fontsize=30)
    handles, labels = ax1.get_legend_handles_labels()
    lgd = plt.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.75, 0.75), borderaxespad=0.)
    ax1.legend(fontsize=14)
    ax2.legend(fontsize=14)
    #ax3.legend(fontsize=14)
    #ax4.legend(fontsize=14)
#plt.show()
plt.savefig(title+".png", dpi=600, bbox_inches='tight')
plt.close()

