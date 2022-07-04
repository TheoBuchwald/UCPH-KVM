import numpy as np
import matplotlib.pyplot as plt

#Copy of LevelSystem5.m by Stine

if __name__ == '__main__':
    #Coupling elements(eV)
    V01L = 1.09116
    V10L = 1.09116
    V0m1L = 1.4449
    Vm10L = 1.4449
    V12L = 1.9374
    V21L = 1.9374
    Vm1m2L = 2.5959
    Vm2m1L = 2.5959

    V01R = 1.09116
    V10R = 1.09116
    V0m1R = 1.4449
    Vm10R = 1.4449
    V12R = 1.9374
    V21R = 1.9374
    Vm1m2R = 2.5959
    Vm2m1R = 2.5959

    #Energies(eV)
    Em2 = -36325.0807329464
    Em1 = -36335.1122096212
    E0 = -36342.3605059506
    E1 = -36346.0032939133
    E2 = -36346.9788220261

    # Uncomment if your energies are in Hartree
    # Em2 /= 0.036749405469679
    # Em1 /= 0.036749405469679
    # E0 /= 0.036749405469679
    # E1 /= 0.036749405469679
    # E2 /= 0.036749405469679

    #Constants
    C = 0.1 # The energy conservation (eV^-1)
    hbar = 6.58211899*10**(-16) # (eV*s)
    T = 0.025 # Temperture in eV
    mu0 = - 5.1 # er det det ??? (eV)
    e = 1.602176565*10**(-19) # Elementary charge (C)
    Wk= (2*np.pi)/hbar*C # Combined constant

    dV = 0.1
    Vmin, Vmax = -15, 15 # Change to vary voltage range calculated
    num_V = int(np.ceil((Vmax - Vmin) / dV))

    V = np.linspace(Vmin,Vmax,num_V)   # Range of source drain voltage


    muR = mu0-0.5*V # Shifts in the electrochemical potential
    muL = mu0+0.5*V

    #Define Fermi function
    if T==0:
        fermi = lambda x: 0.5*(-x/np.abs(x) + 1)
    else:
        fermi = lambda x: 1 / (np.exp(x/T) + 1)

    #Calculations of weights for master equations
    #First to the right
    wm2m1 = (Vm2m1L)**2 * (1 - fermi(Em1 - Em2 - muL)) + (Vm2m1R)**2 * (1 - fermi(Em1 - Em2 - muR))
    wm10 = (Vm10L)**2 * (1 - fermi(E0 - Em1 - muL)) + (Vm10R)**2 * (1 - fermi(E0 - Em1 - muR))
    w01 = (V01L)**2 * (1 - fermi(E1 - E0 - muL)) + (V01R)**2 * (1 - fermi(E1 - E0 - muR))
    w12 = (V12L)**2 * (1 - fermi(E2 - E1 - muL)) + (V12R)**2 * (1 - fermi(E2 - E1 - muR))

    w21 = (V21L)**2 * fermi(E2 - E1 - muL) + (V21R)**2 * fermi(E2 - E1 - muR)
    w10 = (V10L)**2 * fermi(E1 - E0 - muL) + (V10R)**2 * fermi(E1 - E0 - muR)
    w0m1= (V0m1L)**2 * fermi(E0 - Em1 - muL) + (V0m1R)**2 * fermi(E0 - Em1 - muR)
    wm1m2 = (Vm1m2L)**2 * fermi(Em1 - Em2 - muL) + (Vm1m2R)**2 * fermi(Em1 - Em2 - muR)
    #Then to the left
    wm2m1L = (Vm2m1L)**2 * (1 - fermi(Em1 - Em2 - muL))
    wm10L = (Vm10L)**2 * (1 - fermi(E0 - Em1 - muL))
    w01L = (V01L)**2 * (1 - fermi(E1 - E0 - muL))
    w12L = (V12L)**2 * (1 - fermi(E2 - E1 - muL))

    w21L = (V21L)**2 * fermi(E2 - E1 - muL)
    w10L = (V10L)**2 * fermi(E1 - E0 - muL)
    w0m1L = (V0m1L)**2 * fermi(E0 - Em1 - muL)
    wm1m2L = (Vm1m2L)**2 * fermi(Em1 - Em2 - muL)


    #Current I from populations
    I = Wk*10**3*e*((wm1m2L) / (1 + (wm1m2/wm2m1)*(1 + (w0m1/wm10) + (w10/w01) * (w0m1/wm10) + (w21/w12) * (w10/w01) * (w0m1/wm10))) \
        + (w0m1L - wm2m1L) / (1 + (wm2m1/wm1m2) + (w0m1/wm10)*(1 + (w10/w01) + (w21/w12) *(w10/w01))) \
        + (w10L-wm10L) / (1 + (wm10/w0m1) * (1 + (wm2m1/wm1m2)) + (w10/w01) * (1 + (w21/w12))) \
        + (w21L-w01L) / (1+(w01/w10)*(1+ (wm10/w0m1) + (wm2m1/wm1m2) * (wm10/w0m1)) + w21/w12) \
        - (w12L) / (1 + (w12/w21) * (1 + (w01/w10) + (wm10/w0m1) * (w01/w10) + (wm2m1/wm1m2) * (wm10/w0m1) * (w01/w10))))

    I *= 1E9

    plt.rc('font', size=22)
    plt.rc('axes', titlesize=25)

    fig,ax = plt.subplots(figsize=(10,5))
    ax.plot(V, I)

    ax.set_ylabel('I  (nA)', size=22)
    ax.set_xlabel('$V_{SD} (V)$', size=22)
    ax.tick_params(size=6, labelsize=16)

    plt.savefig('Conducting_Magnitude_NBD_SCN_Ag.png', format='png', dpi=600)

    # plt.show()
