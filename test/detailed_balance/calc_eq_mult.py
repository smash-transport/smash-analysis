import numpy as np
from scipy.special import kn
from scipy.integrate import quad

def Afunc(m, g, T, V):
    hbarc = 0.19732  # GeV*fm
    two_pi = 2.*np.pi
    return 2. * two_pi * g * (m**2) * T / (two_pi * hbarc)**3 * kn(2,  m/T) * V

def AGinteg(x, m, gam, g, T, V):
    return Afunc(m + 0.5 * gam * np.tan(x), g, T, V)

def AGfunc(m, Mth, Gam, g, T, V):
    a = np.arctan(2*(m - Mth)/Gam)
    BW_aver = quad(AGinteg, -a, np.pi/2, args = (m, Gam, g, T, V)) / (np.pi/2 + a)
    return BW_aver[0]  # [1] is error

T = 0.330   # GeV
V = 10.**3  # fm**3

mpi = 0.138
mrho = 0.776
gamrho = 0.149

Npi_init = 300.0
Nrho_init = 0.0

gpi = 3.
grho = 9.

N0 = Npi_init + 2* Nrho_init

Api  = Afunc(mpi, gpi, T, V)
#Arho = Afunc(mrho, grho, T,V)
Arho = AGfunc(mrho, 2*mpi, gamrho, grho, T, V)

l = ( - Api + np.sqrt( Api**2 +8 * Arho * N0) ) / (4. * Arho)

print "Should be zero: ", 2* Arho*l*l + Api*l - N0
print "Arho, Arho with width: ", Afunc(mrho, grho, T, V), AGfunc(mrho, 2*mpi, gamrho, grho, T, V)
print "one pion sort multiplicity: ", Api * l / 3.
print "one rho sort multiplicity: ", Arho * l**2 /3.

