import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy import exp, sqrt, coth, Symbol, latex, simplify, solve, nonlinsolve, dsolve
# from math import exp, sqrt

from params import *
import numpy as np

xH2 = xH2_in  #################
p = p_a_in  #############

cH2_agg = p * xH2 / KH2  # "Henry's law hydrogen agglomerate concentration"


# # eta_a = phis1 - phil - E_eq_a  # "Anodic overpotential"
# def eta_a(x, y, phis1, phil):
#     return phis1(x, y) - phil(x, y) - E_eq_a
#
#
# # beta_a = cH2_agg * (1 - exp(-2 * F_const * eta_a / (R_const * T)))  # ""
# def beta_a(x, y, phis1, phil):
#     return cH2_agg * (1 - exp(-2 * F_const * eta_a(x, y, phis1, phil) / (R_const * T)))
#
#
# lda_a = sqrt(i0_a * S * R_agg ** 2 / (2 * F_const * cH2_ref * D_agg))  # "Anodic current density subexpression"
#
#
# # i_a = K * (1 - lda_a * coth(lda_a)) * beta_a  # "Anode current density"
# def i_a(x, y, phis1, phil):
#     return K * (1 - lda_a * coth(lda_a)) * beta_a(x, y, phis1, phil)  # "Anode current density"


# eta_a = phis1 - phil - E_eq_a  # "Anodic overpotential"
def eta_a(eta):
    return eta - E_eq_a


# beta_a = cH2_agg * (1 - exp(-2 * F_const * eta_a / (R_const * T)))  # ""
def beta_a(eta):
    return cH2_agg * (1 - exp(-2 * F_const * eta_a(eta) / (R_const * T)))


lda_a = sqrt(i0_a * S * R_agg ** 2 / (2 * F_const * cH2_ref * D_agg))  # "Anodic current density subexpression"


# i_a = K * (1 - lda_a * coth(lda_a)) * beta_a  # "Anode current density"
def coth(x):
    return (exp(x)+exp(-x))/(exp(x)-exp(-x))
def i_a(eta):
    return K * (1 - lda_a * coth(lda_a)) * beta_a(eta)  # "Anode current density"
    # return K * (1 - lda_a * 1) * beta_a(phis1, phil)  # "Anode current density"

if __name__ == '__main__':
    # print(i_a(0.5,2e-3))

    # phis1 = Symbol('phis1')
    # phil = Symbol('phil')
    # phis3 = Symbol('phis3')
    eta = Symbol('eta')
    ioc = Symbol('ioc')
    fi_a = i_a(eta)
    i_a1 = fi_a.evalf(subs={eta:0.2})
    print(latex(fi_a))
    print(i_a1)
    ff=simplify(fi_a)
    print(ff)
    x = np.linspace(-1,1,1000)
    y = [i_a(xi) for xi in x]
    # plt.plot(x,y)
    # plt.show()
    # print(solve(2-2*exp(-6*eta)-ioc,eta))
    print(nonlinsolve(218940.044783997 - 218940.044783997*exp(-65.7479785482249*eta) - ioc, eta))
