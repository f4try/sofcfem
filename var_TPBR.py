# from sympy import coth, Symbol, latex, sqrt, exp, expand, simplify, N, nsimplify
from math import sqrt, exp

from params import *


# eta_c = phis - phil - E_eq_c  # "Cathodic overvoltage"
from var_TPBL import coth

import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from sympy import exp, sqrt, coth, Symbol, latex, simplify, solve
# from math import exp, sqrt

from params import *
import numpy as np

def eta_c(eta):
    return eta -E_eq_c
xO2 = xO2_in   ##########################
p = p_c_in###################
cO2_agg = p * xO2 / KO2  # "Henry's law oxygen agglomerate concentration"

# lda_c = sqrt(i0_c * S * R_agg ^ 2 * exp(-F_const * eta_c / (2 * R_const * T)) / (
#             4 * F_const * cO2_ref * D_agg))  # "Cathodic current density subexpression"
def lda_c(eta):
    return sqrt(i0_c * S * R_agg ** 2 * exp(-F_const * eta_c(eta) / (2 * R_const * T)) / (
            4 * F_const * cO2_ref * D_agg))
# i_c = -2 * K * (1 - lda_c * coth(lda_c)) * cO2_agg  # "Cathode current density"
# e_const = 2.71828
# def coth0(x):
#     return (exp(x)+exp(-x))/(exp(x)-exp(-x))
def i_c(eta):
    return -2 * K * (1 - lda_c(eta) * coth(lda_c(eta))) * cO2_agg
    # return -2 * K * (1 - lda_c(phis3, phil) * 1) * cO2_agg

if __name__ == '__main__':
    # print(i_c(0.2,2e-3))
    # print(lda_c(0.2,2e-3))
    # print(eta_c(0.2,2e-3))

#
# phis1 = 0.2
# # phil = 2e-3
# # phis3 = 0.5
# #
# # print(eta_c(phis3, phil))
# # print(cO2_agg)
# # print(lda_c(phis3, phil))
# # print(i_c(phis3, phil))
#
#
    # phis1 = Symbol('phis1')
    # phil = Symbol('phil')
    # phis3 = Symbol('phis3')
    # fi_c = i_c(phis3, phil)
    # i_c1 = fi_c.evalf(subs={phis3:0.5,phil:2e-3})
    # print(fi_c)
    # print(i_c1)
    # ff=nsimplify(fi_c,tolerance=0.01)
    # print(ff)

    eta = Symbol('eta')
    ioc = Symbol('ioc')
    fi_a = i_c(eta)
    i_a1 = fi_a.evalf(subs={eta:0.2})
    print(latex(fi_a))
    print(i_a1)
    ff=simplify(fi_a)
    print(ff)
    x = np.linspace(-1,1,1000)
    y = [i_c(xi) for xi in x]
    plt.plot(x,y)
    plt.show()