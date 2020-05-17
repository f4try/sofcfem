#!/usr/bin/env python
"""
First solve the stationary electric conduction problem. Then use its
results to solve the evolutionary heat conduction problem.

Run this example as on a command line::

    $ python <path_to_this_file>/thermal_electric.py
"""
from __future__ import absolute_import
import sys
import numpy as nm

import params
from var_TPBR import i_c

sys.path.append( '.' )
import os
from var_TPBL import i_a
import sfepy.discrete.fem.periodic as per


filename_mesh = 'sofc2d_mesh.mesh'
# Time stepping for the heat conduction problem.
t0 = 0.0
t1 = 0.5
n_step = 11

# Material parameters.
specific_heat = 1.2

##########

cwd = os.path.split(os.path.join(os.getcwd(), __file__))[0]

options = {
    'absolute_mesh_path' : True,
    'output_dir' : os.path.join(cwd, 'output')
}

regions = {
    'Omega' : 'all',
    'Omega1' : 'cells of group 1',
    'Omega2' : 'cells of group 2',
    'Omega3' : 'cells of group 3',
    'Omega2_Surface_L': ('r.Omega1 *v r.Omega2', 'facet'),
    'Omega2_Surface_L1' : ('copy r.Omega2_Surface_L', 'facet', 'Omega1'),
    'Omega2_Surface_L2' : ('copy r.Omega2_Surface_L', 'facet', 'Omega2'),
    'Omega2_Surface_R': ('r.Omega2 *v r.Omega3', 'facet'),
    'Omega2_Surface_R1': ('copy r.Omega2_Surface_R', 'facet', 'Omega2'),
    'Omega2_Surface_R2': ('copy r.Omega2_Surface_R', 'facet', 'Omega3'),
    'Left' : ('vertices in (x < %f)&(y > %f)&(y < %f)' % (0.1e-4,0.49e-3,1.51e-3),'facet'),
    'Right' : ('vertices in (x > %f)&(y > %f)&(y < %f)' % (5.9e-4,0.49e-3,1.51e-3),'facet'),
    # 'Left' : ('vertices in (x < %f)' % (-0.9e-4),'edge'),
    # 'Right' : ('vertices in (x > %f)' % (6.9e-4),'edge'),

}

def get_is(ts, coors, mode=None, **kwargs):
    if mode == 'qp':
        nqp, dim = coors.shape
        phis1 = 0.5
        phil = 2e-3
        phis3 = 0.2
        eta1 = 0
        eta2 = 1
        a = i_a(eta1)
        c = i_c(eta2)
        out = {'i_a': a*nm.ones((nqp,1,1), dtype=nm.float64),
               'i_c': c*nm.ones((nqp, 1, 1), dtype=nm.float64),
               }
        return out

def set_electric_bc(coor):
    y = coor[:,1]
    ymin, ymax = y.min(), y.max()
    val = 2.0 * (((y - ymin) / (ymax - ymin)) - 0.5)
    return val
# def get_circle(coors, domain=None):
#     r = nm.sqrt(coors[:,0]**2.0 + coors[:,1]**2.0)
#     return nm.where(r < 0.2)[0]
def get_shift1(ts, coors, region):
    # val = 0.1 * coors[:, 0]
    val = nm.empty_like(coors[:, 1])
    val.fill(-params.E_eq_c)
    return val

def get_shift2(ts, coors, region):
    val = nm.empty_like(coors[:, 1])
    val.fill(-params.E_eq_a)

    return val
# def mat_fun(ts, coors, mode=None, **kwargs):
#     if mode == 'qp':
#         nqp, dim = coors.shape
#         alpha = nm.zeros((nqp,1,1), dtype=nm.float64)
#         alpha[0:nqp // 2,...] = alpha1
#         alpha[nqp // 2:,...] = alpha2
#         K = nm.eye(dim, dtype=nm.float64)
#         K2 = nm.tile(K, (nqp,1,1))
#         out = {
#             'K' : K2,
#             'f_1': 2800.0 * nm.ones((nqp,1,1), dtype=nm.float64),
#             'f_2': -20.0 * nm.ones((nqp,1,1), dtype=nm.float64),
#             'G_alfa': G_bar * alpha,
#
#             }
#
#         return out
functions = {
    'get_is' : (get_is,),
    'get_shift1' : (get_shift1,),
    'get_shift2' : (get_shift2,),
    'match_y_line' : (per.match_y_line,),
    'match_x_line' : (per.match_x_line,),
    # 'set_electric_bc' : (lambda ts, coor, bc, problem, **kwargs:
    #                      set_electric_bc(coor),),
}

materials = {
    'is': ('get_is'),
    'm' : ({
        'thermal_conductivity' : 2.0,
        'electric_conductivity' : 1.5,
    },),
}

# The fields use the same approximation, so a single field could be used
# instead.
fields = {
    'temperature': ('real', 1, 'Omega', 1),
    'potential_l' : ('real', 1, 'Omega', 1),
    'potential_s' : ('real', 1, 'Omega', 1),
}

variables = {
    # 'T' : ('unknown field', 'temperature', 0, 1),
    # 's' : ('test field', 'temperature', 'T'),
    'phil' : ('unknown field', 'potential_l', 1),
    'psil' : ('test field', 'potential_l', 'phil'),
    'phis': ('unknown field', 'potential_s', 2),
    'psis': ('test field', 'potential_s', 'phis'),
    # 'phi_known_l' : ('parameter field', 'potential_l', '(set-to-None)'),
    # 'phi_known' : ('parameter field', 'potential_s', '(set-to-None)'),
    # 'eta_l' : ('parameter field', 'potential_l', '(set-to-None)'),
    # 'eta_s' : ('parameter field', 'potential_s', '(set-to-None)'),
}

# ics = {
#     'ic' : ('Omega', {'T.0' : 0.0}),
# }

ebcs = {
	'left' : ('Left', {'T.0' : 0.0, 'phis.0' :0.4}), #V_cell000
    'right' : ('Right', {'T.0' : 2.0, 'phis.0' : 0.0}),
    # 'inside' : ('Omega2_Surface_L', {'phi.0' : 'set_electric_bc'}),
    # 'inside1' : ('Omega2_Surface_L', {'phi.0' : 0.5}),
    # 'inside2' : ('Omega2_Surface_R', {'phi.0' : 0.2}),
}

lcbcs = {
    'shifted1' : (('Omega2_Surface_L2', 'Omega2_Surface_L1'),
                  {'phil.all' : 'phis.all'},
                  'match_y_line', 'shifted_periodic',
                  'get_shift1'),
    'shifted2' : (('Omega2_Surface_R1', 'Omega2_Surface_R2'),
                  {'phil.all' : 'phis.all'},
                  'match_y_line', 'shifted_periodic',
                  'get_shift2'),
}
integrals = {
    'i' : 1
}

equations = {
    # '2' : """%.12e * dw_volume_dot.2.Omega( s, dT/dt )
    #          + dw_laplace.2.Omega( m.thermal_conductivity, s, T )
    #          = dw_electric_source.2.Omega( m.electric_conductivity,
    #                                        s, phi_known ) """ % specific_heat,
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_surface_integrate.i.Left(flux.val1,psi)""",
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_electric_source.i.Left(m.electric_conductivity,psi,phi)""",
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_surface_integrate.1.Omega2_Surface_L(flux.val1,psi)""",
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_bc_newton.1.Omega2_Surface_L(exchange.val1,exchange.val1,psi,phi)""",
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_surface_integrate.1.Omega2_Surface_L(flux.val1,psi)""",
    # '1' : """dw_laplace.2.Omega( m.electric_conductivity, psi, phi ) = dw_surface_integrate.1.Omega2_Surface_L(is.i_a,psi)+dw_surface_integrate.1.Omega2_Surface_R(is.i_c,psi)""",
    'eq1' : """dw_laplace.2.Omega( m.electric_conductivity, psil, phil ) = 0""",
    'eq2' : """dw_laplace.2.Omega( m.electric_conductivity, psis, phis ) = 0""",
}

# solvers = {
#     'ls' : ('ls.scipy_direct', {}),
#     'newton' : ('nls.newton', {
#         'i_max'      : 1,
#         'eps_a'      : 1e-10,
#         'problem'   : 'nonlinear',
#     }),
#     'ts' : ('ts.simple', {
#         't0'     : t0,
#         't1'     : t1,
#         'dt'     : None,
#         'n_step' : n_step, # has precedence over dt!
#         'verbose' : 1,
#     }),
# }
solvers = {
    'ls' : ('ls.scipy_direct', {}),
    'newton' : ('nls.newton', {
        'i_max'      : 1,
        'eps_a'      : 1e-10,
    }),
}
