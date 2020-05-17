"""
FEniCS tutorial demo program: Convection-diffusion-reaction for a system
describing the concentration of three species A, B, C undergoing a simple
first-order reaction A + B --> C with first-order decay of C. The velocity
is given by the flow field w from the demo navier_stokes_cylinder.py.

  u_1' + w . nabla(u_1) - div(eps*grad(u_1)) = f_1 - K*u_1*u_2
  u_2' + w . nabla(u_2) - div(eps*grad(u_2)) = f_2 - K*u_1*u_2
  u_3' + w . nabla(u_3) - div(eps*grad(u_3)) = f_3 + K*u_1*u_2 - K*u_3

"""

from __future__ import print_function
from fenics import *
# import dolfin
from mshr import *
from vtkplotter.dolfin import datadir, plot

# from dolfin.cpp.mesh import MeshFunctionSizet

T = 5.0  # final time
num_steps = 500  # number of time steps
dt = T / num_steps  # time step size
eps = 0.01  # diffusion coefficient
K = 10.0  # reaction rate

# Create mesh
rectangle1 = Rectangle(Point(0, 0), Point(2.5e-4, 2e-3))
rectangle2 = Rectangle(Point(2.5e-4, 0), Point(2.5e-4 + 1e-4, 2e-3))
rectangle3 = Rectangle(Point(3.5e-4, 0), Point(3.5e-4 + 2.5e-4, 2e-3))
# domain = rectangle1+rectangle2+rectangle3
# rectangle_arr = []
# for i in range(2):
#     for j in range(2):
#         rectangle_arr.append(Rectangle(Point(-1e-4+j*7e-4, 0+i*1.5e-3), Point(0+j*7e-4, 5e-4+i*1.5e-3)))
#         domain += Rectangle(Point(-1e-4+j*7e-4, 0+i*1.5e-3), Point(0+j*7e-4, 5e-4+i*1.5e-3))
rectangle_arr = [
    [Rectangle(Point(-1e-4 + j * 7e-4, 0 + i * 1.5e-3), Point(0 + j * 7e-4, 5e-4 + i * 1.5e-3)) for j in range(2)] for i
    in range(2)]
union1 = rectangle1 + rectangle_arr[0][0] + rectangle_arr[1][0]
union2 = rectangle3 + rectangle_arr[0][1] + rectangle_arr[1][1]
domain = union1 + rectangle2 + union2
domain.set_subdomain(1, union1)
domain.set_subdomain(2, rectangle2)
domain.set_subdomain(3, union2)
# mesh = generate_mesh(domain, 64)
mesh = generate_mesh(union1, 64)
domain_markers = MeshFunction('size_t', mesh, 2, mesh.domains())
# plot(mesh)

# Define function space for velocity
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
# u_D = Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)
#
# def boundary(x, on_boundary):
#     return on_boundary
#
# bc = DirichletBC(V, u_D, boundary)

u_electric_gnd = Expression('0.0', degree=2)


def boundary_electric_gnd(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol) and 5e-4 < x[1] < 15e-4


def boundary_TPB_L(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 2.5e-4, tol)


def boundary_TPB_R(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 3.5e-4, tol)


bc_electric_gnd = DirichletBC(V, u_electric_gnd, boundary_electric_gnd)
u_electric_potential = Expression('0.7', degree=2)


def boundary_electric_potential(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 6e-4, tol) and 5e-4 < x[1] < 15e-4


bc_electric_potential = DirichletBC(V, u_electric_potential, boundary_electric_potential)

# bc = [bc_electric_gnd, bc_electric_potential]
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)


class BoundaryElectricGnd(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_electric_gnd(x, on_boundary)


class BoundaryTPBL(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_TPB_L(x, on_boundary)


class BoundaryInsulationL(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and not (boundary_electric_gnd(x, on_boundary) or boundary_TPB_L(x, on_boundary))


bx0 = BoundaryElectricGnd()
bx0.mark(boundary_markers, 0)
bx1 = BoundaryTPBL()
bx1.mark(boundary_markers, 1)
bx2 = BoundaryInsulationL()
bx2.mark(boundary_markers, 2)
g = Expression('-4e3', degree=1)

boundary_conditions = {0: {'Dirichlet': u_electric_gnd},
                      # 1: {'Dirichlet': u_electric_potential},
                       1: {'Neumann': g},
                       #  2: {'Dirichlet': u_electric_potential}}
                       2: {'Neumann': 0}}

# bcs_an = []
# for i in boundary_conditions:
#     if 'Dirichlet' in boundary_conditions[i]:
#         bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'],
#                          boundary_markers, i)
#         bcs_an.append(bc)

# bc_tpb_l = DirichletBC(V, u_electric_potential, boundary_TPB_L)
# bcs_an = [bc_electric_gnd, bc_tpb_l]

bcs_an = [bc_electric_gnd]

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
dx = Measure('dx', domain=mesh, subdomain_data=domain_markers)

integrals_N = []
for i in boundary_conditions:
    if 'Neumann' in boundary_conditions[i]:
        if boundary_conditions[i]['Neumann'] != 0:
            g = boundary_conditions[i]['Neumann']
            integrals_N.append(g * v * ds(i))

f = Constant(0.0)
a = dot(grad(u), grad(v)) * dx

L = f * v * dx - sum(integrals_N)
# L = f * v * dx - g * v * ds(1)
# L = f * v * dx


# Compute solution
u = Function(V)
solve(a == L, u, bcs_an)

# Plot solution and mesh
# plot(u)
# plot(mesh)

# Save solution to file in VTK format
vtkfile = File('main/solution.pvd')
vtkfile << u

# Compute error in L2 norm
# error_L2 = errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)
# import numpy as np
# error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))
#
# # Print errors
# print('error_L2  =', error_L2)
# print('error_max =', error_max)
print(u(2.5e-4,1e-3), u(2.5e-4,2e-3))
plot(u)
