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
# from vtkplotter import Plotter, screenshot, exportWindow
from vtkplotter.dolfin import datadir, plot
from params import *
from var_TPBL import *
from var_TPBR import *
import numpy as np

# from dolfin.cpp.mesh import MeshFunctionSizet


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
mesh = generate_mesh(domain, 64)
# mesh1 = generate_mesh(union1, 64)
# mesh2 = generate_mesh(rectangle2, 64)
# mesh3 = generate_mesh(union2, 64)
# Get mesh parts
# domain_markers = MeshFunction('size_t', mesh, 2, mesh.domains())
cell_f = MeshFunction('size_t', mesh, 2, mesh.domains())
cell_f.set_all(0)
for cell in cells(mesh):
    if cell.midpoint().x() < 2.5e-4:
        cell_f[cell] = 1
    elif 2.5e-4 <= cell.midpoint().x() <= 3.5e-4:
        cell_f[cell] = 2
    else:
        cell_f[cell] = 3
mesh1 = SubMesh(mesh, cell_f, 1)
mesh2 = SubMesh(mesh, cell_f, 2)
mesh3 = SubMesh(mesh, cell_f, 3)
# plot(mesh1)
# Define function space for velocity
W = VectorFunctionSpace(mesh, 'P', 2)
# V = FunctionSpace(mesh, 'P', 1)
V1 = FunctionSpace(mesh1, 'P', 1)
V2 = FunctionSpace(mesh2, 'P', 1)
V3 = FunctionSpace(mesh3, 'P', 1)
# P1 = FiniteElement('P', triangle, 1)
# element = MixedElement([P1, P1, P1])
# V = FunctionSpace(mesh, element)
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
    # on_boundary = True
    tol = 1E-14
    return on_boundary and near(x[0], 2.5e-4, tol)


def boundary_TPB_R(x, on_boundary):
    # on_boundary = True
    tol = 1E-14
    return on_boundary and near(x[0], 3.5e-4, tol)


# bc_electric_gnd = DirichletBC(V.sub(0), u_electric_gnd, boundary_electric_gnd)
bc_electric_gnd = DirichletBC(V1, u_electric_gnd, boundary_electric_gnd)

u_electric_potential = Expression('V_cell', V_cell=V_cell, degree=2)


def boundary_electric_potential(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 6e-4, tol) and 5e-4 < x[1] < 15e-4


# bc_electric_potential = DirichletBC(V.sub(2), u_electric_potential, boundary_electric_potential)
bc_electric_potential = DirichletBC(V3, u_electric_potential, boundary_electric_potential)

# bc = [bc_electric_gnd, bc_electric_potential]
# boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)
boundary_markers = MeshFunction("size_t", mesh, 1, 0)
# boundary_markers.set_all(0)

boundary_markers1 = MeshFunction("size_t", mesh1, mesh1.topology().dim() - 1, 0)
boundary_markers2 = MeshFunction("size_t", mesh2, mesh2.topology().dim() - 1, 0)
boundary_markers3 = MeshFunction("size_t", mesh3, mesh3.topology().dim() - 1, 0)


class BoundaryElectricGnd(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_electric_gnd(x, on_boundary)


class BoundaryTPBL(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_TPB_L(x, True)


class BoundaryTPBR(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_TPB_R(x, True)


class BoundaryElectricPontential(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_electric_potential(x, on_boundary)


class BoundaryInsulationL(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and not (boundary_electric_gnd(x, on_boundary) or boundary_TPB_L(x, on_boundary)
                                    or boundary_TPB_R(x, on_boundary) or boundary_electric_potential(x, on_boundary))


# bx1 = BoundaryElectricGnd()
# bx1.mark(boundary_markers, 1)
# bx2 = BoundaryTPBL()
# bx2.mark(boundary_markers, 2)
# bx3 = BoundaryTPBR()
# bx3.mark(boundary_markers, 3)
# bx4 = BoundaryElectricPontential()
# bx4.mark(boundary_markers, 4)
# bx5 = BoundaryInsulationL()
# bx5.mark(boundary_markers, 5)

bx0 = BoundaryElectricGnd()
bx0.mark(boundary_markers1, 0)
bx1 = BoundaryTPBL()
bx1.mark(boundary_markers1, 1)
bx1.mark(boundary_markers2, 1)
bx2 = BoundaryTPBR()
bx2.mark(boundary_markers2, 2)
bx2.mark(boundary_markers3, 2)
bx3 = BoundaryElectricPontential()
bx3.mark(boundary_markers3, 3)
bx4 = BoundaryInsulationL()
bx4.mark(boundary_markers1, 4)
bx4.mark(boundary_markers2, 4)
bx4.mark(boundary_markers3, 4)

phis1 = TrialFunction(V1)
phis3 = TrialFunction(V3)
phil = TrialFunction(V2)
v1 = TestFunction(V1)
v2 = TestFunction(V2)
v3 = TestFunction(V3)
# phix = [phis1, phil, phis3]


# v1, v2, v3 = TestFunctions(V)
# phix = Function(V)
# phis1, phil, phis3 = split(phix)

# g_L = Expression(("-2800.0", "0."),degree=2)
# g_R = Expression(("2800.0", "0."),degree=2)
# g_L = Constant(-2800.0)
# g_R = Constant(2800.0)
# g_L = i_a(phis1, phil)
# g_R = i_c(phis3, phil)
for loop in range(2):
    if loop == 0:
        g_L = Constant(-2800.0)
        g_R = Constant(2800.0)

    else:
        class IAExpression(UserExpression):
            def eval(self, value, x):
                # value[0] = -5600*random()
                # x = [0.00025, 0.0014375]
                a = phis1F(2.5e-4, x[1])
                b = philF(2.5e-4, x[1])
                c = i_a(a, b)
                print(x, type(x))
                print(a, type(a))
                print(b, type(b))
                print(c, type(c))
                value[0] = c


        class ICExpression(UserExpression):
            def eval(self, value, x):
                # x = [0.00035, 0.0014375]
                # value[0] = 5600 * random()
                print(x, type(x))
                value[0] = i_c(phis3F(3.5e-4, x[1]), philF(3.5e-4, x[1]))


        g_L = IAExpression()
        g_R = ICExpression()

    boundary_conditions = {0: {'Dirichlet': u_electric_gnd},
                           1: {'Neumann': g_L},
                           2: {'Neumann': g_R},
                           3: {'Dirichlet': u_electric_potential},
                           4: {'Neumann': 0}}
    # bcs = []
    # for i in boundary_conditions:
    #     if 'Dirichlet' in boundary_conditions[i]:
    #         bc = DirichletBC(V, boundary_conditions[i]['Dirichlet'],
    #                          boundary_markers, i)
    #         bcs.append(bc)

    # bc_tpb_l = DirichletBC(V, u_electric_potential, boundary_TPB_L)
    # # bcs_an = [bc_electric_gnd, bc_tpb_l]
    # bc_TBPL1 = DirichletBC(V.sub(0), 0.3, boundary_TPB_L)
    # bc_TBPL2 = DirichletBC(V.sub(1), 0.4, boundary_TPB_L)
    # bc_TBPR1 = DirichletBC(V.sub(1), 0.5, boundary_TPB_R)
    # bc_TBPR2 = DirichletBC(V.sub(2), 0.6, boundary_TPB_R)
    # bcs = [bc_electric_gnd, bc_electric_potential,bc_TBPL1,bc_TBPL2,bc_TBPR1,bc_TBPR2]
    # bcs = [bc_electric_gnd, bc_electric_potential]
    bc1 = [bc_electric_gnd]
    bc2 = []
    bc3 = [bc_electric_potential]
    # Define variational problem

    # ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)
    # dS = Measure('dS', domain=mesh, subdomain_data=boundary_markers)
    # dx = Measure('dx', domain=mesh, subdomain_data=cell_f)

    ds1 = Measure('ds', domain=mesh1, subdomain_data=boundary_markers1)
    ds2 = Measure('ds', domain=mesh2, subdomain_data=boundary_markers2)
    ds3 = Measure('ds', domain=mesh3, subdomain_data=boundary_markers3)

    # dx = Measure("dx")[cell_f]
    # ds = Measure("ds")[boundary_markers]
    # dS = Measure("dS")[boundary_markers]

    integrals_N1 = []
    integrals_N2 = []
    integrals_N3 = []
    for i in boundary_conditions:
        if 'Neumann' in boundary_conditions[i]:
            if boundary_conditions[i]['Neumann'] != 0:
                g = boundary_conditions[i]['Neumann']
                if i == 1:
                    integrals_N1.append(g * v1 * ds1(i))
                    integrals_N2.append(-g * v2 * ds2(i))
                elif i == 2:
                    integrals_N2.append(-g * v2 * ds2(i))
                    integrals_N3.append(g * v3 * ds3(i))

    integrals_IN1 = []
    integrals_IN2 = []
    integrals_IN3 = []
    for i in boundary_conditions:
        if 'Internal Neumann' in boundary_conditions[i]:
            if boundary_conditions[i]['Internal Neumann'] != 0:
                g = boundary_conditions[i]['Internal Neumann']
                if i == 1:
                    integrals_IN1.append(-g * v1 * ds(i))
                    integrals_IN2.append(g * v2 * ds(i))
                elif i == 2:
                    integrals_IN2.append(g * v2 * ds(i))
                    integrals_IN3.append(-g * v3 * ds(i))

    # # Define function G such that G \cdot n = g
    # class BoundarySource(UserExpression):
    #     def __init__(self, mesh, **kwargs):
    #         self.mesh = mesh
    #         super().__init__(**kwargs)
    #     def eval_cell(self, values, x, ufc_cell):
    #         cell = Cell(self.mesh, ufc_cell.index)
    #         n = cell.normal(ufc_cell.local_facet)
    #         g = sin(4e3*x[1])
    #         values[0] = g*n[0]
    #         # values[1] = g*n[1]
    #     # def value_shape(self):
    #     #     return (1,)
    #
    # G1 = BoundarySource(mesh1, degree=2)
    # G2_L = BoundarySource(mesh2, degree=2)
    # G2_R = BoundarySource(mesh2, degree=2)
    # G3 = BoundarySource(mesh3, degree=2)
    #
    # bc1.append(DirichletBC(V1, G1, boundary_TPB_L))
    # bc2.append(DirichletBC(V2, G2_L, boundary_TPB_L))
    # bc2.append(DirichletBC(V2, G2_R, boundary_TPB_R))
    # bc3.append(DirichletBC(V3, G3, boundary_TPB_R))

    f = Constant(0.0)
    n = FacetNormal(mesh)
    # a = 10 * dot(grad(u), grad(v1)) * dx(0)+ dot(grad(u), grad(v1)) * dx(2)

    # a1 = 1. * inner(grad(phis1), grad(v1)) * dx(0)\
    #      +1.2 * inner(grad(phis1), grad(v1)) * dx(1)\
    #      +1.3 * inner(grad(phis1), grad(v1)) * dx(2)\
    #      +1.4 * inner(grad(phis1), grad(v1)) * dx(3)
    a1 = kappa_s * dot(grad(phis1), grad(v1)) * dx
    L1 = f * v1 * dx - sum(integrals_N1)
    # L1 = f * v1 * dx(0)+f * v1 * dx(1)+f * v1 * dx(2)\
    #      +f * v1 * dx(3)
    # L1 = - sum(integrals_N) - sum(integrals_IN)
    # L1 = f * v1 * dx - g_L("+")*v1("+")*ds(3)
    # L1 = f * v1 * dx - g_L*v1*ds1(1)
    # L1 = f * v1 * dx
    # L = f * v1 * dx - 2800 * v1 * ds(1) + 2800 * v1 * ds(2)
    a2 = kappa_m * dot(grad(phil), grad(v2)) * dx
    L2 = f * v2 * dx - sum(integrals_N2)
    # L2 = f * v2 * dx
    a3 = kappa_s * dot(grad(phis3), grad(v3)) * dx
    L3 = f * v3 * dx - sum(integrals_N3)
    # L3 = f * v3 * dx

    # Compute solution
    phis1FS = Function(V1)
    philFS = Function(V2)
    phis3FS = Function(V3)

    solve(a1 == L1, phis1FS, bc1)
    solve(a2 == L2, philFS, bc2)
    solve(a3 == L3, phis3FS, bc3)

    phis1F = phis1FS
    philF = philFS
    phis3F = phis3FS
    # F = a1 + a2 + a3 - L1 - L2 - L3
    # F =
    # linear_solver.parameters [“ maximum_iterations”] = 500
    # solve(F == 0, phix, bcs)
    # F = assemble(F)

    # # Separate left and right hand sides of equation
    # a, L = lhs(F), rhs(F)
    #
    # # Solve problem
    # solve(a == L, phix, bcs)

    # J = derivative(F, phix)
    # problem = NonlinearVariationalProblem(F, phix, bcs, J)
    # solver = NonlinearVariationalSolver(problem)
    # prm = solver.parameters
    # prm['newton_solver']['absolute_tolerance'] = 1E-8
    # prm['newton_solver']['relative_tolerance'] = 1E-7
    # prm['newton_solver']['maximum_iterations'] = 25
    # prm['newton_solver']['relaxation_parameter'] = 1.0
    # solver.solve()

    # Plot solution and mesh
    # plot(u)
    # plot(mesh)

    # Save solution to file in VTK format
    # phis1F, philF, phis3F = phix.split()
    vtkfile1 = File('main/phis1.pvd')
    vtkfile1 << phis1F
    vtkfile2 = File('main/phil.pvd')
    vtkfile2 << philF,
    vtkfile3 = File('main/phis3.pvd')
    vtkfile3 << phis3F

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
    print(phis1F(2.5e-4, 1e-3), phis1F(2.5e-4, 2e-3))
    print(philF(2.5e-4, 1e-3), philF(2.5e-4, 2e-3))
    print(philF(3.5e-4, 1e-3), philF(3.5e-4, 2e-3))
    print(phis3F(3.5e-4, 1e-3), phis3F(3.5e-4, 2e-3))

# plt1 = plot(phis1F,vmin=0.,vmax=0.7,interactive=False, scalarbar=False)
# # plt2 = plot(phil, add=True, scalarbar='horizontal')
# plt2 = plot(philF, add=True)
# plt3 = plot(phis3F,vmin=0.,vmax=0.7, add=True, scalarbar=False)
# plot()

plt1 = plot(phis1F, interactive=False, scalarbar=False)
# # plt2 = plot(phil, add=True, scalarbar='horizontal')
plt2 = plot(philF, add=True, scalarbar=False)
plt3 = plot(phis3F, add=True)
plot()
