from __future__ import print_function
from dolfin import *
from mshr import *
from var_TPBL import *
from var_TPBR import *

# Create mesh
rectangle1 = Rectangle(Point(0, 0), Point(2.5e-4, 2e-3))
rectangle2 = Rectangle(Point(2.5e-4, 0), Point(2.5e-4 + 1e-4, 2e-3))
rectangle3 = Rectangle(Point(3.5e-4, 0), Point(3.5e-4 + 2.5e-4, 2e-3))
rectangle_arr = [
    [Rectangle(Point(-1e-4 + j * 7e-4, 0 + i * 1.5e-3), Point(0 + j * 7e-4, 5e-4 + i * 1.5e-3)) for j in range(2)] for i
    in range(2)]
union1 = rectangle1 + rectangle_arr[0][0] + rectangle_arr[1][0]
union2 = rectangle3 + rectangle_arr[0][1] + rectangle_arr[1][1]
domain = union1 + rectangle2 + union2
domain.set_subdomain(1, union1)
domain.set_subdomain(2, rectangle2)
domain.set_subdomain(3, union2)
mesh_full = generate_mesh(domain, 64)

# Define domains
domains_markers = MeshFunction("size_t", mesh_full, 2)
domains_markers.set_all(0)

# Create mesh for the left,mid,right domain
for cell in cells(mesh_full):
    if cell.midpoint().x() < 2.5e-4:
        domains_markers[cell] = 1
    elif 2.5e-4 <= cell.midpoint().x() <= 3.5e-4:
        domains_markers[cell] = 2
    else:
        domains_markers[cell] = 3
mesh1 = SubMesh(mesh_full, domains_markers, 1)
mesh2 = SubMesh(mesh_full, domains_markers, 2)
mesh3 = SubMesh(mesh_full, domains_markers, 3)
# plot(mesh_right)

# Define boundaries
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

def boundary_electric_potential(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 6e-4, tol) and 5e-4 < x[1] < 15e-4

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


boundary_markers1 = MeshFunction("size_t", mesh1, 1)
boundary_markers2 = MeshFunction("size_t", mesh2, 1)
boundary_markers3 = MeshFunction("size_t", mesh3, 1)
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

# Define function space for velocity
V1 = FunctionSpace(mesh1, "Lagrange", 1)
V2 = FunctionSpace(mesh2, "Lagrange", 1)
V3 = FunctionSpace(mesh3, "Lagrange", 1)
V = MixedFunctionSpace(V1, V2, V3)

# Create all of the weak form variables
phis1, phil, phis3 = TrialFunctions(V)
v1, v2, v3 = TestFunctions(V)

# Define boundary condition
u_electric_gnd = Expression('0.0', degree=2)
bc_electric_gnd = DirichletBC(V.sub_space(0), u_electric_gnd, boundary_electric_gnd)

u_electric_potential = Expression('V_cell', V_cell=V_cell, degree=2)
bc_electric_potential = DirichletBC(V.sub_space(2), u_electric_potential, boundary_electric_potential)

bcs = [bc_electric_gnd, bc_electric_potential]

# Need to manually specify these.
dx1 = Measure("dx", domain=V.sub_space(0).mesh())
dx2 = Measure("dx", domain=V.sub_space(1).mesh())
dx3 = Measure("dx", domain=V.sub_space(2).mesh())
ds1 = Measure('ds', domain=V.sub_space(0).mesh(), subdomain_data=boundary_markers1)
ds2 = Measure('ds', domain=V.sub_space(1).mesh(), subdomain_data=boundary_markers2)
ds3 = Measure('ds', domain=V.sub_space(2).mesh(), subdomain_data=boundary_markers3)

# g_L = Constant(-2800.0)
# g_R = Constant(2800.0)
#
# class IAExpression(UserExpression):
#     def eval(self, value, x):
#         value[0] = c
#
#
# class ICExpression(UserExpression):
#     def eval(self, value, x):
#         value[0] = i_c(phis3F(3.5e-4, x[1]), philF(3.5e-4, x[1]))
#
#
# g_L = IAExpression()
# g_R = ICExpression()

# boundary_conditions = {0: {'Dirichlet': u_electric_gnd},
#                        1: {'Neumann': g_L},
#                        2: {'Neumann': g_R},
#                        3: {'Dirichlet': u_electric_potential},
#                        4: {'Neumann': 0}}
#
# integrals_N1 = []
# integrals_N2 = []
# integrals_N3 = []
# for i in boundary_conditions:
#     if 'Neumann' in boundary_conditions[i]:
#         if boundary_conditions[i]['Neumann'] != 0:
#             g = boundary_conditions[i]['Neumann']
#             if i == 1:
#                 integrals_N1.append(g * v1 * ds1(i))
#                 integrals_N2.append(-g * v2 * ds2(i))
#             elif i == 2:
#                 integrals_N2.append(-g * v2 * ds2(i))
#                 integrals_N3.append(g * v3 * ds3(i))
f = Constant(0.0)
e_const = 2.71828
def coth0(x):
    return (e_const**x+e_const**(-x))/(e_const**x-e_const**(-x))
# i_a0 = (218940.044783997*e_const**(65.7479785482249*phil - 65.7479785482249*phis1) + 218940.044783997)
i_a0 = 28000 * phil - 28000 * phis1
i_c0 = 28000 * phis3 - 28000 * phil
a1 = kappa_s * dot(grad(phis1), grad(v1)) * dx1
L1 = f * v1 * dx1 - i_a0 * v1*ds1(1)
a2 = kappa_m * dot(grad(phil), grad(v2)) * dx2
L2 = f * v2 * dx2 + i_a0*v2*ds2(1) - 2800.*v2*ds2(2)
a3 = kappa_s * dot(grad(phis3), grad(v3)) * dx3
L3 = f * v3 * dx3 + 2800.*v3*ds3(2)

a = a1 + a2 + a3
L = L1 + L2 + L3
# Compute solution
phix = Function(V)
F = a1 + a2 + a3 - L1 + L2 + L3
J = derivative(F, phix)
problem = NonlinearVariationalProblem(F, phix, bcs, J)
solver = NonlinearVariationalSolver(problem)
prm = solver.parameters
prm['newton_solver']['absolute_tolerance'] = 1E-8
prm['newton_solver']['relative_tolerance'] = 1E-7
prm['newton_solver']['maximum_iterations'] = 25
prm['newton_solver']['relaxation_parameter'] = 1.0
solver.solve()
phis1F, philF, phis3F = phix.split()

# Save solution to file in VTK format
vtkfile1 = File('main/phis1.pvd')
vtkfile1 << phis1F
vtkfile2 = File('main/phil.pvd')
vtkfile2 << philF
vtkfile3 = File('main/phis3.pvd')
vtkfile3 << phis3F

print(phis1F(2.5e-4, 1e-3), phis1F(2.5e-4, 2e-3))
print(philF(2.5e-4, 1e-3), philF(2.5e-4, 2e-3))
print(philF(3.5e-4, 1e-3), philF(3.5e-4, 2e-3))
print(phis3F(3.5e-4, 1e-3), phis3F(3.5e-4, 2e-3))

