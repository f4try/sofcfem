"""
FEniCS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.

  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary

  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
from fenics import *
from vtkplotter.dolfin import plot

# Create mesh and define function space
mesh = UnitSquareMesh(8, 8)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_electric_gnd = Expression('0', degree=2)


def boundary_electric_gnd(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 0, tol)
    # return on_boundary and (near(x[0], 0, tol) or near(x[1], 0, tol) or near(x[1], 1, tol))

bc_electric_gnd = DirichletBC(V, u_electric_gnd, boundary_electric_gnd)

u_electric_potential = Expression('0.7', degree=2)


def boundary_electric_potential(x, on_boundary):
    tol = 1E-14
    return on_boundary and near(x[0], 1, tol)

bc_electric_potential = DirichletBC(V, u_electric_potential, boundary_electric_potential)

bc = [bc_electric_gnd]

markers = MeshFunction('size_t', mesh, 2, mesh.domains())
boundary_markers = MeshFunction("size_t", mesh, mesh.topology().dim() - 1, 0)


class BoundaryElectricGnd(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_electric_gnd(x, on_boundary)


class BoundaryElectricPontential(SubDomain):

    def inside(self, x, on_boundary):
        return boundary_electric_potential(x, on_boundary)


class BoundaryInsulationL(SubDomain):

    def inside(self, x, on_boundary):
        return on_boundary and not (boundary_electric_gnd(x, on_boundary) or boundary_electric_potential(x, on_boundary))

bx0 = BoundaryElectricGnd()
bx0.mark(boundary_markers, 0)
bx1 = BoundaryElectricPontential()
bx1.mark(boundary_markers, 1)
bx2 = BoundaryInsulationL()
bx2.mark(boundary_markers, 2)

g = Expression('1e5', degree=1)

ds = Measure('ds', domain=mesh, subdomain_data=boundary_markers)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0.0)
a = dot(grad(u), grad(v))*dx
g = Expression('-1.0', degree=1)
L = f*v*dx - g*v*ds(1)

# Compute solution
u = Function(V)
solve(a == L, u, bc)


# Save solution to file in VTK format
vtkfile = File('neumann/solution.pvd')
vtkfile << u

print(u(0.5,0.5),u(0.5,0.5),u(1,0.5))
# Plot solution and mesh
plot(u)