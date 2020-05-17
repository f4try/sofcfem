from dolfin import Constant, between, File, Function, inner, Measure, MixedFunctionSpace, \
    MeshFunction, nabla_grad, near, SubDomain, RectangleMesh, Point, MeshView, FunctionSpace, DirichletBC,\
    TestFunctions, TrialFunctions, solve

# Width of total domain
width = 1.0
# Height of total domain
height = 1.0
# Total resolution
resolution = 16
# Forcing function
f = Constant(10)
# Boundary values for the two equations
a_0 = 1
b_0 = 1


# Define SubDomain
class LeftDomain(SubDomain):
    def inside(self, x, on_boundary):
        return between(x[0], (0, height/2))


# Define Boundary
class LeftBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and near(x[0], 0.0)


# Mesh for the entire domain
mesh_full = RectangleMesh(Point(0, 0),
                          Point(width, height),
                          resolution,
                          resolution)

# Define domains
domains = MeshFunction("size_t", mesh_full, 2)
domains.set_all(0)
left_domain = LeftDomain()
left_domain.mark(domains, 1)

# Define boundaries
boundaries_full = MeshFunction("size_t", mesh_full, 1)
boundaries_full.set_all(0)
left_boundary = LeftBoundary()
left_boundary.mark(boundaries_full, 1)

# Create mesh for the left domain
mesh_left = MeshView.create(domains, 1)

# Define boundaries for left domain
boundaries_left = MeshFunction("size_t", mesh_left, 1)
boundaries_left.set_all(0)
left_boundary_2 = LeftBoundary()
left_boundary_2.mark(boundaries_left, 1)

# Create the function spaces
W_a = FunctionSpace(mesh_full, "Lagrange", 1)
W_b = FunctionSpace(mesh_left, "Lagrange", 1)
W = MixedFunctionSpace(W_a, W_b)

# Define boundary conditions
bc_a = DirichletBC(W.sub_space(0), a_0, boundaries_full, 1)
bc_b = DirichletBC(W.sub_space(1), b_0, boundaries_left, 1)
bcs = [bc_a, bc_b]

# Create all of the weak form variables
u = Function(W)
ua, ub = TrialFunctions(W)
va, vb = TestFunctions(W)

# Need to manually specify these.
dxa = Measure("dx", domain=W.sub_space(0).mesh())
dxb = Measure("dx", domain=W.sub_space(1).mesh())

# The weak formulation
a = inner(nabla_grad(ua), nabla_grad(va)) * dxa + inner(nabla_grad(ub), nabla_grad(vb)) * dxb
L = f * va * dxa + f * vb * dxb

solve(a == L, u, bcs)

File("sd/ua.pvd") << u.sub(0)
File("sd/ub.pvd") << u.sub(1)