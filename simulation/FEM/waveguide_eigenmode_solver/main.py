# ------------------------------
# Dependencies
# ------------------------------
from mpi4py import MPI
from petsc4py import PETSc
import numpy as np
import gmsh
import ufl

from dolfinx.fem import (
    Function, functionspace, dirichletbc,
    locate_dofs_topological, form, Constant
)
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import XDMFFile
import dolfinx


# ------------------------------
# Parameters
# ------------------------------
from config import k0, Lx, Ly, poly_deg, n_mode

c = 1
ky = n_mode * np.pi / Ly
kx = np.sqrt(k0**2 - ky**2)
omega = c / k0

# ------------------------------
# Mesh
# ------------------------------
gmsh.initialize()
gmsh.open("mesh/mesh_delta.msh")

physical_name_to_tag = {}
for dim, tag in gmsh.model.getPhysicalGroups():
    name = gmsh.model.getPhysicalName(dim, tag)
    physical_name_to_tag[name] = tag

mesh_data = dolfinx.io.gmsh.model_to_mesh(
    gmsh.model,
    MPI.COMM_WORLD,
    rank=0,
    gdim=2
)
gmsh.finalize()

mesh = mesh_data.mesh
facet_tags = mesh_data.facet_tags

V = functionspace(mesh, ("Lagrange", poly_deg))

# ------------------------------
# Boundaries
# ------------------------------
mirror_tag = physical_name_to_tag["mirrors"]
inlet_tag = physical_name_to_tag["inlet"]
outlet_tag = physical_name_to_tag["outlet"]

mirror_facets = facet_tags.find(mirror_tag)
inlet_facets = facet_tags.find(inlet_tag)
outlet_facets = facet_tags.find(outlet_tag)

dofs_mirrors = locate_dofs_topological(V, mesh.topology.dim-1, mirror_facets)

# ------------------------------
# BCs (ONLY mirrors now)
# ------------------------------
bcs = []

u_mirrors = Function(V)
u_mirrors.x.array[:] = 0.0
bcs.append(dirichletbc(u_mirrors, dofs_mirrors))

# ------------------------------
# PML condition
# ------------------------------
L_PML = 0.05 * Lx
x_coords = mesh.geometry.x
x_start = Lx - L_PML

# Damping profile PML
sigma_max = 0.1 * k0
sigma = np.zeros(x_coords.shape[0])

mask_pml = x_coords[:,0] > x_start
sigma[mask_pml] = sigma_max * ((x_coords[mask_pml,0]-x_start)/L_PML)**2

sigma_f = Function(V)
sigma_f.x.array[:] = sigma

# ------------------------------
# Variational formulation
# ------------------------------
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)

ds = ufl.Measure("ds", domain=mesh, subdomain_data=facet_tags)

# --- Helmholtz operator ---
a = (
    ufl.inner(ufl.grad(u)/(1 + 1j * sigma_f / omega), ufl.grad(v)) * ufl.dx
    - k0**2 * ufl.inner(u, v) * ufl.dx
)

# ------------------------------
# Sommerfeld outlet BC (weak form)
# ------------------------------
#a += 1j * k0 * ufl.inner(u, v) * ds(outlet_tag)

# ------------------------------
# Modal source at inlet
# ------------------------------
g = Function(V)

def inlet_profile(x):
    values = np.zeros(x.shape[1], dtype=np.complex128)
    values[:] = np.sin(ky * x[1])
    return values

g.interpolate(inlet_profile)

L = ufl.inner(g, v) * ds(inlet_tag)

# ------------------------------
# Solve
# ------------------------------
uh = Function(V, name="u")

problem = LinearProblem(
    a,
    L,
    bcs=bcs,
    u=uh,
	petsc_options_prefix = "waveguide_sim",
    petsc_options={
        "ksp_type": "gmres",
        "pc_type": "ilu",
    },
)

problem.solve()

# ------------------------------
# Output
# ------------------------------
with XDMFFile(MPI.COMM_WORLD, "out/helmholtz_scattering.xdmf", "w") as file:
    file.write_mesh(mesh)
    file.write_function(uh)

if MPI.COMM_WORLD.rank == 0:
    print("Helmholtz scattering simulation complete.")
    print(f"Mode ky = {ky}, kx = {kx}")


# ------------------------------
# Robust Reflection / Transmission Analysis
# ------------------------------

def compute_modal_amplitude_slice(uh, V, x_slice, Lx, L_PML, Ly, n_mode, slice_width=0.05):
    """
    Compute modal amplitude at a slice of the waveguide.

    slice_width: fraction of Lx to collect enough DOFs
    """
    coords = V.tabulate_dof_coordinates()
    # Collect DOFs in the slice
    dofs_slice = np.where(
        (coords[:,0] >= x_slice - slice_width/2) & (coords[:,0] <= x_slice + slice_width/2)
    )[0]

    if len(dofs_slice) < 2:
        raise ValueError(f"Slice at x={x_slice:.3f} too thin or no DOFs found. Increase slice_width.")

    y_slice = coords[dofs_slice, 1]
    u_slice = uh.x.array[dofs_slice]

    # Define transverse mode
    phi = np.sin(n_mode * np.pi * y_slice / Ly)

    # Manual trapezoidal integration
    def trapz_manual(f, y):
        return np.sum(0.5 * (f[1:] + f[:-1]) * np.diff(y))

    # Mode projection
    numerator = trapz_manual(u_slice * np.conj(phi), y_slice)
    denominator = trapz_manual(np.abs(phi)**2, y_slice)

    if denominator == 0:
        raise ValueError(f"Mode normalization failed at x={x_slice:.3f}. Check slice DOFs.")

    return numerator / denominator

# --- select slices safely outside PML ---
x_in = 0.05 * Lx                 # just right of inlet
x_out = 0.5 * Lx        # just before PML

A_in = compute_modal_amplitude_slice(uh, V, x_in, Lx, L_PML, Ly, n_mode)
A_out = compute_modal_amplitude_slice(uh, V, x_out, Lx, L_PML, Ly, n_mode)

R = np.abs(A_in)**2
T = np.abs(A_out)**2

print(f"Reflection R = {R:.4f}, Transmission T = {T:.4f}, R+T = {R+T:.4f}")
