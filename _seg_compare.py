import numpy as np
from copy import deepcopy
from initialise import prepare
from MPM import run
from inputs.granular_column import Params
P = Params([])
P.granular_solver, P.granular_tolerance, P.diagnostics = "velocity", 1e-4, True
P.theta = np.arctan(0.75); P.t_f = 0.2
P.G.ns, P.G.s = 2, [0.001, 0.003]
P.S.PHI = np.tile([0.5, 0.5], (P.S.n, 1))
P.segregate, P.c = True, 1e-2
state = prepare(P)
state[0].B.cyclic_lr = state[0].B.roughness = True; state[1].boundary_v[:] = False
cpu, metal = deepcopy(state), deepcopy(state)
metal[0].backend = "metal"
for s in (cpu, metal): run(*s, output=False)
bc, bm = cpu[2].materials[0], metal[2].materials[0]
change = np.abs(bc.phi - state[2].materials[0].phi).max()
print("steps", cpu[0].tstep, metal[0].tstep)
print("phi change (CPU)", change, "| phi difference Metal-CPU", np.abs(bm.phi - bc.phi).max(), "relative to change", np.abs(bm.phi - bc.phi).max() / change)
print("species CPU", (bc.m[:, None] * bc.phi).sum(0) - (state[2].materials[0].m[:, None] * state[2].materials[0].phi).sum(0))
print("species Metal (FP64 sum of snapshot)", (bm.m[:, None] * bm.phi).sum(0) - (state[2].materials[0].m[:, None] * state[2].materials[0].phi).sum(0), "phi sum error", np.abs(bm.phi.sum(1) - 1).max())
print("positions dx/h", np.abs(bm.x - bc.x).max() / cpu[1].dx)
for f in ("phi", "dphi", "pk", "shear_rate", "pressure", "I"):
    a, c = getattr(cpu[1], f), getattr(metal[1], f)
    print("G." + f, "max |CPU|", np.abs(a).max(), "max diff", np.abs(a - c).max())
print("grad_pk", np.abs(cpu[1].grad_pk).max(), np.abs(cpu[1].grad_pk - metal[1].grad_pk).max())
print("limits", cpu[0].timestep_limits, metal[0].timestep_limits)
print("diag species", cpu[2].diagnostics.last["species_mass_error"], metal[2].diagnostics.last["species_mass_error"])
