"""Compare corrected scalar reference and array timesteps on identical state.

Run: python -m benchmarks.compare --output benchmarks/results.json
Reference construction and deepcopy are outside timed intervals.
"""
import argparse
from copy import deepcopy
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
import tracemalloc
from types import SimpleNamespace as NS
import numpy as np
from initialise import prepare
from particles import ParticleSystem

# Test-only scalar reference is deliberately excluded from the installed solver.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tests'))
from reference import ScalarParticles


def make_state(count, law, ns=1, segregate=False):
    rng = np.random.default_rng(41)
    S = NS(n=count, X=rng.uniform(.05, .95, count), Y=rng.uniform(.05, .95, count),
           A=.81/count, rho=1000., law=law, K=1e4, G=4e3, k=500., s=2.5,
           mu_s=10., mu_v=5., rho_s=1666., mu_0=.3, mu_1=.7, delta_mu=.4,
           I_0=.279, eta_max=100., packing=.6, initial_v=[.01, -.02, 0])
    if law != 'rigid':
        S.initial_v = np.column_stack((.01+.1*S.X+.2*S.Y, -.02-.1*S.Y, np.zeros(count)))
    P = NS(dt=1e-6, t_f=1., savetime=1., max_g=-1., pressure='compression', initial_pressure=1.,
           G=NS(nx=41, ny=41, x_m=0., x_M=1., y_m=0., y_M=1., ns=ns, s=np.linspace(.001, .003, ns)),
           B=NS(), O=NS(), S=S, segregate=segregate, c=1e-4)
    return prepare(P)


def stages(P, G, L):
    times = {}
    start = perf_counter()
    G.wipe(P)
    P.update_forces()
    L.get_reference_nodes(P, G)
    L.get_basis_functions(P, G)
    L.get_nodal_mass_momentum(P, G)
    times['transfers'] = perf_counter()-start
    start = perf_counter()
    L.update_stress_strain(P, G)
    times['constitutive_and_diagnostics'] = perf_counter()-start
    start = perf_counter()
    L.get_nodal_forces(P, G)
    times['forces'] = perf_counter()-start
    start = perf_counter()
    G.BCs(P)
    G.update_momentum(P)
    G.calculate_gammadot(P, G)
    times['grid'] = perf_counter()-start
    start = perf_counter()
    if P.segregate:
        G.update_pk(P, G)
        G.calculate_phi_increment(P)
        L.move_grainsize_on_grid(P, G)
    times['segregation'] = perf_counter()-start
    start = perf_counter()
    L.move_material_points(P, G)
    L.process_boundaries(P, G)
    times['movement'] = perf_counter()-start
    times['total'] = sum(times.values())
    return times


def compare(count, law, ns, segregate, repeats):
    initial = make_state(count, law, ns, segregate)
    row = dict(particles=count, law=law, grain_classes=ns, segregation=segregate)
    outcomes = {}
    for implementation in ('scalar', 'array'):
        template = deepcopy(initial)
        if implementation == 'scalar':
            # Same state, alternative kernels; no change in numerical setup.
            for batch in template[2].materials:
                batch.__class__ = ScalarParticles
        stages(*deepcopy(template))
        samples = []
        for _ in range(repeats):
            current = deepcopy(template)
            samples.append(stages(*current))
        row[implementation] = {name: float(np.median([s[name] for s in samples])) for name in samples[0]}
        outcomes[implementation] = current
    for field in ('x', 'v', 'stress', 'phi'):
        np.testing.assert_allclose(getattr(outcomes['array'][2].materials[0], field),
                                   getattr(outcomes['scalar'][2].materials[0], field), rtol=1e-8, atol=1e-12)
    memory_state = deepcopy(initial)
    tracemalloc.start()
    stages(*memory_state)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    row['array_timestep_peak_allocated_bytes'] = peak
    row['particle_state_bytes'] = sum(getattr(initial[2].materials[0], name).nbytes
                                      for name in initial[2].materials[0].persistent_fields)
    row['speedup'] = row['scalar']['total']/row['array']['total']
    return row


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--sizes', type=int, nargs='+', default=[1000, 10000, 100000])
    parser.add_argument('--repeats', type=int, default=5)
    parser.add_argument('--output', default='benchmarks/results.json')
    parser.add_argument('--laws', nargs='+', default=['elastic', 'viscous', 'von_mises', 'pouliquen', 'rigid'])
    args = parser.parse_args()
    report = dict(environment=dict(python=sys.version, platform=platform.platform(),
                                  versions={name: importlib.metadata.version(name) for name in ['numpy', 'scipy', 'matplotlib', 'astropy']}),
                  repeats=args.repeats, warmups=1, cases=[])
    for count in args.sizes:
        for law in args.laws:
            cases = [(1, False)] + ([(2, True), (8, True)] if law == 'pouliquen' else [])
            for ns, segregation in cases:
                row = compare(count, law, ns, segregation, args.repeats)
                report['cases'].append(row)
                Path(args.output).write_text(json.dumps(report, indent=2)+'\n')
                print(f'{count} {law} ns={ns}: {row["speedup"]:.1f}x; array={row["array"]["total"]:.5f}s scalar={row["scalar"]["total"]:.5f}s', flush=True)
