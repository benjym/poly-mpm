"""Measure the existing KT/remapping scheme's composition drift without clipping."""
import argparse
import json
from pathlib import Path
import numpy as np
from initialise import get_parameters
from MPM import run


def measure(ns, steps=100):
    P, G, L = get_parameters(['bi_seg_test', '22', str(ns), '11'])
    b = L.materials[0]
    before = np.sum(b.m[:, None]*b.phi, axis=0)
    run(P, G, L, max_steps=steps, output=False)
    after = np.sum(b.m[:, None]*b.phi, axis=0)
    return dict(grain_classes=ns, steps=P.tstep, time=P.t, normalise_phi=P.normalise_phi,
                minimum_fraction=float(b.phi.min()), maximum_fraction=float(b.phi.max()),
                maximum_fraction_sum_error=float(np.max(np.abs(b.phi.sum(axis=1)-1))),
                initial_species_mass=before.tolist(), final_species_mass=after.tolist(),
                species_mass_change=(after-before).tolist())


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', default='benchmarks/segregation_drift.json')
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps([measure(ns) for ns in [2, 3, 8]], indent=2)+'\n')
