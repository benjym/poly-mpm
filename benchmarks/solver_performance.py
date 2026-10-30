"""Time the updated solver and compare complete numerical snapshots.

Run before edits with --snapshot-dir, then after edits with --reference pointing
at that directory. Output and progress are disabled; diagnostics stay enabled.
"""

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
from time import perf_counter
from types import SimpleNamespace as NS

import numpy as np
from initialise import get_parameters, prepare
from MPM import run

CASES = {"ball": 200, "modern_ball": 200, "apic_10000": 40, "segregation_1000": 40}
SOURCES = (
    "MPM.py",
    "initialise.py",
    "particles.py",
    "grid.py",
    "interpolation.py",
    "transfers.py",
    "kinematics.py",
    "constit.py",
    "stability.py",
    "diagnostics.py",
    "integrators.py",
)


def make_case(name):
    if name in ("ball", "modern_ball"):
        return get_parameters([name])
    segregation = name == "segregation_1000"
    count, ns = (1000, 3) if segregation else (10000, 1)
    rng = np.random.default_rng(41)
    x, y = rng.uniform(0.05, 0.95, (2, count))
    material = NS(
        n=count,
        X=x,
        Y=y,
        A=0.81 / count,
        rho=1000.0,
        law="pouliquen" if segregation else "elastic",
        K=1e4,
        G=4e3,
        rho_s=1666.0,
        mu_s=0.3,
        mu_d=0.7,
        I_0=0.279,
        eta_max=100.0,
        initial_v=np.column_stack((0.01 + 0.1 * x + 0.2 * y, -0.02 - 0.1 * y, np.zeros(count))),
    )
    params = NS(
        interpolation="quadratic",
        transfer="apic",
        diagnostics=True,
        dt=1e-6,
        t_f=1.0,
        savetime=1.0,
        max_g=-1.0,
        pressure="compression",
        initial_pressure=1.0,
        G=NS(nx=41, ny=41, x_m=0.0, x_M=1.0, y_m=0.0, y_M=1.0, ns=ns, s=np.linspace(0.001, 0.003, ns)),
        B=NS(),
        O=NS(),
        S=material,
        segregate=segregation,
        c=1e-4,
    )
    return prepare(params)


def snapshot(P, G, L):
    arrays = {
        f"particle_{i}_{field}": getattr(batch, field).copy()
        for i, batch in enumerate(L.materials)
        for field in batch.persistent_fields
    }
    arrays.update(
        {"grid_" + name: value.copy() for name, value in vars(G).items() if isinstance(value, np.ndarray)}
    )
    arrays["clock"] = np.array([P.t, P.dt, P.tstep])
    return arrays, L.diagnostics.history


def compare(actual, expected, path="diagnostics"):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        return max((compare(actual[k], v, f"{path}.{k}") for k, v in expected.items()), default=0.0)
    if isinstance(expected, list):
        assert len(actual) == len(expected), path
        return max((compare(a, b, path) for a, b in zip(actual, expected)), default=0.0)
    if expected is None or isinstance(expected, (str, bool)):
        assert actual == expected, path
        return 0.0
    np.testing.assert_allclose(actual, expected, rtol=1e-10, atol=1e-12, err_msg=path)
    return abs(actual - expected)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--snapshot-dir", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()
    if args.repeats < 1:
        parser.error("--repeats must be positive")
    args.snapshot_dir.mkdir(parents=True, exist_ok=True)
    root = Path(__file__).resolve().parents[1]
    report = dict(
        environment=dict(
            python=sys.version,
            platform=platform.platform(),
            versions={n: importlib.metadata.version(n) for n in ("numpy", "scipy")},
        ),
        source_sha256={n: hashlib.sha256((root / n).read_bytes()).hexdigest() for n in SOURCES},
        output=False,
        progress=False,
        diagnostics=True,
        warmups=1,
        repeats=args.repeats,
        cases=[],
    )
    for name, steps in CASES.items():
        run(*make_case(name), max_steps=steps, output=False)
        samples = []
        for _ in range(args.repeats):
            P, G, L = make_case(name)
            start = perf_counter()
            run(P, G, L, max_steps=steps, output=False)
            samples.append(perf_counter() - start)
        arrays, history = snapshot(P, G, L)
        np.savez(args.snapshot_dir / f"{name}.npz", **arrays)
        (args.snapshot_dir / f"{name}.json").write_text(json.dumps(history, allow_nan=False))
        row = dict(
            case=name,
            particles=sum(map(len, L.materials)),
            steps=P.tstep,
            simulated_time=P.t,
            dt=P.dt,
            interpolation=P.interpolation,
            transfer=P.transfer,
            samples_seconds=samples,
            median_seconds=float(np.median(samples)),
        )
        if args.reference:
            with np.load(args.reference / f"{name}.npz") as expected:
                assert arrays.keys() == set(expected.files)
                errors = {}
                for field, value in arrays.items():
                    np.testing.assert_allclose(value, expected[field], rtol=1e-10, atol=1e-12, err_msg=field)
                    errors[field] = float(
                        np.max(np.abs(value.astype(float) - expected[field].astype(float)), initial=0)
                    )
            row["numerical_comparison"] = dict(
                rtol=1e-10,
                atol=1e-12,
                arrays=len(arrays),
                max_absolute_error=max(errors.values()),
                max_diagnostic_absolute_error=compare(
                    history, json.loads((args.reference / f"{name}.json").read_text())
                ),
            )
        report["cases"].append(row)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2) + "\n")
        print(f"{name}: {row['median_seconds']:.3f}s ({P.tstep} steps)", flush=True)


if __name__ == "__main__":
    main()
