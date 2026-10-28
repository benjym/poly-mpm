"""Bounded compatibility audit; run with python -m benchmarks.audit_examples.

Small smoke dt is a compatibility check, not a long-run stability claim.
"""

import argparse
import contextlib
import io
import json
from pathlib import Path
from initialise import get_parameters
from MPM import run

ARGUMENTS = {
    "bi_dam_break": ["11", "mix", "0", "2"],
    "bi_double_insets": ["11"],
    "bi_drum": ["11", "2"],
    "bi_seg_test": ["22", "2", "11"],
    "bi_square_drum": ["11", "2"],
    "bi_vertical_chute": ["2", "11"],
    "francois": ["11"],
    "poly_seg_test": ["22"],
    "silo": ["11"],
}


def audit(steps=5):
    rows = []
    for path in sorted((Path(__file__).resolve().parents[1] / "inputs").glob("*.py")):
        if path.stem == "__init__":
            continue
        args = [path.stem] + ARGUMENTS.get(path.stem, [])
        row = dict(example=path.stem, args=args)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                P, G, L = get_parameters(args)
                P.dt = min(P.dt, 1e-6)
                run(P, G, L, max_steps=steps, output=False)
            row.update(status="supported smoke-tested", particles=[len(b) for b in L.materials])
        except Exception as exc:
            row.update(
                status=(
                    "unsupported material law"
                    if "Unsupported material law" in str(exc)
                    else "known legacy setup defect"
                ),
                error=f"{type(exc).__name__}: {exc}",
            )
        rows.append(row)
        print(f"{row['example']}: {row['status']} {row.get('error', '')}", flush=True)
    return rows


if __name__ == "__main__":
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", default="benchmarks/examples.json")
    parser.add_argument("--steps", type=int, default=5)
    args = parser.parse_args()
    Path(args.output).write_text(json.dumps(audit(args.steps), indent=2) + "\n")
