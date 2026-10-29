"""Checkpoint and restart of a prepared run.

A checkpoint holds the evolving state of P, G, every particle batch and the
diagnostics: arrays, numbers, strings and plain containers of them. Inputs,
materials, callbacks and controllers are not stored; restart rebuilds them
from the same input with initialise.get_parameters and then loads the state:

    P, G, L = get_parameters(args)
    checkpoint.load(path, P, G, L)
    run(P, G, L)

On the CPU the continued run is identical to an uninterrupted one. A Metal run
re-uploads the host state synchronised at the checkpoint; in the frictional-wall
column it also reproduced the uninterrupted run exactly.
"""

import pickle
from pathlib import Path

import numpy as np

VERSION = 1
# Configuration objects and functions held by the evolving objects.
SKIP = {"material", "law", "templates", "materials", "diagnostics"}


def _plain(value):
    if value is None or isinstance(value, (bool, int, float, str, np.ndarray, np.generic)):
        return True
    if isinstance(value, (list, tuple)):
        return all(_plain(v) for v in value)
    if isinstance(value, dict):
        return all(isinstance(k, str) and _plain(v) for k, v in value.items())
    return False


def _state(obj):
    return {
        name: value
        for name, value in vars(obj).items()
        if name not in SKIP and not callable(value) and _plain(value)
    }


def _restore(obj, state):
    for name, value in state.items():
        setattr(obj, name, value.copy() if isinstance(value, np.ndarray) else value)


def save(path, P, G, L):
    """Write the evolving state of a prepared run to path (a pickle)."""
    state = dict(
        version=VERSION,
        args=getattr(P, "args", None),
        P=_state(P),
        G=_state(G),
        materials=[_state(batch) for batch in L.materials],
        diagnostics=None if L.diagnostics is None else _state(L.diagnostics),
    )
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    with open(temporary, "wb") as handle:
        pickle.dump(state, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def load(path, P, G, L):
    """Load a checkpoint into a run prepared from the same input and arguments.

    Only load checkpoints you wrote: the file is a pickle.
    """
    with open(path, "rb") as handle:
        state = pickle.load(handle)
    if state.get("version") != VERSION:
        raise ValueError(f"checkpoint version {state.get('version')} is not {VERSION}")
    if len(state["materials"]) != len(L.materials):
        raise ValueError("checkpoint and prepared run have different material counts")
    if (state["diagnostics"] is None) != (L.diagnostics is None):
        raise ValueError("checkpoint and prepared run differ in whether diagnostics are enabled")
    _restore(P, state["P"])
    _restore(G, state["G"])
    for batch, values in zip(L.materials, state["materials"]):
        _restore(batch, values)
    if L.diagnostics is not None:
        _restore(L.diagnostics, state["diagnostics"])
    return P, G, L
