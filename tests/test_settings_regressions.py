"""Settings that once had no effect, or the wrong one."""

import contextlib
import io

import numpy as np
import pytest

import checkpoint
from initialise import get_parameters
from MPM import run


def _parameters(args):
    with contextlib.redirect_stdout(io.StringIO()):
        return get_parameters(args)


def test_traction_ramp_starts_from_zero():
    # Inputs ramp q_v and q_h in update_forces; they once ramped an unread P.q, so the
    # boundary load stayed at max_q from the first step.
    P, G, L = _parameters(["triaxial"])
    assert P.B.vertical_force and P.max_q > 0
    assert P.q_v == 0.0
    run(P, G, L, max_steps=3, output=False)
    assert 0 < P.q_v < P.max_q


@pytest.mark.parametrize("roughness", [False, True])
def test_explicit_wall_friction_acts_without_roughness(roughness):
    P, G, L = _parameters(["ball"])
    P.B.wall_mu, P.B.roughness = 0.3, roughness
    bottom = np.arange(1, P.G.nx - 1)
    G.m[:] = 1.0
    G.q[:] = 0.0
    G.fi[:] = 0.0
    G.fe[:] = 0.0
    G.fe[bottom, 0] = 2.0  # tangential
    G.fe[bottom, 1] = -1.0  # into the floor
    G.update_momentum(P)
    # Coulomb: the tangential force is reduced by at most wall_mu times the normal force.
    assert np.allclose(G.q_dot[bottom, 0], 2.0 - 0.3 * 1.0)
    assert np.allclose(G.q_dot[bottom, 1], 0.0)


def test_unknown_time_stepping_is_rejected():
    from inputs.ball import Params
    from initialise import prepare

    P = Params(["ball"])
    P.time_stepping = "dynamc"
    with pytest.raises(ValueError, match="time_stepping"):
        prepare(P)


def test_restart_keeps_settings_changed_in_the_input(tmp_path):
    P, G, L = _parameters(["ball"])
    run(P, G, L, max_steps=3, output=False)
    path = tmp_path / "checkpoint.pkl"
    checkpoint.save(path, P, G, L)
    restarted = _parameters(["ball"])
    restarted[0].t_f = 0.5 * P.t_f
    checkpoint.load(path, *restarted)
    assert restarted[0].t_f == 0.5 * P.t_f
    assert restarted[0].tstep == 3


def test_restart_rejects_a_checkpoint_from_other_arguments(tmp_path):
    P, G, L = _parameters(["ball"])
    path = tmp_path / "checkpoint.pkl"
    checkpoint.save(path, P, G, L)
    restarted = _parameters(["ball"])
    restarted[0].args = ["ball", "other"]
    with pytest.raises(ValueError, match="written from"):
        checkpoint.load(path, *restarted)
