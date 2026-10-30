"""The settings schema: rejection of unknown names, defaults, checks and translation."""

import re
from pathlib import Path
from types import SimpleNamespace as NS

import numpy as np
import pytest

import settings
from conftest import grid, material, state
from initialise import prepare
from MPM import run


def params(**overrides):
    """Unprepared settings of a small elastic block."""
    P = NS(dt=1e-5, t_f=0.01, savetime=0.01, G=grid(), B=NS(), O=NS(), S=material())
    P.__dict__.update(overrides)
    return P


def test_unknown_setting_is_rejected_with_a_suggestion():
    P = params()
    P.G.thicknes = 0.1
    with pytest.raises(ValueError, match=r"unknown setting G\.thicknes; did you mean 'thickness'\?"):
        prepare(P)


@pytest.mark.parametrize("name", ["max_g", "CFL", "has_bottom", "update_forces"])
def test_internal_names_are_not_settings(name):
    P = params(**{name: 1.0})
    with pytest.raises(ValueError, match=name):
        prepare(P)


def test_private_names_and_class_methods_are_left_alone():
    events = []

    class Output:
        def at_end(self, L, P, G):
            events.append(P.t)

    P = params(_helper=3.0, O=Output())
    P, G, L = prepare(P)
    assert P._helper == 3.0
    run(P, G, L, max_steps=1)
    assert events == [P.t]


def test_missing_required_setting_is_reported():
    P = params()
    del P.G.nx
    with pytest.raises(ValueError, match=r"G\.nx is required"):
        prepare(P)


@pytest.mark.parametrize(
    "where, name, value, message",
    [
        ("P", "time_stepping", "dynamc", "time_stepping must be one of 'static', 'dynamic'"),
        ("B", "wall_contact", "rough", "wall_contact must be one of"),
        ("B", "walls", {"floor"}, "walls must be a collection drawn from"),
        ("O", "figures", {"particles"}, "figures must be a collection drawn from"),
        ("P", "dt", -1.0, r"dt must be > 0"),
    ],
)
def test_invalid_value_is_reported(where, name, value, message):
    P = params()
    setattr(P if where == "P" else getattr(P, where), name, value)
    with pytest.raises(ValueError, match=message):
        prepare(P)


def test_material_field_the_law_does_not_read_is_rejected():
    P = params(S=material("elastic", yield_stress=1.0))
    with pytest.raises(ValueError, match=r"unknown setting S\[0\] \(elastic\)\.yield_stress"):
        prepare(P)


def test_update_forces_on_an_input_is_rejected():
    class Params:
        def __init__(self):
            self.__dict__.update(vars(params()))

        def update_forces(self):
            self.g = -9.81

    with pytest.raises(ValueError, match="update_forces"):
        prepare(Params())


def test_forces_are_functions_of_time():
    P, G, L = state(
        gravity=-10.0,
        gravity_ramp=lambda t: min(t / 1e-4, 1.0),
        gravity_angle=lambda t: 100 * t,
        B=NS(traction_vertical=lambda t: 1e5 * t, traction_horizontal=2.0),
    )
    assert P.g == 0.0 and P.theta == 0.0 and P.q_v == 0.0 and P.q_h == 2.0
    run(P, G, L, max_steps=3, output=False)
    # Forces are evaluated at the start of each step.
    t = P.t - P.dt
    assert np.isclose(P.g, -10.0 * t / 1e-4)
    assert np.isclose(P.theta, 100 * t)
    assert np.isclose(P.q_v, 1e5 * t)
    assert P.q_h == 2.0


def test_boundary_and_output_settings_translate_to_solver_flags():
    P, G, L = state(
        B=NS(walls={"bottom", "left"}, wall_contact="no_slip", obstacles={"centre_wall"}, outlets={"bottom"}),
        O=NS(figures={"continuum", "shear_rate"}, fields={"u"}),
    )
    B, O = P.B, P.O
    assert B.has_bottom and B.has_left and not B.has_top and not B.has_right
    assert B.roughness and not B.no_slip_bottom and not B.wall_mu
    assert B.wall and not B.box and B.outlet_bottom and not B.outlet_left
    assert O.plot_continuum and O.plot_gamma_dot and not O.plot_material_points
    assert O.save_u and not O.save_density


def test_coulomb_walls_need_a_friction_coefficient():
    with pytest.raises(ValueError, match="wall_friction"):
        state(B=NS(walls={"bottom"}, wall_contact="coulomb"))
    P, _, _ = state(B=NS(walls={"bottom"}, wall_contact="coulomb", wall_friction=0.3))
    assert P.B.wall_mu == 0.3 and not P.B.roughness


def test_features_tables_are_generated_from_the_schema():
    text = (Path(__file__).resolve().parents[1] / "docs/features.md").read_text()
    for where in settings.SETTINGS:
        found = re.search(rf"<!-- settings:{where} -->\n(.*?)\n<!-- /settings:{where} -->", text, re.S)
        assert found, where
        assert found.group(1) == settings.markdown(where)
