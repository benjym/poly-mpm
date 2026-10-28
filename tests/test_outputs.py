from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
from numpy.testing import assert_allclose
from conftest import state
from MPM import run
from plotting import Plotting
from initialise import get_parameters


def test_output_flags_and_saved_shapes(tmp_path):
    O = NS(
        plot_continuum=True,
        plot_material_points=True,
        plot_gsd_mp=True,
        plot_gsd_grid=True,
        save_u=True,
        save_s_bar=True,
        save_density=True,
        save_phi_MP=True,
        measure_energy=True,
        measure_stiffness=True,
    )
    P, G, L = state(ns=3, O=O, save_dir=str(tmp_path), savetime=1e-5)
    run(P, G, L, max_steps=2)
    for name in ("u_00000.npy", "v_00000.npy", "density_00000.npy", "s_bar_00000.npy"):
        assert np.load(tmp_path / "data" / name).shape == (5, 5)
    assert np.load(tmp_path / "data/MP_phi_00000.npy").shape == (4, 5)
    assert np.load(tmp_path / "data/energy.npy").shape == (2, 4)
    from PIL import Image

    images = list(tmp_path.rglob("*.png"))
    assert len(images) >= 5
    for path in images:
        with Image.open(path) as image:
            assert image.width > 20 and image.height > 20
            assert np.asarray(image).std() > 0


def test_custom_callbacks_replace_defaults_and_alias(tmp_path):
    events = []
    O = NS(
        at_start=lambda L, P, G: events.append("start"),
        after_ever_timestep=lambda L, P, G: events.append("step"),
        after_every_savetime=lambda L, P, G: events.append("save"),
        at_end=lambda L, P, G: events.append("end"),
        plot_continuum=True,
    )
    P, G, L = state(O=O, save_dir=str(tmp_path), savetime=1e-5)
    run(P, G, L, max_steps=2)
    assert events == ["start", "step", "save", "step", "save", "end"]
    assert not list(tmp_path.rglob("*.png"))


def test_custom_save_creates_directory_without_flags(tmp_path):
    O = NS(after_every_savetime=lambda L, P, G: Plotting().save_u(L, P, G))
    P, G, L = state(O=O, save_dir=str(tmp_path / "nested"), savetime=1e-5)
    run(P, G, L, max_steps=1)
    assert (tmp_path / "nested/data/u_00000.npy").is_file()


def test_ball_custom_callback(tmp_path):
    P, G, L = get_parameters(["ball"])
    P.save_dir = str(tmp_path)
    P.savetime = P.dt
    run(P, G, L, max_steps=1)
    assert list((tmp_path / "Continuum").glob("*.png"))
    assert list((tmp_path / "MP").glob("*.png"))


def test_constitutive_history_callback(tmp_path):
    P, G, L = get_parameters(["pouliquen_unit_test"])
    P.save_dir = str(tmp_path)
    run(P, G, L, max_steps=3)
    with np.load(tmp_path / "data/constitutive_history.npz") as history:
        assert history["mu"].shape == (3, 9)
        assert history["dstrain"].shape == (3, 9, 3, 3)
        assert np.isfinite(history["pressure"]).all()
    assert (tmp_path / "0.png").is_file()
