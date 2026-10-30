from pathlib import Path
from types import SimpleNamespace as NS
import numpy as np
import pytest
from numpy.testing import assert_allclose
from conftest import state
from MPM import run
from plotting import Plotting
from initialise import get_parameters


@pytest.mark.parametrize("regular_plots", [False, True])
def test_grain_and_strain_rate_frames_increment_independently(tmp_path, regular_plots):
    figures = {"gsd_material_points", "gsd_grid", "shear_rate"}
    if regular_plots:
        figures |= {"material_points", "continuum"}
    O = NS(figures=figures)
    P, G, L = state(ns=2, O=O, output_dir=str(tmp_path), savetime=1e-5)
    # One save-time output and the final output must produce separate frames.
    run(P, G, L, max_steps=1)
    assert sorted(p.name for p in (tmp_path / "GSD").glob("MP_*.png")) == [
        "MP_00000.png",
        "MP_00001.png",
    ]
    assert sorted(p.name for p in (tmp_path / "GSD").glob("[0-9]*.png")) == [
        "00000.png",
        "00001.png",
    ]
    assert sorted(p.name for p in (tmp_path / "MP").glob("gammadot_*.png")) == [
        "gammadot_00000.png",
        "gammadot_00001.png",
    ]


def test_output_flags_and_saved_shapes(tmp_path):
    O = NS(
        figures=("continuum", "material_points", "gsd_material_points", "gsd_grid"),
        fields=("u", "s_bar", "density", "phi_material_points"),
        measure_energy=True,
        measure_stiffness=True,
    )
    P, G, L = state(ns=3, O=O, output_dir=str(tmp_path), savetime=1e-5)
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


def test_custom_callbacks_replace_defaults(tmp_path):
    events = []
    O = NS(
        at_start=lambda L, P, G: events.append("start"),
        after_every_timestep=lambda L, P, G: events.append("step"),
        after_every_savetime=lambda L, P, G: events.append("save"),
        at_end=lambda L, P, G: events.append("end"),
        figures=("continuum",),
    )
    P, G, L = state(O=O, output_dir=str(tmp_path), savetime=1e-5)
    run(P, G, L, max_steps=2)
    assert events == ["start", "step", "save", "step", "save", "end"]
    assert not list(tmp_path.rglob("*.png"))


def test_custom_save_creates_directory_without_flags(tmp_path):
    O = NS(after_every_savetime=lambda L, P, G: Plotting().save_u(L, P, G))
    P, G, L = state(O=O, output_dir=str(tmp_path / "nested"), savetime=1e-5)
    run(P, G, L, max_steps=1)
    assert (tmp_path / "nested/data/u_00000.npy").is_file()


def test_ball_custom_callback(tmp_path):
    P, G, L = get_parameters(["ball"])
    P.save_dir = str(tmp_path)
    P.savetime = P.dt
    run(P, G, L, max_steps=1)
    assert list((tmp_path / "Continuum").glob("*.png"))
    assert list((tmp_path / "MP").glob("*.png"))


def test_reused_plots_refresh_data_colors_and_keep_layout(tmp_path):
    import matplotlib.pyplot as plt

    P, G, L = state(output_dir=str(tmp_path))
    run(P, G, L, max_steps=1, output=False)
    plot = Plotting()
    managers = plt.get_fignums()
    plot.draw_continuum(G, P)
    plot.draw_material_points(L, P, G)
    figures = {kind: cached[1] for kind, cached in plot._figures.items()}
    positions = {kind: [ax.get_position().bounds for ax in fig.axes] for kind, fig in figures.items()}

    # Change masks, positions and color ranges, including a constant field.
    G.valid_nodes[:] = True
    G.valid_nodes[0] = False
    G.s_bar[:] = np.linspace(10, 20, len(G.s_bar))
    for b in L.materials:
        b.x[:, 0] += 0.1
        b.v[:, 0] = np.linspace(10, 20, len(b))
        b.pressure[:] = -50
    plot.draw_continuum(G, P)
    plot.draw_material_points(L, P, G)

    # A fresh render must have the same data, masks and normalization as reuse.
    fresh = Plotting()
    fresh.draw_continuum(G, P)
    fresh.draw_material_points(L, P, G)
    for kind, count in (("continuum", 6), ("particles", 3)):
        fig = plot._figures[kind][1]
        assert fig is figures[kind]
        assert len(fig.axes) == 2 * count
        assert_allclose([ax.get_position().bounds for ax in fig.axes], positions[kind])
        for ax, reference in zip(fig.axes[:count], fresh._figures[kind][1].axes[:count]):
            artist, expected = ax.collections[0], reference.collections[0]
            assert_allclose(artist.get_array(), expected.get_array())
            assert np.array_equal(
                np.ma.getmaskarray(artist.get_array()), np.ma.getmaskarray(expected.get_array())
            )
            assert_allclose(artist.get_clim(), expected.get_clim())
            if kind == "particles":
                assert_allclose(artist.get_offsets(), expected.get_offsets())
    assert plt.get_fignums() == managers
    assert (tmp_path / "Continuum/00001.png").is_file()
    assert (tmp_path / "MP/00001.png").is_file()


def test_constitutive_history_callback(tmp_path):
    P, G, L = get_parameters(["pouliquen_unit_test"])
    P.save_dir = str(tmp_path)
    run(P, G, L, max_steps=3)
    with np.load(tmp_path / "data/constitutive_history.npz") as history:
        assert history["mu"].shape == (3, 9)
        assert history["dstrain"].shape == (3, 9, 3, 3)
        assert np.isfinite(history["pressure"]).all()
    assert (tmp_path / "0.png").is_file()
