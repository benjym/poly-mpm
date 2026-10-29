"""Resident grid segregation for the Metal granular step (integrators.segregation_update).

Static face lists are built once on the host with the CPU scheme's own
orientation and padding helpers, so the device stencils, wall closures and
periodic seams are those of integrators._face_transfers. Everything that
evolves (fractions, deposits, gradients, fluxes, limiters and particle
exchanges) is computed on the device.
"""

import numpy as np

from .primitives import MAX, SUM, ColumnReduction

MODULE = "segregation.slang"
NS_MAX = 8


def faces(P, G):
    """(faces[F, 6], geometry[F, 2], node_faces[N, 4]) in integrators._face_transfers order."""
    from integrators import _oriented, _pad

    index = np.arange(len(G.m))
    rows, geometry = [], []
    for axis in (0, 1):
        periodic = axis == 0 and P.B.cyclic_lr
        nodes = _pad(_oriented(index, P, axis), periodic)
        wall = _oriented(G.boundary_v if axis == 0 else G.boundary_h, P, axis).copy()
        if not periodic:
            wall[0] = wall[-1] = False
        wall = _pad(wall, periodic)
        left, right = nodes[1:-2], nodes[2:-1]
        before, after = nodes[:-3], nodes[3:]
        closed = wall[1:-2] | wall[2:-1]
        if axis == 0:
            area = G.dy_plot * P.G.thickness
            spacing = G.dx
        else:
            area = (np.full(P.G.nx, G.dx) if P.B.cyclic_lr else G.dx_plot) * P.G.thickness
            spacing = G.dy
        area = np.broadcast_to(area[None, :], left.shape)
        selection = slice(None, -1) if periodic else slice(1, -1)
        arrays = [a[selection] for a in (left, right, before, after, closed, area)]
        if axis == 1 and P.B.cyclic_lr:
            arrays = [a[:, :-1] for a in arrays]
        left, right, before, after, closed, area = (a.ravel() for a in arrays)
        distinct = left != right
        count = int(distinct.sum())
        rows.append(
            np.column_stack(
                (
                    left[distinct],
                    right[distinct],
                    before[distinct],
                    after[distinct],
                    np.full(count, axis),
                    closed[distinct],
                )
            )
        )
        geometry.append(np.column_stack((np.full(count, spacing), area[distinct])))
    table = np.vstack(rows).astype(np.uint32)
    geometry = np.vstack(geometry)
    incident = np.zeros((len(G.m), 4), dtype=np.int32)
    filled = np.zeros(len(G.m), dtype=int)
    for f, (a, b) in enumerate(table[:, :2]):
        for node, sign in ((a, 1), (b, -1)):
            if filled[node] == 4:
                raise ValueError("a grid node has more than four segregation faces")
            incident[node, filled[node]] = sign * (f + 1)
            filled[node] += 1
    return table, geometry, incident.view(np.uint32)


class Segregation:
    def __init__(self, rt, grid, P, G, count):
        ns = int(P.G.ns)
        if ns > NS_MAX:
            raise ValueError(f"Metal segregation supports at most {NS_MAX} grain sizes")
        self.rt, self.grid, self.ns, self.count = rt, grid, ns, count
        N = grid.nodes
        f32 = np.float32

        def vector(size):
            return rt.buffer(count=max(int(size), 1), dtype=f32)

        table, geometry, incident = faces(P, G)
        self.face_count = len(table)
        self.faces = rt.buffer(table if len(table) else np.zeros((1, 6), np.uint32))
        self.geometry = rt.buffer(geometry if len(geometry) else np.zeros((1, 2)), dtype=f32)
        self.node_faces = rt.buffer(incident)
        self.sizes = rt.buffer(np.asarray(P.G.s, dtype=float), dtype=f32)
        self.nodal_volume = rt.buffer(G.V, dtype=f32)
        self.phi = vector(ns * count)
        self.phi_start = vector(ns * count)
        self.limit = vector(ns * count)
        self.phim, self.phi_nodes = vector(ns * N), vector(ns * N)
        self.pressure, self.inertial = vector(N), vector(N)
        self.grad_pk, self.diffusivity = vector(2 * N), vector(N)
        self.transfer = vector(ns * self.face_count)
        self.ratio, self.node_limit = vector(ns * N), vector(ns * N)
        self.coefficients, self.change = vector(ns * ns * N), vector(ns * N)
        self.columns = vector(3 * N)
        self.species_columns = vector(ns * count)
        self.reductions = [ColumnReduction(rt, N, (MAX, MAX, MAX)), ColumnReduction(rt, count, (SUM,) * ns)]

    def upload(self, b, P):
        self.rt.write(self.phi, b.phi, np.float32)
        self.rt.write(self.sizes, np.asarray(P.G.s, dtype=float), np.float32)

    def encode(
        self, encoder, P, operator, nodal, velocity, particle_mass, stress, gradient, density, status, dt
    ):
        rt, grid, ns, n, N = self.rt, self.grid, self.ns, self.count, self.grid.nodes
        periodic = int(grid.periodic)
        S = P.S[0]
        bins = operator.bins
        rt.dispatch(
            encoder, "seg_copy", ns * n, module=MODULE, count=ns * n, source=self.phi, target=self.phi_start
        )
        rt.dispatch(
            encoder,
            "seg_nodes",
            N,
            module=MODULE,
            node_count=N,
            nx=grid.nx,
            ny=grid.ny,
            periodic=periodic,
            ns=ns,
            rho_s=float(S.rho_s),
            root_rho_s=float(np.sqrt(S.rho_s)),
            cell_offsets=bins.offsets,
            cell_counts=bins.counts,
            permutation=bins.permutation,
            local_x=operator.local_x,
            local_y=operator.local_y,
            particle_mass=particle_mass,
            phi=self.phi,
            sizes=self.sizes,
            stress=stress,
            gradient=gradient,
            density=density,
            phim=self.phim,
            pressure=self.pressure,
            inertial=self.inertial,
        )
        rt.dispatch(
            encoder,
            "seg_fields",
            N,
            module=MODULE,
            node_count=N,
            nx=grid.nx,
            ny=grid.ny,
            ns=ns,
            dx=grid.dx,
            dy=grid.dy,
            m_tol=float(P.M_tol),
            mass=nodal.mass,
            velocity=velocity,
            pressure=self.pressure,
            inertial=self.inertial,
            phim=self.phim,
            sizes=self.sizes,
            phi_nodes=self.phi_nodes,
            grad_pk=self.grad_pk,
            diffusivity=self.diffusivity,
            columns=self.columns,
        )
        if self.face_count:
            rt.dispatch(
                encoder,
                "seg_faces",
                self.face_count,
                module=MODULE,
                face_count=self.face_count,
                ns=ns,
                dt=float(dt),
                c_seg=float(P.c),
                m_tol=float(P.M_tol),
                faces=self.faces,
                face_geometry=self.geometry,
                mass=nodal.mass,
                nodal_volume=self.nodal_volume,
                grad_pk=self.grad_pk,
                diffusivity=self.diffusivity,
                sizes=self.sizes,
                phi_nodes=self.phi_nodes,
                transfer=self.transfer,
            )
            rt.dispatch(
                encoder,
                "seg_outgoing",
                N,
                module=MODULE,
                node_count=N,
                ns=ns,
                node_faces=self.node_faces,
                mass=nodal.mass,
                phi_nodes=self.phi_nodes,
                transfer=self.transfer,
                ratio=self.ratio,
            )
            corners = dict(
                nx=grid.nx,
                periodic=periodic,
                ns=ns,
                m_tol=float(P.M_tol),
                keys=operator.keys,
                local_x=operator.local_x,
                local_y=operator.local_y,
                mass=nodal.mass,
            )
            rt.dispatch(
                encoder,
                "seg_demand",
                n,
                module=MODULE,
                count=n,
                phi=self.phi,
                ratio=self.ratio,
                limit=self.limit,
                **corners,
            )
            rt.dispatch(
                encoder,
                "seg_node_limit",
                N,
                module=MODULE,
                node_count=N,
                nx=grid.nx,
                ny=grid.ny,
                periodic=periodic,
                ns=ns,
                m_tol=float(P.M_tol),
                cell_offsets=bins.offsets,
                cell_counts=bins.counts,
                permutation=bins.permutation,
                local_x=operator.local_x,
                local_y=operator.local_y,
                mass=nodal.mass,
                phi=self.phi,
                phi_nodes=self.phi_nodes,
                limit=self.limit,
                node_limit=self.node_limit,
            )
            rt.dispatch(
                encoder,
                "seg_limit_faces",
                self.face_count,
                module=MODULE,
                face_count=self.face_count,
                ns=ns,
                faces=self.faces,
                node_limit=self.node_limit,
                transfer=self.transfer,
            )
            rt.dispatch(
                encoder,
                "seg_coefficients",
                N,
                module=MODULE,
                node_count=N,
                nx=grid.nx,
                periodic=periodic,
                ns=ns,
                node_faces=self.node_faces,
                mass=nodal.mass,
                phi_nodes=self.phi_nodes,
                transfer=self.transfer,
                ratio=self.ratio,
                coefficients=self.coefficients,
                change=self.change,
            )
            rt.dispatch(
                encoder,
                "seg_particles",
                n,
                module=MODULE,
                count=n,
                normalise=int(bool(P.normalise_phi)),
                ratio=self.ratio,
                coefficients=self.coefficients,
                phi=self.phi,
                status=status,
                **{k: v for k, v in corners.items() if k != "mass"},
                mass=nodal.mass,
            )
        rt.dispatch(
            encoder,
            "seg_species",
            n,
            module=MODULE,
            count=n,
            ns=ns,
            mass=particle_mass,
            phi=self.phi,
            columns=self.species_columns,
        )
        return [
            self.reductions[0].encode(encoder, self.columns, status),
            self.reductions[1].encode(encoder, self.species_columns, status),
        ]

    def particles(self, read):
        """End-of-step fractions and the start-of-step fractions used for I."""
        shape = (self.count, self.ns)
        return read(self.phi, shape), read(self.phi_start, shape)

    def grid_fields(self, P, G, read):
        """Grid fields of the CPU segregation step (start-of-step fractions, dphi)."""
        from transfers import divide

        ns, N = self.ns, self.grid.nodes
        G.phim = read(self.phim, (N, ns))
        G.phi = divide(G.phim, G.m[:, None])
        G.s_bar = G.phi @ P.G.s
        G.dphi = read(self.change, (N, ns))
        G.segregation_m = G.m.copy()
        G.segregation_phi = G.phi.copy()
        G.segregation_particle_dphi = None
