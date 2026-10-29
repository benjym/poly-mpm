"""Finite-domain integration of grid basis functions (point stress per domain)."""

import numpy as np


def _kernel(x, degree):
    a = np.abs(x)
    if degree == 1:
        return np.maximum(1 - a, 0)
    return np.where(a < 0.5, 0.75 - x * x, np.where(a < 1.5, 0.5 * (1.5 - a) ** 2, 0.0))


def _axis(center, half, lo, h, count, degree, periodic):
    radius = 1.0 if degree == 1 else 1.5
    lower = (center - half - lo) / h
    upper = (center + half - lo) / h
    if not np.isfinite(lower).all() or not np.isfinite(upper).all() or np.any(upper - lower > 124):
        raise ValueError("particle domain spans more than 128 grid nodes on an axis")
    first = np.floor(lower - radius).astype(int)
    last = np.ceil(upper + radius).astype(int)
    size = int((last - first + 1).max(initial=0))
    if size > 128:
        raise ValueError("particle domain spans more than 128 grid nodes on an axis")
    raw = first[:, None] + np.arange(size)
    center_offset = (center - lo)[:, None] / h - raw
    radius = half[:, None] / h
    weights = np.zeros_like(center_offset)
    gradients = np.zeros_like(weights)
    knots = (-1.0, 0.0, 1.0) if degree == 1 else (-1.5, -0.5, 0.5, 1.5)

    def boundary(value):
        distance = value - center_offset
        result = np.sign(distance)
        return np.divide(distance, radius, out=result, where=np.abs(distance) < radius)

    for low, high in zip(knots, knots[1:]):
        a, b = boundary(low), boundary(high)
        fraction = 0.5 * (b - a)
        middle = 0.5 * (a + b)
        shift = (b - a) / (2 * np.sqrt(3))
        for xi in (middle - shift, middle + shift):
            distance = center_offset + radius * xi
            weights += 0.5 * fraction * _kernel(distance, degree)
            if degree == 1:
                derivative = np.where(np.abs(distance) < 1, -np.sign(distance), 0)
            else:
                derivative = np.where(
                    np.abs(distance) < 0.5,
                    -2 * distance,
                    np.where(np.abs(distance) < 1.5, -np.sign(distance) * (1.5 - np.abs(distance)), 0),
                )
            gradients += 0.5 * fraction * derivative / h
    if periodic:
        return raw % (count - 1), weights, gradients, raw * h + lo
    start = np.maximum(0, np.minimum(first, count - 2))
    end = np.minimum(count - 1, np.maximum(last, 1))
    width = int((end - start + 1).max(initial=0))
    nodes = np.minimum(start[:, None] + np.arange(width), count - 1)
    w = np.zeros_like(nodes, dtype=float)
    g = np.zeros_like(w)
    row = np.broadcast_to(np.arange(len(center))[:, None], raw.shape)
    active = raw <= last[:, None]
    primary = np.clip(raw, 0, count - 1)
    coefficient = np.where(raw < 0, 1 - raw, np.where(raw >= count, raw - count + 2, 1))
    for result, values in ((w, weights), (g, gradients)):
        np.add.at(result, (row[active], (primary - start[:, None])[active]), (coefficient * values)[active])
        secondary = np.where(raw < 0, 1, count - 2)
        coefficient2 = np.where(raw < 0, raw, np.where(raw >= count, count - 1 - raw, 0))
        # Only exterior entries have a secondary contribution; avoid invalid padding indices.
        outside = ((raw < 0) | (raw >= count)) & active
        np.add.at(
            result, (row[outside], (secondary - start[:, None])[outside]), (coefficient2 * values)[outside]
        )
    return nodes, w, g, nodes * h + lo


def _aligned_stencil(position, F, halfsize, P, G):
    half = halfsize * np.abs(F[:, [0, 1], [0, 1]])
    if np.any(half <= 0) or not np.isfinite(half).all():
        raise ValueError("particle domain has collapsed or nonfinite extents")
    degree = 1 if P.interpolation == "linear" else 2
    ix, wx, gx, xx = _axis(position[:, 0], half[:, 0], P.G.x_m, G.dx, P.G.nx, degree, P.B.cyclic_lr)
    iy, wy, gy, yy = _axis(position[:, 1], half[:, 1], P.G.y_m, G.dy, P.G.ny, degree, False)
    n = len(position)
    kx = ix.shape[1]
    ky = iy.shape[1]
    if kx * ky > getattr(P, "domain_max_support", 256):
        raise ValueError("particle domain exceeds domain_max_support; refine or split the particle domain")
    nodes = (iy[:, :, None] * P.G.nx + ix[:, None, :]).reshape(n, -1)
    w = (wy[:, :, None] * wx[:, None, :]).reshape(n, -1)
    g = np.zeros((*nodes.shape, 3))
    coords = np.zeros_like(g)
    g[:, :, 0] = (wy[:, :, None] * gx[:, None, :]).reshape(n, -1)
    g[:, :, 1] = (gy[:, :, None] * wx[:, None, :]).reshape(n, -1)
    coords[:, :, 0] = np.broadcast_to(xx[:, None, :], (n, ky, kx)).reshape(n, -1)
    coords[:, :, 1] = np.broadcast_to(yy[:, :, None], (n, ky, kx)).reshape(n, -1)
    return nodes, w, g, coords - position[:, None, :]


def _clip(poly, axis, cut, keep_lower):
    result = []
    previous = poly[-1]
    previous_inside = (previous[axis] <= cut) if keep_lower else (previous[axis] >= cut)
    for current in poly:
        inside = (current[axis] <= cut) if keep_lower else (current[axis] >= cut)
        if inside != previous_inside:
            fraction = (cut - previous[axis]) / (current[axis] - previous[axis])
            result.append(previous + fraction * (current - previous))
        if inside:
            result.append(current)
        previous, previous_inside = current, inside
    return np.asarray(result)


def _knots(low, high, origin, h, count, degree, periodic):
    shift = 0.5 if degree == 2 else 0.0
    first = int(np.floor((low - origin) / h - shift)) + 1
    last = int(np.ceil((high - origin) / h - shift))
    if not periodic:
        first = max(first, 0 if degree == 2 else 1)
        last = min(last, count - 1)
    if last - first > 128:
        raise ValueError("particle domain spans more than 128 grid knots")
    indices = np.arange(first, last)
    return origin + (indices + shift) * h


def _points(xy, P, G):
    from interpolation import _quadratic_axis

    axes = []
    degree = 1 if P.interpolation == "linear" else 2
    for axis, (origin, h, count, periodic) in enumerate(
        ((P.G.x_m, G.dx, P.G.nx, P.B.cyclic_lr), (P.G.y_m, G.dy, P.G.ny, False))
    ):
        if degree == 2:
            index, w, g, coordinate = _quadratic_axis(xy[:, axis], origin, h, count, periodic)
        else:
            position = (xy[:, axis] - origin) / h
            cell = np.floor(position).astype(int)
            if not periodic:
                cell = np.clip(cell, 0, count - 2)
            raw = cell[:, None] + [0, 1]
            local = position - cell
            w = np.column_stack((1 - local, local))
            g = np.broadcast_to([-1 / h, 1 / h], w.shape)
            coordinate = origin + raw * h
            index = raw % (count - 1) if periodic else raw
        raw = np.rint((coordinate - origin) / h).astype(int)
        axes.append((index, w, g, coordinate, raw))
    ix, wx, gx, xx, rawx = axes[0]
    iy, wy, gy, yy, rawy = axes[1]
    n = len(xy)
    kx = ix.shape[1]
    ky = iy.shape[1]
    keys = np.stack(
        (np.broadcast_to(rawx[:, None, :], (n, ky, kx)), np.broadcast_to(rawy[:, :, None], (n, ky, kx))),
        axis=-1,
    ).reshape(-1, 2)
    w = (wy[:, :, None] * wx[:, None, :]).reshape(n, -1)
    g = np.stack((wy[:, :, None] * gx[:, None, :], gy[:, :, None] * wx[:, None, :]), axis=-1).reshape(
        n, -1, 2
    )
    return keys, w, g


_BARY = np.array(
    [
        [0.445948490915965, 0.445948490915965, 0.108103018168070],
        [0.445948490915965, 0.108103018168070, 0.445948490915965],
        [0.108103018168070, 0.445948490915965, 0.445948490915965],
        [0.091576213509771, 0.091576213509771, 0.816847572980459],
        [0.091576213509771, 0.816847572980459, 0.091576213509771],
        [0.816847572980459, 0.091576213509771, 0.091576213509771],
    ]
)
_TRIANGLE_WEIGHTS = np.r_[np.full(3, 0.223381589678011), np.full(3, 0.109951743655322)]


def _polygon_stencil(center, F, half, P, G):
    corners = (
        np.array([[-half[0], -half[1]], [half[0], -half[1]], [half[0], half[1]], [-half[0], half[1]]])
        @ F[:2, :2].T
        + center[:2]
    )
    extent = np.ptp(corners, axis=0) / [G.dx, G.dy]
    if not np.isfinite(extent).all() or np.any(extent > 124):
        raise ValueError("particle domain spans more than 128 grid nodes on an axis")
    area = 4 * np.prod(half) * np.linalg.det(F[:2, :2])
    pieces = [corners]
    degree = 1 if P.interpolation == "linear" else 2
    for axis, (origin, h, count, periodic) in enumerate(
        ((P.G.x_m, G.dx, P.G.nx, P.B.cyclic_lr), (P.G.y_m, G.dy, P.G.ny, False))
    ):
        cuts = _knots(corners[:, axis].min(), corners[:, axis].max(), origin, h, count, degree, periodic)
        if len(cuts) > 128:
            raise ValueError("particle domain spans more than 128 grid knots")
        for cut in cuts:
            split = []
            for poly in pieces:
                if poly[:, axis].min() < cut < poly[:, axis].max():
                    split.extend((_clip(poly, axis, cut, True), _clip(poly, axis, cut, False)))
                else:
                    split.append(poly)
            pieces = split
            if len(pieces) > getattr(P, "domain_max_support", 256):
                raise ValueError("particle domain exceeds domain_max_support during knot integration")
    points = []
    weights = []
    for poly in pieces:
        for j in range(1, len(poly) - 1):
            triangle = poly[[0, j, j + 1]]
            u, v = triangle[1] - triangle[0], triangle[2] - triangle[0]
            triangle_area = 0.5 * abs(u[0] * v[1] - u[1] * v[0])
            if triangle_area:
                points.append(_BARY @ triangle)
                weights.append(_TRIANGLE_WEIGHTS * triangle_area / area)
    if not points:
        raise ValueError("particle domain geometry is unresolved at floating-point precision")
    points = np.concatenate(points)
    weights = np.concatenate(weights)
    keys, w, g = _points(points, P, G)
    unique, inverse = np.unique(keys, axis=0, return_inverse=True)
    w = np.bincount(inverse, weights=(weights[:, None] * w).ravel(), minlength=len(unique))
    gradient = np.zeros((len(unique), 3))
    for axis in (0, 1):
        gradient[:, axis] = np.bincount(
            inverse, weights=(weights[:, None] * g[:, :, axis]).ravel(), minlength=len(unique)
        )
    active = (w != 0) | np.any(gradient != 0, axis=1)
    unique, w, gradient = unique[active], w[active], gradient[active]
    nodes = unique[:, 0] % (P.G.nx - 1) if P.B.cyclic_lr else unique[:, 0]
    nodes = nodes + P.G.nx * unique[:, 1]
    coordinates = np.column_stack(
        (P.G.x_m + unique[:, 0] * G.dx, P.G.y_m + unique[:, 1] * G.dy, np.zeros(len(unique)))
    )
    return nodes, w, gradient, coordinates - center


def _clip_batch(vertices, count, axis, bound, keep_lower):
    """Sutherland-Hodgman clip of many convex polygons against one half-plane each.

    vertices is (m, 8, 2) with count valid vertices per row; a convex polygon
    clipped by a rectangle keeps at most eight. Rows emptied by the clip keep
    count 0.
    """
    slots = vertices.shape[1]
    index = np.arange(slots)
    valid = index < count[:, None]
    previous = (index - 1) % np.maximum(count, 1)[:, None]
    before = np.take_along_axis(vertices, previous[:, :, None], axis=1)
    side = (
        (lambda p: p[..., axis] <= bound[:, None])
        if keep_lower
        else (lambda p: p[..., axis] >= bound[:, None])
    )
    inside, before_inside = side(vertices), side(before)
    crossing = (inside != before_inside) & valid
    span = vertices[..., axis] - before[..., axis]
    fraction = np.divide(bound[:, None] - before[..., axis], span, out=np.zeros_like(span), where=crossing)
    candidates = np.stack((before + fraction[..., None] * (vertices - before), vertices), axis=2)
    keep = np.stack((crossing, inside & valid), axis=2).reshape(len(vertices), -1)
    order = np.argsort(~keep, axis=1, kind="stable")[:, :slots]
    clipped = np.take_along_axis(candidates.reshape(len(vertices), -1, 2), order[:, :, None], axis=1)
    return clipped, keep.sum(axis=1)


def _intervals(low, high, origin, h, count, periodic, degree):
    """Knot intervals covering [low, high] for many domains at once (as _knots).

    Returns the interval ends, flattened over domains, with each domain's first
    interval index and interval count.
    """
    shift = 0.5 if degree == 2 else 0.0
    first = np.floor((low - origin) / h - shift).astype(int) + 1
    last = np.ceil((high - origin) / h - shift).astype(int)
    if not periodic:
        first = np.maximum(first, 0 if degree == 2 else 1)
        last = np.minimum(last, count - 1)
    if np.any(last - first > 128):
        raise ValueError("particle domain spans more than 128 grid knots")
    knots = np.maximum(last - first, 0)
    domain = np.repeat(np.arange(len(low)), knots)
    step = np.arange(knots.sum()) - np.repeat(np.cumsum(knots) - knots, knots)
    cuts = origin + (first[domain] + step + shift) * h
    inside = (cuts > low[domain]) & (cuts < high[domain])
    domain, cuts = domain[inside], cuts[inside]
    # Boundaries per domain: low, interior cuts, high; sorted within each domain.
    owner = np.r_[np.arange(len(low)), domain, np.arange(len(low))]
    value = np.r_[low, cuts, high]
    rank = np.r_[np.zeros(len(low)), np.ones(len(cuts)), np.full(len(low), 2.0)]
    order = np.lexsort((value, rank, owner))
    owner, value = owner[order], value[order]
    pair = owner[:-1] == owner[1:]
    lows, highs = value[:-1][pair], value[1:][pair]
    counts = np.bincount(owner[:-1][pair], minlength=len(low))
    return lows, highs, np.cumsum(counts) - counts, counts


def _polygon_stencils(center, F, half, P, G):
    """All rotated/sheared domains at once; same integrals as _polygon_stencil.

    Each domain is cut into its intersections with the basis knot cells, which
    are clipped together, fan-triangulated and integrated by the same
    degree-four triangle rule. Returns one (nodes, weights, gradients, offsets)
    row per domain, as _polygon_stencil does.
    """
    unit = np.array([[-1.0, -1.0], [1.0, -1.0], [1.0, 1.0], [-1.0, 1.0]])
    corners = np.einsum("kd,nd,njd->nkj", unit, half, F[:, :2, :2]) + center[:, None, :2]
    spacing = np.array([G.dx, G.dy])
    extent = np.ptp(corners, axis=1) / spacing
    if not np.isfinite(extent).all() or np.any(extent > 124):
        raise ValueError("particle domain spans more than 128 grid nodes on an axis")
    area = 4 * np.prod(half, axis=1) * np.linalg.det(F[:, :2, :2])
    degree = 1 if P.interpolation == "linear" else 2
    intervals = [
        _intervals(corners[:, :, axis].min(axis=1), corners[:, :, axis].max(axis=1), *grid, degree)
        for axis, grid in enumerate(((P.G.x_m, G.dx, P.G.nx, P.B.cyclic_lr), (P.G.y_m, G.dy, P.G.ny, False)))
    ]
    (x_lo, x_hi, x_start, x_count), (y_lo, y_hi, y_start, y_count) = intervals
    total = x_count * y_count
    if np.any(total > getattr(P, "domain_max_support", 256)):
        raise ValueError("particle domain exceeds domain_max_support during knot integration")
    owner = np.repeat(np.arange(len(center)), total)
    local = np.arange(total.sum()) - np.repeat(np.cumsum(total) - total, total)
    xi = x_start[owner] + local // y_count[owner]
    yi = y_start[owner] + local % y_count[owner]
    cells = np.column_stack((x_lo[xi], x_hi[xi], y_lo[yi], y_hi[yi]))
    vertices = np.zeros((len(owner), 8, 2))
    vertices[:, :4] = corners[owner]
    count = np.full(len(owner), 4)
    for axis, column, keep_lower in ((0, 0, False), (0, 1, True), (1, 2, False), (1, 3, True)):
        vertices, count = _clip_batch(vertices, count, axis, cells[:, column], keep_lower)
    # Fan triangles (0, j, j + 1) of each clipped piece.
    j = np.arange(1, 7)
    present = j[None, :] + 1 < count[:, None]
    piece, fan = np.nonzero(present)
    a = vertices[piece, 0]
    b = vertices[piece, fan + 1]
    c = vertices[piece, fan + 2]
    u, v = b - a, c - a
    triangle_area = 0.5 * np.abs(u[:, 0] * v[:, 1] - u[:, 1] * v[:, 0])
    nonzero = triangle_area > 0
    piece, a, b, c, triangle_area = piece[nonzero], a[nonzero], b[nonzero], c[nonzero], triangle_area[nonzero]
    particle = owner[piece]
    points = (
        _BARY[None, :, 0, None] * a[:, None]
        + _BARY[None, :, 1, None] * b[:, None]
        + _BARY[None, :, 2, None] * c[:, None]
    ).reshape(-1, 2)
    weights = (_TRIANGLE_WEIGHTS[None, :] * (triangle_area / area[particle])[:, None]).ravel()
    particle = np.repeat(particle, len(_TRIANGLE_WEIGHTS))
    if np.any(np.bincount(particle, minlength=len(center)) == 0):
        raise ValueError("particle domain geometry is unresolved at floating-point precision")
    keys, w, g = _points(points, P, G)
    support = w.shape[1]
    # One integer per (particle, knot x, knot y): a 1-D sort, not a row sort.
    low = keys.min(axis=0)
    width = keys.max(axis=0) - low + 1
    labels = (np.repeat(particle, support) * width[1] + keys[:, 1] - low[1]) * width[0] + keys[:, 0] - low[0]
    codes, inverse = np.unique(labels, return_inverse=True)
    rest, kx = np.divmod(codes, width[0])
    owner, ky = np.divmod(rest, width[1])
    unique = np.column_stack((owner, kx + low[0], ky + low[1]))
    w = np.bincount(inverse, weights=(weights[:, None] * w).ravel(), minlength=len(unique))
    gradient = np.zeros((len(unique), 3))
    for axis in (0, 1):
        gradient[:, axis] = np.bincount(
            inverse, weights=(weights[:, None] * g[:, :, axis]).ravel(), minlength=len(unique)
        )
    active = (w != 0) | np.any(gradient != 0, axis=1)
    unique, w, gradient = unique[active], w[active], gradient[active]
    nodes = unique[:, 1] % (P.G.nx - 1) if P.B.cyclic_lr else unique[:, 1]
    nodes = nodes + P.G.nx * unique[:, 2]
    coordinates = np.column_stack(
        (P.G.x_m + unique[:, 1] * G.dx, P.G.y_m + unique[:, 2] * G.dy, np.zeros(len(unique)))
    )
    offsets = coordinates - center[unique[:, 0]]
    starts = np.searchsorted(unique[:, 0], np.arange(len(center) + 1))
    return [(nodes[s:e], w[s:e], gradient[s:e], offsets[s:e]) for s, e in zip(starts[:-1], starts[1:])]


def domain_stencil(position, F, halfsize, P, G):
    """Integrate over F-convected rectangles; polynomial pieces are integrated exactly.

    Axis-aligned domains use exact piecewise polynomial Gauss integration. Rotated/sheared domains
    are clipped at basis knots and integrated by a degree-four triangle rule.
    """
    if not len(position):
        k = 4 if P.interpolation == "linear" else 9
        return np.zeros((0, k), dtype=int), np.zeros((0, k)), np.zeros((0, k, 3)), np.zeros((0, k, 3))
    if not np.isfinite(F).all() or not np.isfinite(halfsize).all() or np.any(halfsize <= 0):
        raise ValueError("particle domain geometry must be finite and positive")
    if (
        np.any(np.abs(F[:, :2, 2]) > 1e-12)
        or np.any(np.abs(F[:, 2, :2]) > 1e-12)
        or np.any(np.abs(F[:, 2, 2] - 1) > 1e-12)
    ):
        raise ValueError("domain quadrature requires planar deformation and fixed thickness")
    determinant = np.linalg.det(F[:, :2, :2])
    if np.any(determinant <= 0):
        raise ValueError("particle domain has nonpositive area")
    tolerance = 32 * np.finfo(float).eps
    aligned = (np.abs(F[:, 0, 1]) * halfsize[:, 1] <= tolerance * G.dx) & (
        np.abs(F[:, 1, 0]) * halfsize[:, 0] <= tolerance * G.dy
    )
    if aligned.all():
        return _check(_aligned_stencil(position, F, halfsize, P, G), P)
    rows = {}
    if aligned.any():
        result = _aligned_stencil(position[aligned], F[aligned], halfsize[aligned], P, G)
        for j, i in enumerate(np.flatnonzero(aligned)):
            active = (result[1][j] != 0) | np.any(result[2][j] != 0, axis=1)
            rows[i] = tuple(a[j, active] for a in result)
    general = np.flatnonzero(~aligned)
    # "reference" keeps the per-domain clipping path that the batched one is tested against.
    if getattr(P, "domain_polygon", "batched") == "reference":
        for i in general:
            rows[i] = _polygon_stencil(position[i], F[i], halfsize[i], P, G)
    else:
        rows.update(zip(general, _polygon_stencils(position[general], F[general], halfsize[general], P, G)))
    k = max(len(r[0]) for r in rows.values())
    if k > getattr(P, "domain_max_support", 256):
        raise ValueError("particle domain exceeds domain_max_support; refine or split the particle domain")
    nodes = np.zeros((len(position), k), dtype=int)
    w = np.zeros_like(nodes, dtype=float)
    gradient = np.zeros((*nodes.shape, 3))
    offsets = np.zeros_like(gradient)
    for i, (n, weight, g, r) in rows.items():
        size = len(n)
        nodes[i, :size], w[i, :size], gradient[i, :size], offsets[i, :size] = n, weight, g, r
    return _check((nodes, w, gradient, offsets), P)


def _check(result, P):
    nodes, w, gradient, offsets = result
    if w.shape[1] > getattr(P, "domain_max_support", 256):
        raise ValueError("particle domain exceeds domain_max_support; refine or split the particle domain")
    if not np.isfinite(w).all() or not np.isfinite(gradient).all():
        raise FloatingPointError("nonfinite integrated particle basis")
    if np.any(w < -128 * np.finfo(float).eps):
        raise ValueError(
            "particle domain extends too far beyond a physical grid boundary (negative integrated weights)"
        )
    if not np.allclose(w.sum(axis=1), 1, rtol=0, atol=2e-11):
        raise FloatingPointError(
            "particle domain integration lost partition of unity; geometry is unresolved"
        )
    # Negative values only at roundoff are zeroed; no weight renormalisation.
    w[w < 0] = 0
    return nodes, w, gradient, offsets
