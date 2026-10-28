# Maintaining the documentation

The current manual starts at [getting started](getting_started.rst); the root
README provides a short first-run route. Sphinx builds the `.rst`
and `.md` sources together, so the model audit, example matrix, and conservation
report appear in the same navigation as the tutorial and API reference.

From the checkout root:

```sh
uv sync --extra dev --extra docs --locked
uv run --locked --extra docs python -m sphinx -b html -W --keep-going docs /tmp/poly-mpm-docs/html
```

Open `/tmp/poly-mpm-docs/html/index.html` in a browser. On Windows use a writable
local output directory instead of `/tmp`. The equivalent Makefile invocation,
after activating an environment with the docs dependencies installed, is
`make -C docs html`. Its default output is `docs/build/html/`.

The repository retains generated HTML under `docs/build/html/` and the
`docs/index.html` redirect used by the existing static-site layout. To refresh
that copy, use the same strict build with `docs/build/html` as its output directory.
Review the generated diff before committing. Building or committing locally does
not publish a website; deployment is a separate action.

Keep documentation versions sourced from `pyproject.toml` in `conf.py`; do not
maintain a separate hand-written release number. Build warnings are errors.
Do not commit Sphinx doctree caches. The API reference imports solver modules,
so install runtime dependencies as well as the docs extra.

When changing commands or examples:

1. Execute the first-run preview and the tutorial input in a temporary output directory.
2. Check saved PNGs/NPY shapes and the stated output location.
3. Run the strict Sphinx build and check local links in the generated HTML.
4. If material or example support changes, regenerate `benchmarks/examples.json`
   and update the example matrix and model guidance together. Preserve historical
   performance reports as dated measurements rather than relabelling their environment.

The optional external-link check is `make -C docs linkcheck`; it needs network
access and external sites can rate-limit it. This is separate from the offline
HTML build and internal-link verification.
