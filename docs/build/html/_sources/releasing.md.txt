# Building and publishing packages

The project builds a pure-Python wheel and source distribution. A GitHub Actions
job builds both, checks their metadata/README with Twine, and installs the wheel
in an isolated environment for a bounded solver run. CI validates artifacts; the separate manual publishing workflow uploads a release only when dispatched. The repository version is currently `1.0.0`; confirm that version and
project-name ownership before the first public release.

From a clean checkout:

```sh
uv sync --locked --extra dev --extra docs --extra release
uv run --locked black --check .
uv run --locked python -m pytest -q -W error::RuntimeWarning --cov --cov-report=xml
uv run --locked python -m sphinx -b html -W --keep-going docs build/manual
uv build --no-sources
uv run --locked --extra release python -m twine check --strict dist/*
```

Use an empty `dist/` directory for each release so old artifacts cannot be
uploaded accidentally. Install both distributions in fresh environments outside
the checkout, check imports and `poly-mpm modern_ball`, and inspect the rendered
package description. A full example produces output; use a temporary directory.
CI checks the wheel; the publishing workflow checks both the wheel and sdist.

## Trusted Publishing workflow

`.github/workflows/pypi.yml` builds, checks, and publishes using PyPI Trusted
Publishing. The configured pending publisher uses `benjym/poly-mpm`, workflow
`pypi.yml`, and **no GitHub environment**. The publishing job deliberately has
no `environment` field so its OIDC identity matches that configuration.
Only that job receives `id-token: write`; it downloads the tested distributions
without checking out or executing project source. No API token is required.

To publish an approved version:

1. Set the intended version in `pyproject.toml`, update `uv.lock`, and push the
   tested release source to `main`. Inspect the package description and contents;
   check documentation links as rendered on PyPI.
2. Open **Actions → Publish to PyPI → Run workflow**, select `main`, and enter
   the exact version from `pyproject.toml`.
3. The workflow checks formatting/tests, builds wheel and sdist, checks metadata,
   and installs/runs each distribution outside the checkout. Publication runs
   only after those checks pass. A dispatch on another branch is skipped.
4. Verify installation from PyPI and create release notes for that version.

The workflow runs only by manual dispatch. Pushing commits or creating tags does
not publish. PyPI versions cannot be overwritten; increment the version for a new
release. Adding this workflow does not itself create or upload a release.

Before the first release, confirm PyPI project ownership and metadata, review the
README rendering/links, and consider a TestPyPI rehearsal with a separately
configured publisher. The existing pending publisher is for production PyPI.
If you add a protected GitHub environment later, update both the publisher and
workflow to the same environment name.

See [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
for the OIDC identity requirements. `CODECOV_TOKEN` is used only by coverage CI
and provides no package-publishing access.
