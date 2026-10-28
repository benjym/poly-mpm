# Building and publishing packages

The project builds a pure-Python wheel and source distribution. CI checks
formatting, tests, documentation, metadata, and installed-wheel execution.
The publishing workflow also installs and runs the sdist before uploading.
Version tags trigger publication; manual publication from `main` is available.

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
Publishing. The configured publisher uses `benjym/poly-mpm`, workflow
`pypi.yml`, and **no GitHub environment**. The publishing job deliberately has
no `environment` field so its OIDC identity matches that configuration.
Only that job receives `id-token: write`; it downloads the tested distributions
without checking out or executing project source. No API token is required.

## Publish a version tag

1. Set the intended version in `pyproject.toml`, run `uv lock`, and commit/push
   the tested source to `main`. Inspect metadata, package contents and README
   links before tagging. Every published version must be unique.
2. Create an annotated tag matching that version and push it. For example, after
   setting the project version to `1.0.1`:

   ```sh
   git tag -a v1.0.1 -m "Release poly-mpm 1.0.1"
   git push origin refs/tags/v1.0.1
   ```

3. The tag push starts **Publish to PyPI**. Tags matching `v[0-9]*` trigger the
   workflow, and its version check requires the exact name `v<project version>`.
   The checked-out source, tests, distributions and published package all come
   from that tagged commit. A mismatched tag fails before publication.
4. After success, verify installation from PyPI and add release notes for the tag.

Pushing ordinary commits to `main` runs CI without publishing. Creating a GitHub
release for an already-pushed tag does not start another upload. PyPI versions
cannot be overwritten: increment the version and create a new tag for a new
release.

## Manual publication

Open **Actions → Publish to PyPI → Run workflow**, select `main`, and enter the
exact version from `pyproject.toml`. The workflow runs the same checks before
uploading. Manual dispatch on another branch is skipped. Use this option for an
intentional release without a new tag; do not dispatch an already-published
version. A pending publisher becomes active after its first successful upload.

The publisher uses production PyPI. A TestPyPI rehearsal requires a separately
configured publisher and destination. If you add a protected GitHub environment
later, update both the publisher and workflow to the same environment name.

See [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/using-a-publisher/)
for the OIDC identity requirements. `CODECOV_TOKEN` is used only by coverage CI
and provides no package-publishing access.
