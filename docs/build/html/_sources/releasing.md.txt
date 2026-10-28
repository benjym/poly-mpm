# Packaging and future PyPI releases

The project builds a pure-Python wheel and source distribution. A GitHub Actions
job builds both, checks their metadata/README with Twine, and installs the wheel
in an isolated environment for a bounded solver run. These checks do not upload
a release. The repository version is currently `1.0.0`; confirm that version and
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
CI performs a short wheel run; manually verify the sdist installation as well.

Before publishing:

- Confirm access to the `poly-mpm` project on PyPI and TestPyPI. A local project
  name does not reserve it. Review version, author/contact, license, package
  contents, Python requirement, and dependency lower bounds.
- Review README links as rendered on PyPI; repository-relative guide links may
  need absolute URLs for the published description.
- Test installation from TestPyPI in a fresh environment, resolving numerical
  dependencies from PyPI. Check the console entry point and a bounded simulation.
- Configure [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
  for the exact repository, release workflow, and protected GitHub environment.
  Use an approved release/tag workflow with `id-token: write` on the publishing
  job and an environment approval if desired. Keep upload permissions out of CI.
- Publish the already-tested artifacts, verify installation from PyPI, and
  attach release notes describing numerical changes and compatibility settings.

The existing workflow has read-only repository permissions. `CODECOV_TOKEN`
is used only to upload coverage reports; it provides no package-publishing access.
No PyPI upload workflow or credentials are configured by this change.
