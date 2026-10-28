"""Build the current user manual from reStructuredText and Markdown sources."""

from pathlib import Path
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())["project"]
project = metadata["name"]
author = "Benjy Marks"
copyright = "2026, Benjy Marks"
release = version = metadata["version"]
extensions = ["sphinx.ext.autodoc", "sphinx.ext.mathjax", "myst_parser"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}
root_doc = "index"
exclude_patterns = ["build"]
html_theme = "alabaster"
html_title = "poly-mpm user guide"
html_theme_options = {"description": "2D material point method with grain-size segregation"}
html_sidebars = {"**": ["about.html", "navigation.html", "searchbox.html"]}
myst_heading_anchors = 3
