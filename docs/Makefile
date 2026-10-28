# Run after installing the optional docs dependencies; see docs/README.md.
SPHINXBUILD ?= sphinx-build
SPHINXOPTS ?= -W --keep-going
BUILDDIR ?= build

.PHONY: html linkcheck
html:
	$(SPHINXBUILD) -b html $(SPHINXOPTS) . $(BUILDDIR)/html

linkcheck:
	$(SPHINXBUILD) -b linkcheck $(SPHINXOPTS) . $(BUILDDIR)/linkcheck
