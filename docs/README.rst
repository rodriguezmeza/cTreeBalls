Building the cTreeBalls Documentation
=====================================

The maintained Sphinx source is this directory. Start with ``index.rst``.
Historical copies under ``addons/docs_v1.0.0`` are not current documentation.

Use Python 3.11 or newer for the pinned documentation toolchain.
From the repository root::

   python3 -m pip install -r docs/requirements.txt
   python3 -m sphinx -E -a -n -W --keep-going -b html docs docs/_build/html

Open ``docs/_build/html/index.html``. No C executable, compiled cyballs
extension, or input catalog is required for the documentation build.

Other builders::

   make -C docs latex
   make -C docs latexpdf
   make -C docs linkcheck

``latexpdf`` additionally needs a TeX installation. ``linkcheck`` accesses
external sites and may report network restrictions independently of source
errors. The strict HTML build resolves internal references without networking.

The API reference is maintained from the Cython source rather than importing
a machine-specific extension. Update it whenever public methods change.
Keep the root README, active addon READMEs, and the Sphinx pages consistent.
Do not edit generated files or archived documentation.

Read the Docs uses ``.readthedocs.yaml`` and ``docs/requirements.txt``.
Changes become public only after the appropriate repository/deployment update.
