cTreeBalls documentation

See README.rst in this directory for build instructions.
Maintained entry point: docs/index.rst.
Strict HTML check, from the source root:

    python3 -m sphinx -E -a -n -W --keep-going -b html docs docs/_build/html

Open docs/_build/html/index.html after building.
Edit .rst sources, not generated HTML, LaTeX, or the historical man pages.
The current numerical contract is in docs/3pcf.rst; engine selection is in
docs/search_methods.rst. Old documentation under addons/docs_v1.0.0 is archival.
