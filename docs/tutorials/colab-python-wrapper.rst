Google Colab Python Wrapper Example
===================================

An independent Google Colab notebook is available for running the Python
wrapper without cloning or building the repository manually:

`Open the cTreeBalls Colab notebook <https://colab.research.google.com/github/rodriguezmeza/cTreeBalls/blob/main/examples/cTreeBalls_minimal_colab.ipynb>`_

The notebook installs the published PyPI package, imports ``cyballs``, runs a
compact synthetic-catalog calculation based on :doc:`python-wrapper`, copies the
returned arrays into NumPy, and plots the radial bins, 2PCF-style values, and
histogram counts.

It is intended as a smoke test and starting point.  For production catalog
workflows, move on to :doc:`catalog-workflow` after the Colab example runs
successfully.
