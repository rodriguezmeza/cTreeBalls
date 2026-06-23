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

Notebook Flow
-------------

The notebook follows this sequence:

1. install native build tools and ``cTreeBalls==1.0.1`` from PyPI;
2. import ``cyballs.cballs``;
3. run a small generated-catalog calculation;
4. copy ``rBins``, ``histXi2pcf``, and ``histNN`` into NumPy arrays;
5. plot the 2PCF-style output and histogram counts;
6. optionally download the generated output directory.

This mirrors the :doc:`python-wrapper` tutorial but keeps everything inside the
Colab runtime.
