API Reference
=============

The supported user interfaces are the ``cballs`` executable and the
``cyballs.cballs`` Cython class.  C entry points are documented here for
contributors; they are not a separately versioned stable library API.

Command-Line Interface
----------------------

Executable
    ``./cballs``

Help
    ``./cballs --help``

Parameter file
    ``./cballs path/to/parameters_file``

See :doc:`user/command-line` for syntax and common parameters.

Python Module
-------------

.. py:module:: cyballs

.. py:class:: cballs(default=True)

   Wrapper around the compiled cTreeBalls C lifecycle.

   .. py:method:: set(*parameters, **kwargs)

      Update runtime parameters from one mapping or keyword arguments.

   .. py:method:: Run(level=["MainLoop"])

      Execute dependencies through the requested stage.  The default performs
      a full search and leaves result arrays available to getters.

   .. py:method:: clean()

      Clear Python-side parameters.

   .. py:method:: struct_cleanup()

      Release C-owned structures according to allocation flags.

   .. py:method:: clean_all()

      Release C state and clear Python parameters.

   .. py:method:: getrBins()

      Return the radial-bin array.

   .. py:method:: getHistNN()

      Return the neighbor-count histogram.

   .. py:method:: getHistXi2pcf()

      Return the two-point correlation function.

   .. py:method:: getHistZetaMsincos(m, type)

      Return one 3PCF multipole/component matrix.

Exceptions
~~~~~~~~~~

``cBallsError`` is the base wrapper exception.  ``cBallsSevereError`` reports
invalid parameters and setup failures; ``cBallsComputationError`` reports C
stage failures.

C Lifecycle
-----------

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Function
     - Responsibility
   * - ``StartRun``
     - Parse process arguments and initialize command/global state.
   * - ``StartRun_Common``
     - Shared validation, input loading, directory setup, and allocation.
   * - ``PrintParameterFile``
     - Write the used-values parameter record.
   * - ``SetNumberThreads``
     - Apply OpenMP thread control.
   * - ``MainLoop``
     - Build trees, dispatch the selected search, and evaluate histograms.
   * - ``EvalHist``
     - Select a concrete search method and write enabled products.
   * - ``MakeTree``
     - Construct and prepare the catalog tree.
   * - ``EndRun``
     - Close outputs and release process-owned state.

See :doc:`code_structure` for module ownership and extension points.
