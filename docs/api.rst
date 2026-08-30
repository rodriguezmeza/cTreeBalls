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

   Python wrapper around the compiled cTreeBalls C lifecycle.

   If ``default`` is true, construction calls :py:meth:`set_default`.
   Construction also checks that Cython and the linked C library agree on the
   sizes of ``cmdline_data`` and ``global_data``.

   .. py:method:: set_default(**overrides)

      Reset the wrapper to the compiled C defaults by calling the C
      ``input_default_params`` path. Existing C allocations are released first.
      Keyword overrides are applied with :py:meth:`set`.

   .. py:method:: set(*parameters, **kwargs)

      Stage runtime parameters. Accepts either one mapping or keyword
      arguments::

         model.set({"nbody": 4096})
         model.set(nbody=4096, rootDir="Output_python")

      Unknown parameters are rejected when :py:meth:`Run` sends them through
      the C parser.

   .. py:method:: Run(level=["MainLoop"])

      Execute the requested C lifecycle stages and return CPU time per thread.

      Valid levels are ``input``, ``StartRun_Common``, ``PrintParameterFile``,
      ``SetNumberThreads``, ``Initial``, ``MainLoop``, and ``EndRun``.
      The default ``["MainLoop"]`` leaves histogram arrays alive for Python
      getters. Do not request ``EndRun`` before reading histogram getters.

   .. py:method:: clean()

      Clear Python-side staged parameters. This does not release C allocations.

   .. py:method:: struct_cleanup()

      Release C-owned allocations when a run is active.

   .. py:method:: clean_all()

      Release C-owned allocations and clear Python-side staged parameters.

   .. py:method:: abi_sizes()

      Return the linked C sizes of ``cmdline_data`` and ``global_data``.

Scalar Getters
~~~~~~~~~~~~~~

``getNThreads()``, ``getnMultipoles()``, ``getTheta()``, ``getrsmooth()``,
``getCPUTime()``, ``getsizeHistN()``, ``getVersion()``, ``getRootDir()``, and
``getNBody()`` expose selected C command/global state.

Allocation-state getters include ``getTreeAllocated()``, ``getAllocated2()``,
``getBodytableAllocated()``, ``getHistogramsAllocated()``, ``getGDAllocated()``,
and ``getCMDAllocated()``.

Histogram Getters
~~~~~~~~~~~~~~~~~

Histogram getters require live arrays. Call ``Run(level=["MainLoop"])`` first,
read the arrays, then call ``clean_all()``.

.. py:method:: cballs.getrBins()

   Return radial-bin values as a NumPy array.

.. py:method:: cballs.getHistNN()

   Return the neighbor-count histogram.

.. py:method:: cballs.getHistCF()

   Return the count-derived correlation histogram.

.. py:method:: cballs.getHistXi2pcf()

   Return the 2PCF histogram.

.. py:method:: cballs.getHistXi2pcf12()
.. py:method:: cballs.getHistXi2pcf13()

   Return cross-catalog 2PCF histograms when available in the active build.

.. py:method:: cballs.getHistZetaMsincos(m, type)

   Return one 3PCF multipole matrix. ``type`` values are:
   ``1=cos``, ``2=sin``, ``3=sincos``, ``4=cossin``.

.. py:method:: cballs.getHistZetaM_EE(m)

   Return the edge-corrected 3PCF multipole matrix when enabled.

Shear Getters
~~~~~~~~~~~~~

The following getters require an ``OCTREESHEAROMP`` build and a successful
``octree-shear-omp`` run:

``getShearXiPlus()``, ``getShearXiMinus()``, and ``getShearXiWeight()``
    Return complex 2PCFs and their real pair-weight denominator, each with
    shape ``(B,)``.

``getShearUpsilonXMultipoles()`` and ``getShearGammaXMultipoles()``
    Return raw and edge-corrected Porth x-projection natural-component
    multipoles with shape ``(4, 2*nmax+1, B, B)``.

``getShearUpsilonMultipoles()`` and ``getShearGammaMultipoles()``
    Backward-compatible aliases for the projection-explicit getters above.

``getShearWindowMultipoles()``
    Return ``N_n`` with shape ``(4*nmax+1, B, B)``.

``getShearGammaX()``
    Return angular Porth x-projection natural components with shape
    ``(4, P, B, B)``.

``getShearGamma()``
    Backward-compatible alias for ``getShearGammaX()``.

``getShearMultipoleOrders()`` and ``getShearPhiCenters()``
    Return the coordinate axes for the multipole and angular products.

Exceptions
~~~~~~~~~~

``CosmoError``
    Base wrapper exception.

``CosmoSevereError``
    Invalid input, ABI mismatch, parser rejection, allocation failure, or unsafe
    getter call.

``CosmoComputationError``
    Failure reported by one of the C computation stages.

``CosmoSevereErrorDummy``
    Internal wrapper used when a lower-level C getter reports failure.

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
