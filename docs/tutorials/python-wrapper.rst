Python Wrapper Tutorial
=======================

This tutorial runs a reduced synthetic calculation and copies the 2PCF into
NumPy before releasing C-owned memory.

1. Build and Import
-------------------

.. code-block:: bash

   make clean
   make PYTHON=python3 all
   python3 -c "from cyballs import cballs; print(cballs)"

2. Create a Script
------------------

Save the following as ``run_cyballs.py`` in the repository root:

.. code-block:: python

   import numpy as np
   from cyballs import cballs

   model = cballs()
   model.set(
       nbody=4096,
       sizeHistN=12,
       mChebyshev=3,
       rootDir="Output_python_tutorial",
       numberThreads=1,
       verbose=0,
       verbose_log=0,
   )

   cpu_time = model.Run()
   radius = np.array(model.getrBins(), copy=True)
   xi = np.array(model.getHistXi2pcf(), copy=True)
   model.clean_all()

   print(f"CPU time: {cpu_time:.3f} s")
   print(np.column_stack([radius, xi]))

3. Run and Check
----------------

.. code-block:: bash

   python3 run_cyballs.py
   cat Output_python_tutorial/cyballs_param.txt-usedvalues

The parameter filename may differ slightly with wrapper/build revisions; use
the generated ``*-usedvalues`` file as the authoritative record.

4. Continue to Catalog Data
---------------------------

Replace the generated-model parameters with ``infile``, ``infileformat``, and
``iCatalogs`` only after the compact run works.  A HEALPix example is provided
in ``tests/python/kappa_corr.py`` and summarized in
:doc:`catalog-workflow`.
