
Transform gadget multi-columns-ascii halos to columns-ascii

$ time cballs search=tree-omp-sincos rootDir=Output rangeN=150 rminHist=0.1 sizeHistN=30 nthreads=16 stepState=1000 computeTPCF=false useLogHist=false usePeriodic=true verb=2 in=halos_0.0.ascii infmt=multi-columns-ascii columns=9,10,11 options=stop,only-pos o=halos ofmt=columns-ascii


