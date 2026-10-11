# Timing examples supplied during the scientific discussion

These shell timings were supplied by the user. They are historical examples,
not new benchmarks performed by the tutor. The runs used the DES/Takahashi
shear workflow, realization 2, tomographic bin 1, region 1, at most 5,000,000
objects, the spherical ball-tree OpenMP engine, 3PCF, local-east-north shear,
Sofia figure-1 bins, seven retained multipoles, 32 phi bins, tree theta 1,
no smoothing, nsmooth 1, pivot reuse, phase tolerance 3 and 16 threads.

| Bin allowance | Profile | Wall seconds | User CPU seconds | System CPU seconds |
| --- | --- | ---: | ---: | ---: |
| 0.1 | enabled | 535.365 | 7258.336 | 22.222 |
| 0.1 | disabled | 533.647 | 7185.865 | 26.777 |
| 0.2 | disabled | 337.904 | 4510.271 | 15.913 |

The two profile-disabled runs differ by 1.5793 times in wall time, or 36.6802%
less elapsed time for bin allowance 0.2. This is not a factor of two. Larger bin
allowance can accept more grouped interactions and reduce descent, but these
timings alone contain no measurement of the approximation error. Compare the
saved observables against a matched reference to obtain numerical errors.

For the profile-enabled example, (user + system) / real is about 13.60 CPU
cores used on average. CPU time adds work across threads; 121 CPU minutes
therefore need not mean 121 minutes of elapsed time. The separate profile
on/off samples do not isolate profiling overhead from normal run variation.

The later measured approximation comparisons and their reference definitions
are recorded in `docs/APPROXIMATION_DEFAULTS.md`. Those measurements are also
specific to their catalogue samples and estimator settings.
