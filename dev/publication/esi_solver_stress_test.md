# ESI: Solver stress test for deeply temporalized technosphere systems

This supplementary experiment documents the computational reason for the
graph-matrix hybrid strategy used in `TRAILS`. The input was a Frictionless data
package generated with `premise`, extended with temporal technosphere exchanges,
and intended to be loaded by `TRAILS`. `TRAILS` interpolated the package to an
annual grid from 2004 to 2101, giving 98 year-specific technosphere systems of
size 42,300 x 42,300. The loaded three-dimensional technosphere array contained
49.2 million non-zero entries and occupied 1.28 GiB in sparse coordinate form.

The benchmark runner is `dev/publication/solver_stress_test.py`. It uses the
publication-local data package `dev/publication/trails_remind_SSP2-PkBudg1000.zip`
and writes `dev/publication/monolithic_vs_trails_remind_SSP2-PkBudg1000.csv`.
This path layout matters for reproducibility. During repository restructuring,
the old helper script location and absolute defaults under `dev/` became stale;
the runner should stay next to the large publication ZIPs or receive explicit
`--package` and `--output` paths. The annual interpolation cache used in this
run was outside the repository at
`~/Library/Application Support/trails/cache/interp_0536f3a1c684`; a clean
machine can rebuild it, but load-time memory and timing will differ. The ZIP
inputs are large publication artifacts and are not part of the Python package
itself, so they should either remain documented as local artifacts or be moved
to a durable external data location before final archiving.

We compared two mathematically comparable solution pathways. In the monolithic
pathway, annual technosphere matrices were stitched along a block diagonal and
temporal technosphere exchanges were inserted as time-expanded activity-year
links, following the matrix logic used by tools such as `bw_timex` (Diepers et
al., 2026). The resulting sparse matrix was solved with UMFPACK. In the
`TRAILS` pathway, the same temporal system was handled by routing temporal
exchanges through the graph and solving ordinary year-specific technosphere
systems sequentially. For the reduced one-, two-, three-, and four-year
comparisons, `TRAILS` first performed the same depth-four routing from the
functional unit and then solved the routed frontier demands that fell within
the selected horizon. The benchmark used activity index 31,387, demand 1, and
start year 2025. Benchmarks were run on a MacBook Pro with an Apple M1 Pro chip,
8 CPU cores, 16 GB unified memory, macOS 26.3, and arm64 architecture.

The temporal-entry row reports temporal-distribution entries retained inside
the equivalent time-expanded matrix horizon. It therefore includes same-year
remainders as well as cross-year links; the one-year case is not zero. In the
memory row, paired rows report LU numeric object size / peak factorization
memory. The full `TRAILS` row reports the additional process RSS peak during
the full routing and LCA step after the package had already been loaded. For
the 98-year monolithic system, memory is a conservative extrapolated lower
bound from the four-year UMFPACK failure.

<table>
  <thead>
    <tr>
      <th rowspan="2">Metric</th>
      <th colspan="2">1 year</th>
      <th colspan="2">2 years</th>
      <th colspan="2">3 years</th>
      <th colspan="2">4 years</th>
      <th colspan="2">98 years</th>
    </tr>
    <tr>
      <th>Monolithic LU</th>
      <th>`TRAILS`</th>
      <th>Monolithic LU</th>
      <th>`TRAILS`</th>
      <th>Monolithic LU</th>
      <th>`TRAILS`</th>
      <th>Monolithic LU</th>
      <th>`TRAILS`</th>
      <th>Monolithic LU</th>
      <th>`TRAILS`</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>Years</td>
      <td>2025</td>
      <td>2025</td>
      <td>2024-2025</td>
      <td>2024-2025</td>
      <td>2023-2025</td>
      <td>2023-2025</td>
      <td>2022-2025</td>
      <td>2022-2025</td>
      <td>2004-2101</td>
      <td>2004-2101</td>
    </tr>
    <tr>
      <td>Matrix / solve size</td>
      <td>42,300 x 42,300; 482,387 nnz</td>
      <td>42,300 x 42,300</td>
      <td>84,600 x 84,600; 997,882 nnz</td>
      <td>42,300 x 42,300</td>
      <td>126,900 x 126,900; 1,517,500 nnz</td>
      <td>42,300 x 42,300</td>
      <td>169,200 x 169,200; 2,037,270 nnz</td>
      <td>42,300 x 42,300</td>
      <td>4,145,400 x 4,145,400; approximately 59.3 million nnz</td>
      <td>98 systems; 49,162,344 nnz in `TRAILS`.A</td>
    </tr>
    <tr>
      <td>Temporal entries retained</td>
      <td>1,144</td>
      <td>1,144</td>
      <td>14,978</td>
      <td>14,978</td>
      <td>32,935</td>
      <td>32,935</td>
      <td>51,044</td>
      <td>51,044</td>
      <td>10,440,416</td>
      <td>10,440,416</td>
    </tr>
    <tr>
      <td>Sparse matrix storage</td>
      <td>5.7 MiB</td>
      <td>5.9 MiB</td>
      <td>11.7 MiB</td>
      <td>5.9 MiB</td>
      <td>17.9 MiB</td>
      <td>5.9 MiB</td>
      <td>24.0 MiB</td>
      <td>5.9 MiB</td>
      <td>approximately 694 MiB</td>
      <td>1.28 GiB</td>
    </tr>
    <tr>
      <td>LU / memory footprint</td>
      <td>62 MiB / 223 MiB</td>
      <td>88 MiB / 458 MiB</td>
      <td>206 MiB / 732 MiB</td>
      <td>99 MiB / 477 MiB</td>
      <td>383 MiB / 1.8 GiB</td>
      <td>99 MiB / 478 MiB</td>
      <td>15.3 GiB / 19.5 GiB</td>
      <td>99 MiB / 478 MiB</td>
      <td>&ge;375 GiB / &ge;478 GiB</td>
      <td>615 MiB RSS peak increase after load</td>
    </tr>
    <tr>
      <td>Routing time</td>
      <td>n.a.</td>
      <td>0.98 s; 6,214 nodes</td>
      <td>n.a.</td>
      <td>0.18 s; 6,214 nodes</td>
      <td>n.a.</td>
      <td>0.25 s; 6,214 nodes</td>
      <td>n.a.</td>
      <td>0.17 s; 6,214 nodes</td>
      <td>n.a.</td>
      <td>0.21 s; 6,214 nodes</td>
    </tr>
    <tr>
      <td>Solving time</td>
      <td>1.4 s</td>
      <td>3.4 s</td>
      <td>10.7 s</td>
      <td>6.5 s</td>
      <td>29.7 s</td>
      <td>9.6 s</td>
      <td>failed before solve</td>
      <td>11.8 s</td>
      <td>not attempted</td>
      <td>72.1 s</td>
    </tr>
    <tr>
      <td>Outcome</td>
      <td>solved</td>
      <td>solved</td>
      <td>solved</td>
      <td>solved</td>
      <td>solved</td>
      <td>solved</td>
      <td>failed</td>
      <td>solved</td>
      <td>estimated impractical</td>
      <td>solved</td>
    </tr>
  </tbody>
</table>

The important observation is that sparse storage of the time-expanded matrix is
not the limiting factor. The four-year monolithic matrix occupied only 24 MiB
as a CSC matrix, yet its LU factorization required enough fill-in for UMFPACK
to estimate a 15.3 GiB numeric object and a 19.5 GiB peak memory requirement,
after which factorization failed. `TRAILS` avoids this fill-in by never
factorizing the global time-expanded system. Its memory footprint is therefore
governed by the largest ordinary annual technosphere solve, the temporal
routing data structures, and the loaded annual matrix stack, while the full
98-year deeply temporalized case remains solvable on the test machine.

A conservative linear extrapolation from the four-year UMFPACK failure gives an
approximate lower bound of 375 GiB for the LU numeric object and 478 GiB peak
memory for a 98-year monolithic direct solve. This should be read as a lower
bound rather than a prediction: the observed LU fill-in grew much faster than
sparse matrix storage between the one-, two-, three-, and four-year cases.

Diepers et al. (2026) claim that embedding time directly in one technosphere
and biosphere matrix is, by itself, a computational advantage. A unified
time-expanded matrix can be convenient for integration with existing matrix LCA
software and for representing flexible temporal resolution. However, for deep
temporalization, where temporal exchanges can occur throughout the background
supply chain, the global matrix becomes impractical before the full system is
represented. In this benchmark, the limiting factor is not the number of sparse
entries in the time-expanded matrix, but the fill-in and memory required by
direct factorization. `TRAILS` therefore treats the temporal graph and the
annual matrix solves as separate but coupled operations: temporal routing
carries the time dimension, while bounded year-specific solves retain the
numerical advantages of conventional matrix LCA.

Reference: Diepers et al., (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621, https://doi.org/10.21105/joss.09621
