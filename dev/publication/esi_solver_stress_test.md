# ESI: Solver stress test for deeply temporalized technosphere systems

This supplementary experiment documents the computational reason for the
graph-matrix hybrid strategy used in `TRAILS`. The input was a Frictionless data
package generated with `premise`, extended with temporal technosphere exchanges,
and intended to be loaded by `TRAILS`. `TRAILS` interpolated the package to an
annual grid from 2004 to 2101, giving 98 year-specific technosphere systems of
size 41,792 x 41,792. The loaded three-dimensional technosphere array contained
48.7 million non-zero entries and occupied 1.27 GiB in sparse coordinate form.

We compared two mathematically comparable solution pathways. In the monolithic
pathway, annual technosphere matrices were stitched along a block diagonal and
temporal technosphere exchanges were inserted as off-diagonal links, following
the time-expanded matrix logic used by tools such as `bw_timex` (Diepers et al.,
2026). The resulting sparse matrix was solved with UMFPACK. In the `TRAILS`
pathway, the same temporal system was handled by routing temporal exchanges
through the graph and solving ordinary year-specific technosphere systems
sequentially. For the reduced
one-, two-, and three-year comparisons, `TRAILS` first performed the same
depth-four routing from the functional unit and then solved the routed frontier
demands that fell within the selected horizon. The benchmark used activity
index 31,387, demand 1, and start year 2025. Benchmarks were run on a MacBook Pro with an Apple M1 Pro chip, 8 CPU cores, 16 GB unified memory,
macOS 26.3, and arm64 architecture. The temporal links row reports only
cross-year entries inserted in the equivalent time-expanded matrix. It excludes
same-year remnants of temporal distributions and offsets outside the selected
horizon, so the one-year case has no temporal links in this table. In the
memory column, paired rows report LU numeric object size / peak factorization
memory; the full `TRAILS` row reports process RSS peak increase. For the 98-year
monolithic system, memory is a conservative extrapolated lower bound.

<table>
  <thead>
    <tr>
      <th rowspan="2">Metric</th>
      <th colspan="2">1 year</th>
      <th colspan="2">2 years</th>
      <th colspan="2">3 years</th>
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
      <td>2004-2101</td>
      <td>2004-2101</td>
    </tr>
    <tr>
      <td>Matrix / solve size</td>
      <td>41,792 x 41,792; 478,782 nnz</td>
      <td>41,792 x 41,792</td>
      <td>83,584 x 83,584; 993,357 nnz</td>
      <td>41,792 x 41,792</td>
      <td>125,376 x 125,376; 1,507,926 nnz</td>
      <td>41,792 x 41,792</td>
      <td>4,095,616 x 4,095,616; 57,300,688 nnz</td>
      <td>98 systems; 48,712,196 nnz in `TRAILS`.A</td>
    </tr>
    <tr>
      <td>Cross-year temporal entries</td>
      <td>0</td>
      <td>0</td>
      <td>17,313</td>
      <td>17,313</td>
      <td>34,620</td>
      <td>34,620</td>
      <td>11,004,530</td>
      <td>11,004,530</td>
    </tr>
    <tr>
      <td>Sparse matrix storage</td>
      <td>5.6 MiB</td>
      <td>5.9 MiB</td>
      <td>12 MiB</td>
      <td>5.9 MiB</td>
      <td>18 MiB</td>
      <td>5.9 MiB</td>
      <td>671 MiB</td>
      <td>1.3 GiB</td>
    </tr>
    <tr>
      <td>LU / memory footprint</td>
      <td>70 MiB / 313 MiB</td>
      <td>108 MiB / 493 MiB</td>
      <td>228 MiB / 1.5 GiB</td>
      <td>108 MiB / 493 MiB</td>
      <td>10 GiB / 13 GiB</td>
      <td>108 MiB / 493 MiB</td>
      <td>&ge;327 GiB / &ge;433 GiB</td>
      <td>1.1 GiB</td>
    </tr>
    <tr>
      <td>Routing time</td>
      <td>n.a.</td>
      <td>4.4 s; 50,576 nodes</td>
      <td>n.a.</td>
      <td>3.3 s; 50,576 nodes</td>
      <td>n.a.</td>
      <td>3.6 s; 50,576 nodes</td>
      <td>n.a.</td>
      <td>3.8 s; 50,576 nodes</td>
    </tr>
    <tr>
      <td>Solving time</td>
      <td>1.8 s</td>
      <td>3.4 s</td>
      <td>14 s</td>
      <td>6.6 s</td>
      <td>failed before solve</td>
      <td>10 s</td>
      <td>not attempted</td>
      <td>173 s</td>
    </tr>
    <tr>
      <td>Outcome</td>
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
not the limiting factor. Even the three-year monolithic matrix occupied only
18 MiB as a CSC matrix, yet its LU factorization required enough fill-in for
UMFPACK to estimate a 10 GiB numeric object and a 13 GiB peak memory
requirement, after which factorization failed. `TRAILS` avoids this fill-in by
never factorizing the global time-expanded system. Its memory footprint is
therefore governed by the largest ordinary annual technosphere solve and by the
temporal routing data structures, while the full 98-year deeply temporalized
case remains solvable on the test machine. A conservative linear extrapolation
from the three-year UMFPACK failure gives an approximate lower bound of 327 GiB
for the LU numeric object and 433 GiB peak memory for a 98-year monolithic
direct solve. This should be read as a lower bound rather than a prediction:
the observed LU fill-in already grew faster than sparse matrix storage between
the two- and three-year cases.

Diepers et al. (2026) claim that embedding time directly in one technosphere
and biosphere matrix is, by itself, a computational advantage. A unified
time-expanded matrix can be convenient for integration with existing matrix LCA
software and for representing flexible temporal resolution. However, for deep
temporalization, where temporal exchanges can occur throughout the background
supply chain, the global matrix becomes impractical before the full system is
represented. In this benchmark, the limiting factor is not the number of sparse
entries in the time-expanded matrix, but the fill-in and memory required by
direct factorization. `TRAILS` therefore treats the temporal graph and the annual
matrix solves as separate but coupled operations: temporal routing carries the
time dimension, while bounded year-specific solves retain the numerical
advantages of conventional matrix LCA.

Reference: Diepers et al., (2026). bw_timex: A Python Package for Time-Explicit Life Cycle Assessment. Journal of Open Source Software, 11(120), 9621, https://doi.org/10.21105/joss.09621
