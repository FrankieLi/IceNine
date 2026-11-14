<!-- table:cpp_split -->
| run | whole (s) | start-up and exit (s) | global searches (s) | local optimizations (s) |
|---|---|---|---|---|
| C0 (not gated: no preflight or quiet check) | 1247.03 | 14.0 | 1142 | 91 |
| C1 | 1236.87 | 14.9 | 1137 | 85 |
| C2 | 1231.87 | 13.9 | 1134 | 84 |
<!-- /table:cpp_split -->

<!-- table:cpu -->
| run | user / real | (user + sys) / real |
|---|---|---|
| Python A1 | 1.00 | 1.00 |
| Python A2 | 1.00 | 1.00 |
| C++ C1 | 1.77 | 2.00 |
| C++ C2 | 1.77 | 2.00 |
<!-- /table:cpu -->

<!-- table:paired -->
| wrong in both | wrong only in C++ | wrong only in Python | wrong in neither |
|---|---|---|---|
| 12 | 404 | 12 | 1572 |
<!-- /table:paired -->

<!-- table:quiet -->
| run | load samples | over 2.5 | over 3.5 (limit) | max load | median load | strict check passed |
|---|---|---|---|---|---|---|
| C0 | 43 | 20 | 10 | 33.56 | n/a | not checked |
| C1 | 43 | 1 | 0 | 2.52 | 2.16 | True |
| C2 | 43 | 7 | 3 | 6.22 | 2.14 | False |
<!-- /table:quiet -->

<!-- table:trace -->
| count in run C1 | voxels |
|---|---|
| local-optimization records in the client log | 1921 |
| voxels marked REFIT (hit ratio / best at most 0.9) | 393 |
| of those, wrong | 393 |
| voxels accepted (above 0.9) that are wrong | 23 |
| voxels with no record (seeds), of which wrong | 79 |
| wrong among the voxels with no record | 0 |
| wrong voxels in total | 416 |
| wrong voxels within 1 degree of another grain's truth | 416 |
| records whose Optimized orientation equals Original (printed digits) | 1921 |
| distinct orientations in the C++ output map | 80 |
<!-- /table:trace -->

<!-- table:unit -->
| side | seconds per global search | BFS time (s) | microseconds per evaluation |
|---|---|---|---|
| Python arm A1 | 20.9 | 1620 | 63.3 |
| C++ C1 | 13.4 | 1222 | n/a |
<!-- /table:unit -->

<!-- table:walls -->
| run | side | whole process (s) | label |
|---|---|---|---|
| A1 | Python | 1621.02 | single-worker; strict quiet check failed (busy: mediaanalysisd, searchpartyd) |
| A2 | Python | 1627.40 | single-worker; strict quiet check failed (load max 4.30 > 2.5 in 3 of 56 samples; no process over 20% CPU caught) |
| C1 | C++ one client | 1236.87 | 2 MPI ranks (1 fitting client + 1 server using about one core); strictly quiet at limit 3.5 only: 1 of 43 load samples are over 2.5 (max 2.52); the limit came from an assumed idle load of 1.5, not a measured one |
| C2 | C++ one client | 1231.87 | 2 MPI ranks (1 fitting client + 1 server using about one core); strict quiet check failed: load samples above 2.5 (7 of 43; above 3.5: 3) and no process over 20% CPU caught by the sampler (processes sampled every 300 s, load every 30 s); condition (ii) of N1 held only because nothing was caught; reported under R4b (D-20261114-N1) |
<!-- /table:walls -->

<!-- table:work -->
| side | global searches | local optimizations | evaluations | wrong voxels |
|---|---|---|---|---|
| Python arm A (item 14) | 69 | not counted here | 25602383 | 24 of 2000 (0.8% to 1.8%) |
| C++ one client (C1; C0 and C2 equal) | 85 | 1921 | not measured | 416 of 2000 (19.1% to 22.6%) |
<!-- /table:work -->
