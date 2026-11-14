<!-- table:cpp_split -->
| run | whole (s) | start-up and exit (s) | global searches (s) | local optimizations (s) | user / real | (user + sys) / real |
|---|---|---|---|---|---|---|
| C0 (untimed) | 1247.03 | 14.0 | 1142 | 91 | 1.77 | 2.00 |
| C1 | 1236.87 | 14.9 | 1137 | 85 | 1.77 | 2.00 |
| C2 | 1231.87 | 13.9 | 1134 | 84 | 1.77 | 2.00 |
<!-- /table:cpp_split -->

<!-- table:paired -->
| wrong in both | wrong only in C++ | wrong only in Python | wrong in neither |
|---|---|---|---|
| 12 | 404 | 12 | 1572 |
<!-- /table:paired -->

<!-- table:quiet -->
| run | load samples | over 2.5 | over 3.5 (limit) | max load | median load | strict check passed |
|---|---|---|---|---|---|---|
| C1 | 43 | 1 | 0 | 2.52 | 2.16 | True |
| C2 | 43 | 7 | 3 | 6.22 | 2.14 | False |
<!-- /table:quiet -->

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
| C1 | C++ one client | 1236.87 | 2 MPI ranks: 1 fitting client + 1 busy-waiting server; strictly quiet (limit 3.5) |
| C2 | C++ one client | 1231.87 | 2 MPI ranks: 1 fitting client + 1 busy-waiting server; strict quiet check failed (load over 3.5 in 3 of 43 samples, over 2.5 in 7; no process over 20% CPU caught by the sampler); reported under R4b |
<!-- /table:walls -->

<!-- table:work -->
| side | global searches | local optimizations | evaluations | wrong voxels |
|---|---|---|---|---|
| Python arm A (item 14) | 69 | not counted here | 25602383 | 24 of 2000 |
| C++ one client (C1; C0 and C2 equal) | 85 | 1921 | not measured | 416 of 2000 (19.1% to 22.6%) |
<!-- /table:work -->
