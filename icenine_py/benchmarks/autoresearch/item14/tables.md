<!-- table:e0_arms -->
| arm | wrong | strict evaluations per voxel | smooth evaluations per voxel |
|---|---|---|---|
| off | 41 (15.5% to 26.6%) | 48700 | 0 |
| F | 41 (15.5% to 26.6%) | 48698 | 0 |
| S | 7 (1.7% to 7.0%) | 49215 | 252 |
| FS | 7 (1.7% to 7.0%) | 49215 | 252 |
<!-- /table:e0_arms -->

<!-- table:e0_pairs -->
| pair | fixed | broken |
|---|---|---|
| F vs off | 0 | 0 |
| FS vs S | 0 | 0 |
| FS vs F | 34 | 0 |
| S vs off | 34 | 0 |
<!-- /table:e0_pairs -->

<!-- table:pairs_time -->
| arm | first run, whole process (s) | repeat (s) | difference (%) | median BFS wall_total_s (s) | timing label |
|---|---|---|---|---|---|
| A | 1621.02 | 1627.40 | 0.39 | 1623.0 | single-worker; strict quiet check failed (busy: mediaanalysisd, searchpartyd); strict quiet check failed (load max 4.30 > 2.5 in 3 of 56 samples; no process over 20% CPU caught) |
| B | 4484.94 | 4481.69 | 0.07 | 4482.1 | single-worker; strict quiet check failed (load max 4.46 > 2.5 in 4 of 151 samples; no process over 20% CPU caught); strict quiet check failed (load max 4.69 > 2.5 in 2 of 151 samples; busy: mediaanalysisd, photolibraryd) |
<!-- /table:pairs_time -->

<!-- table:pilot_arms -->
| arm | wrong | full searches | strict evaluations | of these, full searches | of these, local optimization | smooth evaluations |
|---|---|---|---|---|---|---|
| off | 37 (1.3% to 2.5%) | 69 | 25,547,935 | 22,882,865 | 2,665,070 | 0 |
| F | 24 (0.8% to 1.8%) | 69 | 25,602,383 | 22,904,126 | 2,698,257 | 0 |
| S | 34 (1.2% to 2.4%) | 65 | 24,490,331 | 21,810,503 | 2,679,828 | 980,594 |
| FS | 26 (0.9% to 1.9%) | 69 | 25,845,740 | 23,156,710 | 2,689,030 | 1,042,020 |
<!-- /table:pilot_arms -->

<!-- table:pilot_pairs -->
| pair | fixed | broken |
|---|---|---|
| FS vs S | 14 | 6 |
| FS vs F | 6 | 8 |
| F vs off | 20 | 7 |
| S vs off | 14 | 11 |
<!-- /table:pilot_pairs -->

<!-- table:quiet_runs -->
| run | max load | median load | load samples | preflight load | busy processes caught |
|---|---|---|---|---|---|
| A1 | 1.77 | 1.33 | 56 | 1.08 | searchpartyd, mediaanalysisd |
| A2 | 4.30 | 1.32 | 56 | 1.07 | none caught |
| B1 | 4.46 | 1.26 | 151 | 1.26 | none caught |
| B2 | 4.69 | 1.30 | 151 | 1.13 | mediaanalysisd, photolibraryd |
<!-- /table:quiet_runs -->

<!-- table:split -->
| voxel | arm | wall (s) | screen (%) | smooth (%) | strict (%) | rest (%) | strict evaluations | smooth evaluations | timing label |
|---|---|---|---|---|---|---|---|---|---|
| 15901 | F | 21.0 | 11.2 | 0.0 | 48.2 | 40.6 | 328,182 | 0 | single-worker; quiet |
| 15901 | FS | 45.7 | 5.2 | 54.0 | 22.2 | 18.6 | 329,256 | 14,816 | single-worker; quiet |
| 17034 | F | 21.9 | 10.9 | 0.0 | 48.5 | 40.5 | 340,242 | 0 | single-worker; quiet |
| 17034 | FS | 46.2 | 5.2 | 52.4 | 23.2 | 19.3 | 341,081 | 15,199 | single-worker; quiet |
<!-- /table:split -->

<!-- table:timing_runs -->
| run | arm | whole process (s) | BFS wall_total_s | full searches (s) | local optimization (s) | label |
|---|---|---|---|---|---|---|
| A1 | A (F) | 1621.0 | 1619.8 | 1445.2 | 174.2 | strict quiet check failed (busy: mediaanalysisd, searchpartyd) |
| A2 | A (F) | 1627.4 | 1626.2 | 1450.4 | 175.4 | strict quiet check failed (load max 4.30 > 2.5 in 3 of 56 samples; no process over 20% CPU caught) |
| B1 | B (FS) | 4484.9 | 4483.7 | 4308.0 | 175.3 | strict quiet check failed (load max 4.46 > 2.5 in 4 of 151 samples; no process over 20% CPU caught) |
| B2 | B (FS) | 4481.7 | 4480.4 | 4304.7 | 175.3 | strict quiet check failed (load max 4.69 > 2.5 in 2 of 151 samples; busy: mediaanalysisd, photolibraryd) |
<!-- /table:timing_runs -->
