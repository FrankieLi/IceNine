<!-- table:libcheck -->
Reading line: one row per pair of runs on `mid`, before the boundary pass; the library run is compared with the script run it must equal. A counter is compared when the script run has it; the script counts its rescue evaluations inside the grid fallback's count, so that count is compared with the library's grid plus rescue evaluations.

| pair | map before the pass: voxels that differ | evaluations before the pass | counters (of the script run) equal | records equal | rescues: tried, rescued, tries | rescues equal | equal |
|---|---|---|---|---|---|---|---|
| LibG vs Gl | 0 | 10,064,185 (library), 10,064,185 (script) | yes | yes | 13, 10, 30 (library); 13, 10, 30 (script) | yes | yes |
| LibS11 vs S11 (control) | 0 | 13,034,526 (library), 13,034,526 (script) | yes | yes | - | - | yes |
<!-- /table:libcheck -->

<!-- table:libpass -->
Reading line: the boundary pass, by where it was made. The library run makes it inside the run, from the map in memory (unrounded angles). The script makes it after the run, from the stored map (a .mic file stores rounded angles), with the library's single generator. The last column is the library's pass function run on that same stored map.

| | library run (LibG): pass in the run | script (Gl): pass after the run | library pass function on the stored map of Gl |
|---|---|---|---|
| pass tries | 50,026 | 50,027 | 50,027 |
| pass evaluations | 3,051,586 | 3,051,647 | 3,051,647 |
| total with the pass | 13,115,771 | 13,115,832 | - |
| map after the pass against the script's: voxels more than 0.001 degree apart / more than 1 degree apart / largest angle | 1,683 / 8 / 60.2 degrees | - | same tries (every field) and the same changed voxels |
| map after the pass against the script's: voxels that differ in any bit (also storage differences below 0.001 degree) | 1,991 | - | - |
| wrong voxels after the pass (more than 1 degree from the truth) | 37 | 37 | - |
| by voxel, script pass as the base: fixed / broken / wrong in both | 4 / 4 / 33 | - | - |
| median error after the pass, all 24,570 voxels (degrees) | 0.02913 | 0.02913 | - |
| control, key off: LibS11 / stored S11, pass tries | 49,968 | 49,970 | - |
| control, key off: LibS11 / stored S11, pass evaluations | 3,048,048 | 3,048,170 | - |
<!-- /table:libpass -->

<!-- table:libtries -->
Reading line: the in-run pass of the library against the script's pass, try by try, one row per pair. Both passes try the same voxels and candidates in the same order at first. Columns 3 to 7 count, over those tries, where the two differ (ends compared as float32, the script's stored type; the angle between the two ends of a try, without symmetry). The last column counts, over the whole pass, the voxels that got a different number of tries. Same input maps would give an equal first try (try 0); a shifted generator would instead keep the early tries equal.

| pair | tries in the same (voxel, candidate) order from the start | of these: ends differ | costs differ | evaluation counts differ | first try whose end differs | angle between the two ends: median / largest (degrees) | voxels with another number of tries |
|---|---|---|---|---|---|---|---|
| LibG vs Gl | 11,632 | 11,626 | 3,030 | 0 | 0 | 0.0057 / 1.29 | 1 |
| LibS11 vs S11 (control) | 4,872 | 4,870 | 1,307 | 0 | 0 | 0.0057 / 0.74 | 4 |
<!-- /table:libtries -->
