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
| map after the pass against the script's: voxels that differ / by more than 0.001 degree / by more than 1 degree | 1991 / 1683 / 8 | - | same tries (every field) and the same changed voxels |
| control, key off: LibS11 / stored S11, pass tries | 49,968 | 49,970 | - |
| control, key off: LibS11 / stored S11, pass evaluations | 3,048,048 | 3,048,170 | - |
<!-- /table:libpass -->
