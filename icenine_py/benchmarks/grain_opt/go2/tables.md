## go2_As50

args: `{"name": "go2_As50", "pipeline": "Astrict/50", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0141 | 0.0446 | 28 of 1082 (1.8 to 3.7%) | 6 of 1082 (0.3 to 1.2%) | 0.0838 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0148 | 0.0432 | 20 of 853 (1.5 to 3.6%) | 3 of 853 (0.1 to 1.0%) | 0.0829 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0125 | 0.0291 | 2 of 120 (0.5 to 5.9%) | 0 of 120 (0.0 to 3.1%) | 0.0769 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0155 | 0.0439 | 18 of 733 (1.6 to 3.8%) | 3 of 733 (0.1 to 1.2%) | 0.0838 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0186 | 0.0378 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0868 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0236 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0698 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0191 | 0.0478 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0848 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0107 | 0.1879 | 19 of 77 (16.4 to 35.4%) | 2 of 77 (0.7 to 9.0%) | 0.0980 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/50; accepted per stage [16]; renders per voxel 63.0; stage wall time 2.0 s (1.8 ms per voxel); timing: contended: timing preflight refused (duetexpertd at 97% CPU, Python at 100% CPU); single run, one worker

## go2_As100

args: `{"name": "go2_As100", "pipeline": "Astrict/100", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0143 | 0.0367 | 8 of 1082 (0.4 to 1.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0837 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0146 | 0.0350 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0819 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0120 | 0.0279 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0757 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0151 | 0.0364 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0829 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0186 | 0.0378 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0868 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0236 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0698 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0183 | 0.0446 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0850 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0095 | 0.0352 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0869 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/100; accepted per stage [17]; renders per voxel 113.0; stage wall time 3.3 s (3.1 ms per voxel); timing: contended: timing preflight refused (Python at 100% CPU); single run, one worker

## go2_As200

args: `{"name": "go2_As200", "pipeline": "Astrict/200", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0112 | 0.0348 | 8 of 1082 (0.4 to 1.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0818 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0111 | 0.0325 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0802 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0084 | 0.0199 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0756 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0115 | 0.0337 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0810 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0073 | 0.0183 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0807 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0159 | 0.0379 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0855 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0276 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0695 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0141 | 0.0322 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0833 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0093 | 0.0335 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0859 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/200; accepted per stage [17]; renders per voxel 213.0; stage wall time 5.7 s (5.2 ms per voxel); timing: contended: timing preflight refused (Python at 100% CPU); single run, one worker

## go2_As400

args: `{"name": "go2_As400", "pipeline": "Astrict/400", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0062 | 0.0222 | 6 of 1082 (0.3 to 1.2%) | 6 of 1082 (0.3 to 1.2%) | 0.0763 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0062 | 0.0218 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0757 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0049 | 0.0157 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0734 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0063 | 0.0225 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0761 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0082 | 0.0245 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0780 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0060 | 0.0221 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0808 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0056 | 0.0176 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0662 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0070 | 0.0172 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0772 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0031 | 0.0184 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0773 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/400; accepted per stage [20]; renders per voxel 413.0; stage wall time 12.0 s (11.1 ms per voxel); timing: contended: timing preflight refused (corespotlightd at 227% CPU, Python at 46% CPU, Slack at 33% CPU); single run, one worker

## go2_As200o2

args: `{"name": "go2_As200o2", "pipeline": "Astrict/200/o2", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0165 | 0.0504 | 28 of 1082 (1.8 to 3.7%) | 6 of 1082 (0.3 to 1.2%) | 0.0859 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0170 | 0.0491 | 20 of 853 (1.5 to 3.6%) | 3 of 853 (0.1 to 1.0%) | 0.0847 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0192 | 0.0361 | 2 of 120 (0.5 to 5.9%) | 0 of 120 (0.0 to 3.1%) | 0.0794 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0167 | 0.0493 | 18 of 733 (1.6 to 3.8%) | 3 of 733 (0.1 to 1.2%) | 0.0856 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0211 | 0.0402 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0895 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0126 | 0.0303 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0728 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0198 | 0.0516 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0858 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0156 | 0.1879 | 19 of 77 (16.4 to 35.4%) | 2 of 77 (0.7 to 9.0%) | 0.1010 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/200/o2; accepted per stage [17]; renders per voxel 207.0; stage wall time 6.2 s (5.7 ms per voxel); timing: contended: timing preflight refused (Python at 99% CPU); single run, one worker

## go2_As100_V10

args: `{"name": "go2_As100_V10", "pipeline": "Astrict/100,V/10", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0143 | 0.0367 | 8 of 1082 (0.4 to 1.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0837 |
| all | s1 | 1082 | 0.0145 | 0.0344 | 7 of 1082 (0.3 to 1.3%) | 6 of 1082 (0.3 to 1.2%) | 0.0780 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0146 | 0.0350 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0819 |
| voxels with a C fit | s1 | 853 | 0.0151 | 0.0338 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0771 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0120 | 0.0279 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0757 |
| interior (C voxels) | s1 | 120 | 0.0140 | 0.0290 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0737 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0151 | 0.0364 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0829 |
| boundary (C voxels) | s1 | 733 | 0.0152 | 0.0346 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0777 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s1 | 209 | 0.0230 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0795 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0186 | 0.0378 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0868 |
| level 0.1 (C voxels) | s1 | 222 | 0.0164 | 0.0373 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0820 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0236 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0698 |
| level 0.2 (C voxels) | s1 | 201 | 0.0107 | 0.0316 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0668 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0183 | 0.0446 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0850 |
| level 0.5 (C voxels) | s1 | 144 | 0.0169 | 0.0366 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0788 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0095 | 0.0352 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0869 |
| level 1.0 (C voxels) | s1 | 77 | 0.0099 | 0.0463 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0804 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/100,V/10; accepted per stage [17, 20]; renders per voxel 128.0; stage wall time 6.1 s (5.6 ms per voxel); timing: contended: timing preflight refused (Python at 99% CPU); single run, one worker

## go2_As100_V30

args: `{"name": "go2_As100_V30", "pipeline": "Astrict/100,V/30", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0143 | 0.0367 | 8 of 1082 (0.4 to 1.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0837 |
| all | s1 | 1082 | 0.0124 | 0.0328 | 7 of 1082 (0.3 to 1.3%) | 6 of 1082 (0.3 to 1.2%) | 0.0755 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0146 | 0.0350 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0819 |
| voxels with a C fit | s1 | 853 | 0.0130 | 0.0317 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0750 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0120 | 0.0279 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0757 |
| interior (C voxels) | s1 | 120 | 0.0123 | 0.0287 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0726 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0151 | 0.0364 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0829 |
| boundary (C voxels) | s1 | 733 | 0.0132 | 0.0329 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0754 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s1 | 209 | 0.0171 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0782 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0186 | 0.0378 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0868 |
| level 0.1 (C voxels) | s1 | 222 | 0.0154 | 0.0364 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0803 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0236 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0698 |
| level 0.2 (C voxels) | s1 | 201 | 0.0096 | 0.0306 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0647 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0183 | 0.0446 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0850 |
| level 0.5 (C voxels) | s1 | 144 | 0.0140 | 0.0353 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0758 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0095 | 0.0352 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0869 |
| level 1.0 (C voxels) | s1 | 77 | 0.0099 | 0.0317 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0772 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/100,V/30; accepted per stage [17, 20]; renders per voxel 149.0; stage wall time 7.7 s (7.2 ms per voxel); timing: contended: timing preflight refused (Python at 100% CPU); single run, one worker

## go2_As100_V61

args: `{"name": "go2_As100_V61", "pipeline": "Astrict/100,V/61", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0143 | 0.0367 | 8 of 1082 (0.4 to 1.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0837 |
| all | s1 | 1082 | 0.0100 | 0.0306 | 7 of 1082 (0.3 to 1.3%) | 6 of 1082 (0.3 to 1.2%) | 0.0740 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0146 | 0.0350 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0819 |
| voxels with a C fit | s1 | 853 | 0.0114 | 0.0302 | 3 of 853 (0.1 to 1.0%) | 3 of 853 (0.1 to 1.0%) | 0.0739 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0120 | 0.0279 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0757 |
| interior (C voxels) | s1 | 120 | 0.0116 | 0.0287 | 0 of 120 (0.0 to 3.1%) | 0 of 120 (0.0 to 3.1%) | 0.0722 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0151 | 0.0364 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0829 |
| boundary (C voxels) | s1 | 733 | 0.0114 | 0.0302 | 3 of 733 (0.1 to 1.2%) | 3 of 733 (0.1 to 1.2%) | 0.0742 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s1 | 209 | 0.0143 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0773 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0186 | 0.0378 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0868 |
| level 0.1 (C voxels) | s1 | 222 | 0.0146 | 0.0346 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0792 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0095 | 0.0236 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0698 |
| level 0.2 (C voxels) | s1 | 201 | 0.0088 | 0.0306 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0634 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0183 | 0.0446 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0850 |
| level 0.5 (C voxels) | s1 | 144 | 0.0116 | 0.0306 | 0 of 144 (0.0 to 2.6%) | 0 of 144 (0.0 to 2.6%) | 0.0746 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0095 | 0.0352 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0869 |
| level 1.0 (C voxels) | s1 | 77 | 0.0077 | 0.0317 | 2 of 77 (0.7 to 9.0%) | 2 of 77 (0.7 to 9.0%) | 0.0754 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages Astrict/100,V/61; accepted per stage [17, 20]; renders per voxel 175.1; stage wall time 10.2 s (9.4 ms per voxel); timing: contended: timing preflight refused (Claude at 26% CPU, Claude at 43% CPU, Python at 100% CPU); single run, one worker

## go2_V20

args: `{"name": "go2_V20", "pipeline": "V/20", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0195 | 0.0955 | 52 of 1082 (3.7 to 6.2%) | 6 of 1082 (0.3 to 1.2%) | 0.0810 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0235 | 0.1016 | 45 of 853 (4.0 to 7.0%) | 3 of 853 (0.1 to 1.0%) | 0.0820 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0279 | 0.1183 | 9 of 120 (4.0 to 13.6%) | 0 of 120 (0.0 to 3.1%) | 0.0804 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0219 | 0.0964 | 36 of 733 (3.6 to 6.7%) | 3 of 733 (0.1 to 1.2%) | 0.0823 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0194 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0793 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0181 | 0.0451 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0822 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0184 | 0.0553 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0702 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0400 | 0.1012 | 9 of 144 (3.3 to 11.5%) | 0 of 144 (0.0 to 2.6%) | 0.0862 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0742 | 0.2062 | 35 of 77 (34.8 to 56.5%) | 2 of 77 (0.7 to 9.0%) | 0.1122 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages V/20; accepted per stage [20]; renders per voxel 22.0; stage wall time 3.7 s (3.4 ms per voxel); timing: contended: timing preflight refused (Claude at 25% CPU, Claude at 50% CPU, Python at 97% CPU); single run, one worker

## go2_V61

args: `{"name": "go2_V61", "pipeline": "V/61", "sample": "grad", "start": "paint", "grains": 20, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 1082 | 0.0285 | 0.1302 | 94 of 1082 (7.2 to 10.5%) | 6 of 1082 (0.3 to 1.2%) | 0.0948 |
| all | s0 | 1082 | 0.0140 | 0.0753 | 41 of 1082 (2.8 to 5.1%) | 6 of 1082 (0.3 to 1.2%) | 0.0767 |
| all | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| voxels with a C fit | start | 853 | 0.0349 | 0.1354 | 84 of 853 (8.0 to 12.0%) | 3 of 853 (0.1 to 1.0%) | 0.0986 |
| voxels with a C fit | s0 | 853 | 0.0165 | 0.0829 | 34 of 853 (2.9 to 5.5%) | 3 of 853 (0.1 to 1.0%) | 0.0772 |
| voxels with a C fit | C | 853 | 0.0181 | 0.0551 | 18 of 853 (1.3 to 3.3%) | 3 of 853 (0.1 to 1.0%) | 0.0781 |
| interior (C voxels) | start | 120 | 0.0512 | 0.1498 | 14 of 120 (7.1 to 18.6%) | 0 of 120 (0.0 to 3.1%) | 0.0968 |
| interior (C voxels) | s0 | 120 | 0.0243 | 0.0903 | 6 of 120 (2.3 to 10.5%) | 0 of 120 (0.0 to 3.1%) | 0.0764 |
| interior (C voxels) | C | 120 | 0.0227 | 0.0555 | 3 of 120 (0.9 to 7.1%) | 0 of 120 (0.0 to 3.1%) | 0.0755 |
| boundary (C voxels) | start | 733 | 0.0324 | 0.1305 | 70 of 733 (7.6 to 11.9%) | 3 of 733 (0.1 to 1.2%) | 0.0989 |
| boundary (C voxels) | s0 | 733 | 0.0155 | 0.0788 | 28 of 733 (2.7 to 5.5%) | 3 of 733 (0.1 to 1.2%) | 0.0773 |
| boundary (C voxels) | C | 733 | 0.0176 | 0.0545 | 15 of 733 (1.2 to 3.3%) | 3 of 733 (0.1 to 1.2%) | 0.0786 |
| level 0.0 (C voxels) | start | 209 | 0.0165 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0843 |
| level 0.0 (C voxels) | s0 | 209 | 0.0146 | 0.0279 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0773 |
| level 0.0 (C voxels) | C | 209 | 0.0148 | 0.0337 | 0 of 209 (0.0 to 1.8%) | 0 of 209 (0.0 to 1.8%) | 0.0794 |
| level 0.1 (C voxels) | start | 222 | 0.0334 | 0.0649 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0958 |
| level 0.1 (C voxels) | s0 | 222 | 0.0155 | 0.0408 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0797 |
| level 0.1 (C voxels) | C | 222 | 0.0187 | 0.0463 | 0 of 222 (0.0 to 1.7%) | 0 of 222 (0.0 to 1.7%) | 0.0818 |
| level 0.2 (C voxels) | start | 201 | 0.0392 | 0.0764 | 4 of 201 (0.8 to 5.0%) | 1 of 201 (0.1 to 2.8%) | 0.0942 |
| level 0.2 (C voxels) | s0 | 201 | 0.0120 | 0.0415 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0642 |
| level 0.2 (C voxels) | C | 201 | 0.0162 | 0.0354 | 1 of 201 (0.1 to 2.8%) | 1 of 201 (0.1 to 2.8%) | 0.0664 |
| level 0.5 (C voxels) | start | 144 | 0.0684 | 0.1497 | 36 of 144 (18.6 to 32.7%) | 0 of 144 (0.0 to 2.6%) | 0.1112 |
| level 0.5 (C voxels) | s0 | 144 | 0.0233 | 0.0957 | 7 of 144 (2.4 to 9.7%) | 0 of 144 (0.0 to 2.6%) | 0.0799 |
| level 0.5 (C voxels) | C | 144 | 0.0216 | 0.0556 | 1 of 144 (0.1 to 3.8%) | 0 of 144 (0.0 to 2.6%) | 0.0792 |
| level 1.0 (C voxels) | start | 77 | 0.1269 | 0.2487 | 44 of 77 (46.0 to 67.6%) | 2 of 77 (0.7 to 9.0%) | 0.1330 |
| level 1.0 (C voxels) | s0 | 77 | 0.0460 | 0.2053 | 26 of 77 (24.2 to 44.9%) | 2 of 77 (0.7 to 9.0%) | 0.0983 |
| level 1.0 (C voxels) | C | 77 | 0.0247 | 0.2014 | 16 of 77 (13.2 to 31.1%) | 2 of 77 (0.7 to 9.0%) | 0.0927 |

grains: 20; stages V/61; accepted per stage [20]; renders per voxel 62.4; stage wall time 7.1 s (6.6 ms per voxel); timing: contended: timing preflight refused (Claude at 24% CPU, Claude at 41% CPU, Python at 100% CPU); single run, one worker

