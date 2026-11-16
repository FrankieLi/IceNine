## go2s_As100_V20

args: `{"name": "go2s_As100_V20", "pipeline": "Astrict/100,V/20", "sample": "grad", "start": "paint", "grains": 10, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 620 | 0.0294 | 0.1056 | 37 of 620 (4.4 to 8.1%) | 2 of 620 (0.1 to 1.2%) | 0.0938 |
| all | s0 | 620 | 0.0152 | 0.0402 | 4 of 620 (0.3 to 1.6%) | 2 of 620 (0.1 to 1.2%) | 0.0815 |
| all | s1 | 620 | 0.0143 | 0.0339 | 3 of 620 (0.2 to 1.4%) | 2 of 620 (0.1 to 1.2%) | 0.0761 |
| all | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| voxels with a C fit | start | 494 | 0.0358 | 0.1097 | 34 of 494 (5.0 to 9.5%) | 1 of 494 (0.0 to 1.1%) | 0.0981 |
| voxels with a C fit | s0 | 494 | 0.0165 | 0.0373 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0812 |
| voxels with a C fit | s1 | 494 | 0.0153 | 0.0323 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0759 |
| voxels with a C fit | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| interior (C voxels) | start | 108 | 0.0519 | 0.1375 | 11 of 108 (5.8 to 17.3%) | 0 of 108 (0.0 to 3.4%) | 0.0950 |
| interior (C voxels) | s0 | 108 | 0.0128 | 0.0279 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0739 |
| interior (C voxels) | s1 | 108 | 0.0145 | 0.0300 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0714 |
| interior (C voxels) | C | 108 | 0.0235 | 0.0479 | 1 of 108 (0.2 to 5.1%) | 0 of 108 (0.0 to 3.4%) | 0.0726 |
| boundary (C voxels) | start | 386 | 0.0317 | 0.1039 | 23 of 386 (4.0 to 8.8%) | 1 of 386 (0.0 to 1.5%) | 0.0990 |
| boundary (C voxels) | s0 | 386 | 0.0169 | 0.0407 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0832 |
| boundary (C voxels) | s1 | 386 | 0.0153 | 0.0334 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0772 |
| boundary (C voxels) | C | 386 | 0.0182 | 0.0494 | 2 of 386 (0.1 to 1.9%) | 1 of 386 (0.0 to 1.5%) | 0.0784 |
| level 0.0 (C voxels) | start | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s0 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s1 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0814 |
| level 0.0 (C voxels) | C | 116 | 0.0279 | 0.0340 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0809 |
| level 0.1 (C voxels) | start | 131 | 0.0391 | 0.0671 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.1068 |
| level 0.1 (C voxels) | s0 | 131 | 0.0184 | 0.0407 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0946 |
| level 0.1 (C voxels) | s1 | 131 | 0.0154 | 0.0381 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0880 |
| level 0.1 (C voxels) | C | 131 | 0.0195 | 0.0543 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0892 |
| level 0.2 (C voxels) | start | 129 | 0.0361 | 0.0711 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0868 |
| level 0.2 (C voxels) | s0 | 129 | 0.0097 | 0.0243 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0645 |
| level 0.2 (C voxels) | s1 | 129 | 0.0117 | 0.0307 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0603 |
| level 0.2 (C voxels) | C | 129 | 0.0147 | 0.0341 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0607 |
| level 0.5 (C voxels) | start | 79 | 0.0635 | 0.1588 | 22 of 79 (19.2 to 38.6%) | 0 of 79 (0.0 to 4.6%) | 0.0999 |
| level 0.5 (C voxels) | s0 | 79 | 0.0226 | 0.0480 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0753 |
| level 0.5 (C voxels) | s1 | 79 | 0.0170 | 0.0341 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0676 |
| level 0.5 (C voxels) | C | 79 | 0.0218 | 0.0609 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0690 |
| level 1.0 (C voxels) | start | 39 | 0.0847 | 0.1640 | 12 of 39 (18.6 to 46.4%) | 1 of 39 (0.5 to 13.2%) | 0.1404 |
| level 1.0 (C voxels) | s0 | 39 | 0.0066 | 0.0114 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0903 |
| level 1.0 (C voxels) | s1 | 39 | 0.0080 | 0.0214 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0877 |
| level 1.0 (C voxels) | C | 39 | 0.0160 | 0.1443 | 3 of 39 (2.7 to 20.3%) | 1 of 39 (0.5 to 13.2%) | 0.0961 |

grains: 10; stages Astrict/100,V/20; accepted per stage [8, 10]; objective evaluations per voxel 135.0; stage wall time 3.4 s (5.6 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 1.50 > 1.5); single run, one worker

## go2s_Aski150_As100_V20

args: `{"name": "go2s_Aski150_As100_V20", "pipeline": "Asinki/150,Astrict/100,V/20", "sample": "grad", "start": "paint", "grains": 10, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 620 | 0.0294 | 0.1056 | 37 of 620 (4.4 to 8.1%) | 2 of 620 (0.1 to 1.2%) | 0.0938 |
| all | s0 | 620 | 0.0217 | 0.0526 | 4 of 620 (0.3 to 1.6%) | 2 of 620 (0.1 to 1.2%) | 0.0868 |
| all | s1 | 620 | 0.0164 | 0.0425 | 4 of 620 (0.3 to 1.6%) | 2 of 620 (0.1 to 1.2%) | 0.0834 |
| all | s2 | 620 | 0.0150 | 0.0360 | 3 of 620 (0.2 to 1.4%) | 2 of 620 (0.1 to 1.2%) | 0.0768 |
| all | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| voxels with a C fit | start | 494 | 0.0358 | 0.1097 | 34 of 494 (5.0 to 9.5%) | 1 of 494 (0.0 to 1.1%) | 0.0981 |
| voxels with a C fit | s0 | 494 | 0.0258 | 0.0518 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0853 |
| voxels with a C fit | s1 | 494 | 0.0165 | 0.0409 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0821 |
| voxels with a C fit | s2 | 494 | 0.0155 | 0.0351 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0763 |
| voxels with a C fit | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| interior (C voxels) | start | 108 | 0.0519 | 0.1375 | 11 of 108 (5.8 to 17.3%) | 0 of 108 (0.0 to 3.4%) | 0.0950 |
| interior (C voxels) | s0 | 108 | 0.0161 | 0.0473 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0770 |
| interior (C voxels) | s1 | 108 | 0.0137 | 0.0279 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0743 |
| interior (C voxels) | s2 | 108 | 0.0153 | 0.0295 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0715 |
| interior (C voxels) | C | 108 | 0.0235 | 0.0479 | 1 of 108 (0.2 to 5.1%) | 0 of 108 (0.0 to 3.4%) | 0.0726 |
| boundary (C voxels) | start | 386 | 0.0317 | 0.1039 | 23 of 386 (4.0 to 8.8%) | 1 of 386 (0.0 to 1.5%) | 0.0990 |
| boundary (C voxels) | s0 | 386 | 0.0279 | 0.0532 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0876 |
| boundary (C voxels) | s1 | 386 | 0.0174 | 0.0427 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0843 |
| boundary (C voxels) | s2 | 386 | 0.0156 | 0.0371 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0777 |
| boundary (C voxels) | C | 386 | 0.0182 | 0.0494 | 2 of 386 (0.1 to 1.9%) | 1 of 386 (0.0 to 1.5%) | 0.0784 |
| level 0.0 (C voxels) | start | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s0 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s1 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s2 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0814 |
| level 0.0 (C voxels) | C | 116 | 0.0279 | 0.0340 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0809 |
| level 0.1 (C voxels) | start | 131 | 0.0391 | 0.0671 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.1068 |
| level 0.1 (C voxels) | s0 | 131 | 0.0251 | 0.0559 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0986 |
| level 0.1 (C voxels) | s1 | 131 | 0.0220 | 0.0480 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0947 |
| level 0.1 (C voxels) | s2 | 131 | 0.0156 | 0.0390 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0887 |
| level 0.1 (C voxels) | C | 131 | 0.0195 | 0.0543 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0892 |
| level 0.2 (C voxels) | start | 129 | 0.0361 | 0.0711 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0868 |
| level 0.2 (C voxels) | s0 | 129 | 0.0094 | 0.0165 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0648 |
| level 0.2 (C voxels) | s1 | 129 | 0.0099 | 0.0206 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0646 |
| level 0.2 (C voxels) | s2 | 129 | 0.0112 | 0.0307 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0603 |
| level 0.2 (C voxels) | C | 129 | 0.0147 | 0.0341 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0607 |
| level 0.5 (C voxels) | start | 79 | 0.0635 | 0.1588 | 22 of 79 (19.2 to 38.6%) | 0 of 79 (0.0 to 4.6%) | 0.0999 |
| level 0.5 (C voxels) | s0 | 79 | 0.0483 | 0.0555 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0899 |
| level 0.5 (C voxels) | s1 | 79 | 0.0155 | 0.0223 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0769 |
| level 0.5 (C voxels) | s2 | 79 | 0.0155 | 0.0349 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0678 |
| level 0.5 (C voxels) | C | 79 | 0.0218 | 0.0609 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0690 |
| level 1.0 (C voxels) | start | 39 | 0.0847 | 0.1640 | 12 of 39 (18.6 to 46.4%) | 1 of 39 (0.5 to 13.2%) | 0.1404 |
| level 1.0 (C voxels) | s0 | 39 | 0.0140 | 0.0393 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0977 |
| level 1.0 (C voxels) | s1 | 39 | 0.0140 | 0.0393 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0977 |
| level 1.0 (C voxels) | s2 | 39 | 0.0157 | 0.0474 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0904 |
| level 1.0 (C voxels) | C | 39 | 0.0160 | 0.1443 | 3 of 39 (2.7 to 20.3%) | 1 of 39 (0.5 to 13.2%) | 0.0961 |

grains: 10; stages Asinki/150,Astrict/100,V/20; accepted per stage [6, 8, 10]; objective evaluations per voxel 298.0; stage wall time 862.6 s (1391.3 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 1.64 > 1.5); single run, one worker

## go2s_Als150_As100_V20

args: `{"name": "go2s_Als150_As100_V20", "pipeline": "Als/150,Astrict/100,V/20", "sample": "grad", "start": "paint", "grains": 10, "true_grains": null, "group": "flood", "sigma0": 0.05, "b_lr": 0.01, "b_mu": 1.0}`

| group | stage | n | median | p95 | above 0.1 | wrong | strict |
|---|---|---|---|---|---|---|---|
| all | start | 620 | 0.0294 | 0.1056 | 37 of 620 (4.4 to 8.1%) | 2 of 620 (0.1 to 1.2%) | 0.0938 |
| all | s0 | 620 | 0.0132 | 0.0402 | 4 of 620 (0.3 to 1.6%) | 2 of 620 (0.1 to 1.2%) | 0.0815 |
| all | s1 | 620 | 0.0132 | 0.0402 | 4 of 620 (0.3 to 1.6%) | 2 of 620 (0.1 to 1.2%) | 0.0815 |
| all | s2 | 620 | 0.0127 | 0.0313 | 3 of 620 (0.2 to 1.4%) | 2 of 620 (0.1 to 1.2%) | 0.0761 |
| all | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| voxels with a C fit | start | 494 | 0.0358 | 0.1097 | 34 of 494 (5.0 to 9.5%) | 1 of 494 (0.0 to 1.1%) | 0.0981 |
| voxels with a C fit | s0 | 494 | 0.0149 | 0.0373 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0810 |
| voxels with a C fit | s1 | 494 | 0.0149 | 0.0373 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0810 |
| voxels with a C fit | s2 | 494 | 0.0136 | 0.0308 | 1 of 494 (0.0 to 1.1%) | 1 of 494 (0.0 to 1.1%) | 0.0758 |
| voxels with a C fit | C | 494 | 0.0196 | 0.0493 | 3 of 494 (0.2 to 1.8%) | 1 of 494 (0.0 to 1.1%) | 0.0771 |
| interior (C voxels) | start | 108 | 0.0519 | 0.1375 | 11 of 108 (5.8 to 17.3%) | 0 of 108 (0.0 to 3.4%) | 0.0950 |
| interior (C voxels) | s0 | 108 | 0.0107 | 0.0279 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0737 |
| interior (C voxels) | s1 | 108 | 0.0107 | 0.0279 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0737 |
| interior (C voxels) | s2 | 108 | 0.0123 | 0.0289 | 0 of 108 (0.0 to 3.4%) | 0 of 108 (0.0 to 3.4%) | 0.0714 |
| interior (C voxels) | C | 108 | 0.0235 | 0.0479 | 1 of 108 (0.2 to 5.1%) | 0 of 108 (0.0 to 3.4%) | 0.0726 |
| boundary (C voxels) | start | 386 | 0.0317 | 0.1039 | 23 of 386 (4.0 to 8.8%) | 1 of 386 (0.0 to 1.5%) | 0.0990 |
| boundary (C voxels) | s0 | 386 | 0.0164 | 0.0407 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0831 |
| boundary (C voxels) | s1 | 386 | 0.0164 | 0.0407 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0831 |
| boundary (C voxels) | s2 | 386 | 0.0137 | 0.0312 | 1 of 386 (0.0 to 1.5%) | 1 of 386 (0.0 to 1.5%) | 0.0771 |
| boundary (C voxels) | C | 386 | 0.0182 | 0.0494 | 2 of 386 (0.1 to 1.9%) | 1 of 386 (0.0 to 1.5%) | 0.0784 |
| level 0.0 (C voxels) | start | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s0 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s1 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0855 |
| level 0.0 (C voxels) | s2 | 116 | 0.0279 | 0.0279 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0814 |
| level 0.0 (C voxels) | C | 116 | 0.0279 | 0.0340 | 0 of 116 (0.0 to 3.2%) | 0 of 116 (0.0 to 3.2%) | 0.0809 |
| level 0.1 (C voxels) | start | 131 | 0.0391 | 0.0671 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.1068 |
| level 0.1 (C voxels) | s0 | 131 | 0.0138 | 0.0407 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0941 |
| level 0.1 (C voxels) | s1 | 131 | 0.0138 | 0.0407 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0941 |
| level 0.1 (C voxels) | s2 | 131 | 0.0116 | 0.0288 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0876 |
| level 0.1 (C voxels) | C | 131 | 0.0195 | 0.0543 | 0 of 131 (0.0 to 2.8%) | 0 of 131 (0.0 to 2.8%) | 0.0892 |
| level 0.2 (C voxels) | start | 129 | 0.0361 | 0.0711 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0868 |
| level 0.2 (C voxels) | s0 | 129 | 0.0097 | 0.0243 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0645 |
| level 0.2 (C voxels) | s1 | 129 | 0.0097 | 0.0243 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0645 |
| level 0.2 (C voxels) | s2 | 129 | 0.0117 | 0.0307 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0603 |
| level 0.2 (C voxels) | C | 129 | 0.0147 | 0.0341 | 0 of 129 (0.0 to 2.9%) | 0 of 129 (0.0 to 2.9%) | 0.0607 |
| level 0.5 (C voxels) | start | 79 | 0.0635 | 0.1588 | 22 of 79 (19.2 to 38.6%) | 0 of 79 (0.0 to 4.6%) | 0.0999 |
| level 0.5 (C voxels) | s0 | 79 | 0.0226 | 0.0480 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0753 |
| level 0.5 (C voxels) | s1 | 79 | 0.0226 | 0.0480 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0753 |
| level 0.5 (C voxels) | s2 | 79 | 0.0170 | 0.0341 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0676 |
| level 0.5 (C voxels) | C | 79 | 0.0218 | 0.0609 | 0 of 79 (0.0 to 4.6%) | 0 of 79 (0.0 to 4.6%) | 0.0690 |
| level 1.0 (C voxels) | start | 39 | 0.0847 | 0.1640 | 12 of 39 (18.6 to 46.4%) | 1 of 39 (0.5 to 13.2%) | 0.1404 |
| level 1.0 (C voxels) | s0 | 39 | 0.0066 | 0.0114 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0903 |
| level 1.0 (C voxels) | s1 | 39 | 0.0066 | 0.0114 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0903 |
| level 1.0 (C voxels) | s2 | 39 | 0.0080 | 0.0214 | 1 of 39 (0.5 to 13.2%) | 1 of 39 (0.5 to 13.2%) | 0.0877 |
| level 1.0 (C voxels) | C | 39 | 0.0160 | 0.1443 | 3 of 39 (2.7 to 20.3%) | 1 of 39 (0.5 to 13.2%) | 0.0961 |

grains: 10; stages Als/150,Astrict/100,V/20; accepted per stage [7, 8, 10]; objective evaluations per voxel 298.0; stage wall time 69.3 s (111.8 ms per voxel); timing: contended: timing preflight refused (machine not quiet for timing: 1-min load 1.59 > 1.5); single run, one worker

