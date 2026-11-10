<!-- table:s3_batched_12 -->
**ThreeVoxels Q-max 8, default rule**: evaluate 349 us (min 344), K = 112 reflections, 11 sync-forcing calls per evaluate

| device | variant | B | us | us_min | speedup | torch_us | host_us | stage_d_us | differs | counts_differ | peak_differs | cost_differs | max_dcost | med_dcost | pool | syncs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cpu | fixed | 1 | 286.3 | 271.4 | 1.22 | 269.6 | 17.9 | 14.0 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 4 | 145.8 | 139.1 | 2.40 | 135.0 | 6.7 | 11.2 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 16 | 109.4 | 106.2 | 3.19 | 94.9 | 3.8 | 10.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 64 | 98.9 | 97.1 | 3.53 | 84.5 | 3.6 | 10.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 256 | 99.9 | 97.8 | 3.50 | 82.9 | 3.4 | 13.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 1024 | 95.9 | 95.8 | 3.64 | 81.4 | 3.1 | 11.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | compact | 1 | 249.8 | 242.9 | 1.40 | 235.6 | 2.1 | 13.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 4 | 111.0 | 109.0 | 3.15 | 98.9 | 0.7 | 10.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 16 | 72.4 | 71.0 | 4.83 | 61.7 | 0.4 | 9.9 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 64 | 63.4 | 62.3 | 5.51 | 54.0 | 0.3 | 11.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 256 | 65.4 | 63.3 | 5.34 | 50.3 | 0.3 | 13.0 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 1024 | 62.4 | 62.0 | 5.60 | 50.4 | 0.3 | 11.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| mps | fixed | 1 | 2492.5 | 2195.5 | 0.14 | 2136.3 | 20.3 | 14.4 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
| mps | fixed | 4 | 566.0 | 524.0 | 0.62 | 541.2 | 7.6 | 10.7 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
| mps | fixed | 16 | 151.4 | 145.7 | 2.31 | 139.8 | 4.3 | 9.7 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
| mps | fixed | 64 | 62.3 | 51.8 | 5.60 | 38.9 | 3.6 | 9.9 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
| mps | fixed | 256 | 33.9 | 32.3 | 10.29 | 20.5 | 3.1 | 11.3 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
| mps | fixed | 1024 | 31.6 | 30.6 | 11.07 | 16.5 | 3.2 | 11.3 | 118 | 118 | 1 | 6 | 0.0035 | 0.00258 | 1024 | 1 |
<!-- /table:s3_batched_12 -->

<!-- table:s3_batched_13 -->
**ThreeVoxels Q-max 8, C++ rule**: evaluate 340 us (min 335), K = 112 reflections, 11 sync-forcing calls per evaluate

| device | variant | B | us | us_min | speedup | torch_us | host_us | stage_d_us | differs | counts_differ | peak_differs | cost_differs | max_dcost | med_dcost | pool | syncs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cpu | fixed | 1 | 279.8 | 275.3 | 1.22 | 256.8 | 15.2 | 12.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 4 | 141.8 | 139.4 | 2.40 | 126.9 | 6.4 | 10.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 16 | 107.3 | 105.0 | 3.17 | 94.1 | 3.7 | 10.2 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 64 | 98.7 | 97.2 | 3.45 | 84.2 | 3.4 | 10.2 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 256 | 99.7 | 98.2 | 3.41 | 82.1 | 3.2 | 14.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 1024 | 96.0 | 95.2 | 3.54 | 80.5 | 3.4 | 11.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | compact | 1 | 248.3 | 243.0 | 1.37 | 227.5 | 2.0 | 12.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 4 | 110.8 | 109.2 | 3.07 | 100.2 | 0.7 | 10.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 16 | 72.7 | 71.4 | 4.68 | 63.5 | 0.4 | 10.2 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 64 | 64.0 | 62.5 | 5.32 | 53.5 | 0.3 | 9.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 256 | 63.9 | 63.3 | 5.32 | 51.6 | 0.3 | 14.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 1024 | 64.6 | 64.1 | 5.26 | 51.4 | 0.3 | 12.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| mps | fixed | 1 | 2220.4 | 2098.7 | 0.15 | 2664.6 | 38.5 | 25.4 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
| mps | fixed | 4 | 621.1 | 565.6 | 0.55 | 546.7 | 7.3 | 11.0 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
| mps | fixed | 16 | 151.8 | 148.0 | 2.24 | 137.5 | 4.0 | 9.8 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
| mps | fixed | 64 | 61.3 | 49.6 | 5.56 | 46.0 | 3.7 | 11.5 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
| mps | fixed | 256 | 36.0 | 32.7 | 9.47 | 25.0 | 3.6 | 14.4 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
| mps | fixed | 1024 | 31.5 | 29.9 | 10.81 | 16.3 | 3.2 | 11.4 | 118 | 118 | 1 | 5 | 0.0035 | 0.00282 | 1024 | 1 |
<!-- /table:s3_batched_13 -->

<!-- table:s3_batched_14 -->
**ThreeVoxels all reflections (Simulation.config), default rule**: evaluate 886 us (min 831), K = 868 reflections, 11 sync-forcing calls per evaluate

| device | variant | B | us | us_min | speedup | torch_us | host_us | stage_d_us | differs | counts_differ | peak_differs | cost_differs | max_dcost | med_dcost | pool | syncs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cpu | fixed | 1 | 900.7 | 871.0 | 0.98 | 809.7 | 32.6 | 49.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 4 | 749.4 | 732.5 | 1.18 | 678.6 | 26.7 | 47.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 16 | 713.9 | 694.1 | 1.24 | 646.2 | 27.1 | 73.4 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 64 | 744.2 | 721.7 | 1.19 | 643.0 | 24.8 | 78.4 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 256 | 712.8 | 705.9 | 1.24 | 648.3 | 25.0 | 73.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 1024 | 702.2 | 695.1 | 1.26 | 614.3 | 22.5 | 65.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | compact | 1 | 628.1 | 606.8 | 1.41 | 570.1 | 3.6 | 50.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 4 | 477.4 | 467.6 | 1.86 | 427.1 | 2.1 | 48.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 16 | 443.2 | 435.8 | 2.00 | 397.2 | 2.2 | 68.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 64 | 483.1 | 463.9 | 1.83 | 384.8 | 2.1 | 82.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 256 | 455.6 | 443.6 | 1.94 | 375.1 | 2.0 | 70.0 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 1024 | 508.7 | 445.0 | 1.74 | 383.1 | 2.0 | 68.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| mps | fixed | 1 | 2378.8 | 2230.5 | 0.37 | 2211.8 | 41.7 | 52.7 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
| mps | fixed | 4 | 633.6 | 615.9 | 1.40 | 547.5 | 28.9 | 49.4 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
| mps | fixed | 16 | 274.3 | 263.9 | 3.23 | 201.9 | 28.0 | 60.1 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
| mps | fixed | 64 | 233.7 | 219.7 | 3.79 | 132.4 | 25.9 | 80.5 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
| mps | fixed | 256 | 227.6 | 222.3 | 3.89 | 127.5 | 24.4 | 68.1 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
| mps | fixed | 1024 | 231.1 | 228.9 | 3.83 | 140.6 | 23.0 | 62.7 | 391 | 388 | 2 | 31 | 0.0010 | 0.00042 | 1024 | 1 |
<!-- /table:s3_batched_14 -->

<!-- table:s3_batched_15 -->
**E0 clean voxel 0, default rule**: evaluate 326 us (min 326), K = 112 reflections, 11 sync-forcing calls per evaluate

| device | variant | B | us | us_min | speedup | torch_us | host_us | stage_d_us | differs | counts_differ | peak_differs | cost_differs | max_dcost | med_dcost | pool | syncs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cpu | fixed | 1 | 274.7 | 263.0 | 1.19 | 248.4 | 14.0 | 12.7 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 4 | 145.4 | 139.6 | 2.24 | 129.3 | 6.1 | 11.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 16 | 108.7 | 105.4 | 3.00 | 93.5 | 3.4 | 10.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 64 | 100.9 | 96.1 | 3.23 | 85.2 | 3.1 | 10.9 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 256 | 96.2 | 94.8 | 3.39 | 82.3 | 2.9 | 11.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 1024 | 94.8 | 94.3 | 3.44 | 81.3 | 2.9 | 10.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | compact | 1 | 246.2 | 240.4 | 1.33 | 232.7 | 2.2 | 13.0 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 4 | 106.4 | 104.0 | 3.07 | 96.3 | 0.7 | 11.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 16 | 66.9 | 66.0 | 4.88 | 56.6 | 0.3 | 10.4 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 64 | 58.1 | 57.6 | 5.62 | 47.5 | 0.2 | 10.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 256 | 55.2 | 54.8 | 5.90 | 45.4 | 0.2 | 11.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 1024 | 54.9 | 53.7 | 5.94 | 44.1 | 0.2 | 10.5 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| mps | fixed | 1 | 2320.6 | 2159.2 | 0.14 | 2273.0 | 26.0 | 17.4 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
| mps | fixed | 4 | 567.6 | 531.2 | 0.57 | 541.4 | 7.0 | 11.0 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
| mps | fixed | 16 | 150.9 | 147.4 | 2.16 | 137.8 | 3.7 | 10.0 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
| mps | fixed | 64 | 56.1 | 54.2 | 5.81 | 39.6 | 2.9 | 9.9 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
| mps | fixed | 256 | 32.7 | 31.5 | 9.99 | 19.2 | 2.7 | 10.4 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
| mps | fixed | 1024 | 28.9 | 28.5 | 11.29 | 15.5 | 2.6 | 10.4 | 112 | 112 | 0 | 19 | 0.0064 | 0.00091 | 1024 | 1 |
<!-- /table:s3_batched_15 -->

<!-- table:s3_batched_16 -->
**E0 clean voxel 0, C++ rule**: evaluate 325 us (min 324), K = 112 reflections, 11 sync-forcing calls per evaluate

| device | variant | B | us | us_min | speedup | torch_us | host_us | stage_d_us | differs | counts_differ | peak_differs | cost_differs | max_dcost | med_dcost | pool | syncs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| cpu | fixed | 1 | 278.3 | 264.3 | 1.17 | 249.0 | 14.6 | 13.2 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 4 | 143.3 | 140.1 | 2.27 | 131.5 | 5.9 | 11.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 16 | 108.8 | 106.6 | 2.99 | 94.3 | 3.4 | 10.4 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 64 | 100.0 | 97.5 | 3.25 | 85.6 | 3.2 | 11.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 256 | 96.6 | 95.1 | 3.36 | 82.7 | 3.0 | 11.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | fixed | 1024 | 95.6 | 94.9 | 3.40 | 81.2 | 3.0 | 16.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 1 |
| cpu | compact | 1 | 248.7 | 238.5 | 1.31 | 231.3 | 2.0 | 12.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 4 | 106.4 | 104.2 | 3.05 | 95.1 | 0.7 | 11.1 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 16 | 67.5 | 66.0 | 4.82 | 56.4 | 0.3 | 10.4 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 64 | 58.3 | 57.7 | 5.57 | 47.4 | 0.2 | 10.3 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 256 | 55.4 | 54.9 | 5.86 | 44.8 | 0.2 | 10.6 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| cpu | compact | 1024 | 55.4 | 54.7 | 5.87 | 44.5 | 0.2 | 10.8 | 0 | 0 | 0 | 0 | 0.0000 | 0.00000 | 1024 | 9 |
| mps | fixed | 1 | 2145.9 | 2097.8 | 0.15 | 2134.5 | 19.7 | 14.2 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
| mps | fixed | 4 | 549.4 | 521.2 | 0.59 | 534.9 | 7.0 | 11.0 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
| mps | fixed | 16 | 153.2 | 147.6 | 2.12 | 137.6 | 3.8 | 10.0 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
| mps | fixed | 64 | 49.3 | 48.2 | 6.60 | 39.5 | 2.8 | 9.8 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
| mps | fixed | 256 | 34.2 | 31.6 | 9.50 | 18.6 | 2.6 | 10.6 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
| mps | fixed | 1024 | 28.7 | 28.7 | 11.33 | 15.6 | 2.6 | 10.6 | 112 | 112 | 2 | 15 | 0.0111 | 0.00123 | 1024 | 1 |
<!-- /table:s3_batched_16 -->

<!-- table:s3_model_e0_11 -->
**e0: ReconstructQ8.config (E0 clean), cpp rule False**: modelled quick-MC time per run vs the measured sequential quick MC (trial-rotation build c_rot us per orientation); whole_search = Amdahl estimate from the quick-MC share of the wall time

| device | variant | scheme | k | c_rot | modelled_ms | sequential_ms | ratio | whole_search |
|---|---|---|---|---|---|---|---|---|
| cpu | compact | single | 1 | 14.0 | 188 | 986 | 5.24 | 1.67 |
| cpu | compact | single | 1 | 0.0 | 157 | 986 | 6.30 | 1.72 |
| cpu | compact | plan | 1 | 14.0 | 3855 | 986 | 0.26 | 0.41 |
| cpu | compact | plan | 1 | 0.0 | 3111 | 986 | 0.32 | 0.48 |
| cpu | compact | plan | 2 | 14.0 | 3842 | 986 | 0.26 | 0.41 |
| cpu | compact | plan | 2 | 0.0 | 3082 | 986 | 0.32 | 0.49 |
| cpu | compact | plan | 4 | 14.0 | 3972 | 986 | 0.25 | 0.40 |
| cpu | compact | plan | 4 | 0.0 | 3177 | 986 | 0.31 | 0.48 |
| cpu | compact | memo | 4 | 14.0 | 1090 | 986 | 0.90 | 0.95 |
| cpu | compact | memo | 4 | 0.0 | 883 | 986 | 1.12 | 1.05 |
| cpu | compact | memo | 16 | 14.0 | 1272 | 986 | 0.78 | 0.87 |
| cpu | compact | memo | 16 | 0.0 | 1021 | 986 | 0.97 | 0.98 |
| cpu | fixed | single | 1 | 14.0 | 300 | 986 | 3.28 | 1.53 |
| cpu | fixed | single | 1 | 0.0 | 269 | 986 | 3.67 | 1.57 |
| cpu | fixed | plan | 1 | 14.0 | 6100 | 986 | 0.16 | 0.28 |
| cpu | fixed | plan | 1 | 0.0 | 5355 | 986 | 0.18 | 0.31 |
| cpu | fixed | plan | 2 | 14.0 | 6102 | 986 | 0.16 | 0.28 |
| cpu | fixed | plan | 2 | 0.0 | 5342 | 986 | 0.18 | 0.31 |
| cpu | fixed | plan | 4 | 14.0 | 6299 | 986 | 0.16 | 0.27 |
| cpu | fixed | plan | 4 | 0.0 | 5504 | 986 | 0.18 | 0.31 |
| cpu | fixed | memo | 4 | 14.0 | 1725 | 986 | 0.57 | 0.73 |
| cpu | fixed | memo | 4 | 0.0 | 1518 | 986 | 0.65 | 0.79 |
| cpu | fixed | memo | 16 | 14.0 | 2018 | 986 | 0.49 | 0.66 |
| cpu | fixed | memo | 16 | 0.0 | 1766 | 986 | 0.56 | 0.72 |
| mps | fixed | single | 1 | 14.0 | 190 | 986 | 5.18 | 1.67 |
| mps | fixed | single | 1 | 0.0 | 159 | 986 | 6.21 | 1.71 |
| mps | fixed | plan | 1 | 14.0 | 3835 | 986 | 0.26 | 0.41 |
| mps | fixed | plan | 1 | 0.0 | 3095 | 986 | 0.32 | 0.48 |
| mps | fixed | plan | 2 | 14.0 | 3044 | 986 | 0.32 | 0.49 |
| mps | fixed | plan | 2 | 0.0 | 2272 | 986 | 0.43 | 0.61 |
| mps | fixed | plan | 4 | 14.0 | 2775 | 986 | 0.36 | 0.53 |
| mps | fixed | plan | 4 | 0.0 | 1976 | 986 | 0.50 | 0.67 |
| mps | fixed | memo | 4 | 14.0 | 1075 | 986 | 0.92 | 0.96 |
| mps | fixed | memo | 4 | 0.0 | 873 | 986 | 1.13 | 1.06 |
| mps | fixed | memo | 16 | 14.0 | 898 | 986 | 1.10 | 1.05 |
| mps | fixed | memo | 16 | 0.0 | 646 | 986 | 1.53 | 1.21 |
<!-- /table:s3_model_e0_11 -->

<!-- table:s3_model_three_1 -->
**three: Example2.Simulation.config, cpp rule False**: modelled quick-MC time per run vs the measured sequential quick MC (trial-rotation build c_rot us per orientation); whole_search = Amdahl estimate from the quick-MC share of the wall time

| device | variant | scheme | k | c_rot | modelled_ms | sequential_ms | ratio | whole_search |
|---|---|---|---|---|---|---|---|---|
| cpu | compact | single | 1 | 14.0 | 3105 | 5325 | 1.72 | 1.42 |
| cpu | compact | single | 1 | 0.0 | 3031 | 5325 | 1.76 | 1.44 |
| cpu | compact | plan | 1 | 14.0 | 127195 | 5325 | 0.04 | 0.06 |
| cpu | compact | plan | 1 | 0.0 | 123528 | 5325 | 0.04 | 0.06 |
| cpu | compact | plan | 2 | 14.0 | 132671 | 5325 | 0.04 | 0.06 |
| cpu | compact | plan | 2 | 0.0 | 128970 | 5325 | 0.04 | 0.06 |
| cpu | compact | plan | 4 | 14.0 | 139535 | 5325 | 0.04 | 0.05 |
| cpu | compact | plan | 4 | 0.0 | 135754 | 5325 | 0.04 | 0.05 |
| cpu | fixed | single | 1 | 14.0 | 4574 | 5325 | 1.16 | 1.11 |
| cpu | fixed | single | 1 | 0.0 | 4500 | 5325 | 1.18 | 1.12 |
| cpu | fixed | plan | 1 | 14.0 | 192233 | 5325 | 0.03 | 0.04 |
| cpu | fixed | plan | 1 | 0.0 | 188565 | 5325 | 0.03 | 0.04 |
| cpu | fixed | plan | 2 | 14.0 | 191996 | 5325 | 0.03 | 0.04 |
| cpu | fixed | plan | 2 | 0.0 | 188295 | 5325 | 0.03 | 0.04 |
| cpu | fixed | plan | 4 | 14.0 | 194847 | 5325 | 0.03 | 0.04 |
| cpu | fixed | plan | 4 | 0.0 | 191066 | 5325 | 0.03 | 0.04 |
| mps | fixed | single | 1 | 14.0 | 1592 | 5325 | 3.34 | 1.99 |
| mps | fixed | single | 1 | 0.0 | 1519 | 5325 | 3.51 | 2.03 |
| mps | fixed | plan | 1 | 14.0 | 64631 | 5325 | 0.08 | 0.11 |
| mps | fixed | plan | 1 | 0.0 | 60963 | 5325 | 0.09 | 0.12 |
| mps | fixed | plan | 2 | 14.0 | 64960 | 5325 | 0.08 | 0.11 |
| mps | fixed | plan | 2 | 0.0 | 61259 | 5325 | 0.09 | 0.12 |
| mps | fixed | plan | 4 | 14.0 | 66425 | 5325 | 0.08 | 0.11 |
| mps | fixed | plan | 4 | 0.0 | 62644 | 5325 | 0.09 | 0.12 |
<!-- /table:s3_model_three_1 -->

<!-- table:s3_model_three_3 -->
**three: ReconstructQ8.config, cpp rule False**: modelled quick-MC time per run vs the measured sequential quick MC (trial-rotation build c_rot us per orientation); whole_search = Amdahl estimate from the quick-MC share of the wall time

| device | variant | scheme | k | c_rot | modelled_ms | sequential_ms | ratio | whole_search |
|---|---|---|---|---|---|---|---|---|
| cpu | compact | single | 1 | 14.0 | 488 | 2261 | 4.63 | 1.88 |
| cpu | compact | single | 1 | 0.0 | 413 | 2261 | 5.47 | 1.96 |
| cpu | compact | plan | 1 | 14.0 | 32927 | 2261 | 0.07 | 0.11 |
| cpu | compact | plan | 1 | 0.0 | 27091 | 2261 | 0.08 | 0.13 |
| cpu | compact | plan | 2 | 14.0 | 32647 | 2261 | 0.07 | 0.11 |
| cpu | compact | plan | 2 | 0.0 | 26772 | 2261 | 0.08 | 0.13 |
| cpu | compact | plan | 4 | 14.0 | 32672 | 2261 | 0.07 | 0.11 |
| cpu | compact | plan | 4 | 0.0 | 26724 | 2261 | 0.08 | 0.13 |
| cpu | fixed | single | 1 | 14.0 | 708 | 2261 | 3.19 | 1.70 |
| cpu | fixed | single | 1 | 0.0 | 633 | 2261 | 3.57 | 1.76 |
| cpu | fixed | plan | 1 | 14.0 | 47341 | 2261 | 0.05 | 0.08 |
| cpu | fixed | plan | 1 | 0.0 | 41505 | 2261 | 0.05 | 0.09 |
| cpu | fixed | plan | 2 | 14.0 | 46945 | 2261 | 0.05 | 0.08 |
| cpu | fixed | plan | 2 | 0.0 | 41071 | 2261 | 0.06 | 0.09 |
| cpu | fixed | plan | 4 | 14.0 | 47008 | 2261 | 0.05 | 0.08 |
| cpu | fixed | plan | 4 | 0.0 | 41060 | 2261 | 0.06 | 0.09 |
| mps | fixed | single | 1 | 14.0 | 345 | 2261 | 6.56 | 2.03 |
| mps | fixed | single | 1 | 0.0 | 270 | 2261 | 8.38 | 2.11 |
| mps | fixed | plan | 1 | 14.0 | 22944 | 2261 | 0.10 | 0.15 |
| mps | fixed | plan | 1 | 0.0 | 17107 | 2261 | 0.13 | 0.20 |
| mps | fixed | plan | 2 | 14.0 | 20655 | 2261 | 0.11 | 0.17 |
| mps | fixed | plan | 2 | 0.0 | 14781 | 2261 | 0.15 | 0.23 |
| mps | fixed | plan | 4 | 14.0 | 19843 | 2261 | 0.11 | 0.18 |
| mps | fixed | plan | 4 | 0.0 | 13895 | 2261 | 0.16 | 0.25 |
<!-- /table:s3_model_three_3 -->

<!-- table:s3_model_three_5 -->
**three: ReconstructQ8.config, cpp rule True**: modelled quick-MC time per run vs the measured sequential quick MC (trial-rotation build c_rot us per orientation); whole_search = Amdahl estimate from the quick-MC share of the wall time

| device | variant | scheme | k | c_rot | modelled_ms | sequential_ms | ratio | whole_search |
|---|---|---|---|---|---|---|---|---|
| cpu | compact | single | 1 | 14.0 | 333 | 1542 | 4.63 | 1.76 |
| cpu | compact | single | 1 | 0.0 | 282 | 1542 | 5.46 | 1.82 |
| cpu | compact | plan | 1 | 14.0 | 16510 | 1542 | 0.09 | 0.16 |
| cpu | compact | plan | 1 | 0.0 | 13570 | 1542 | 0.11 | 0.19 |
| cpu | compact | plan | 2 | 14.0 | 16647 | 1542 | 0.09 | 0.16 |
| cpu | compact | plan | 2 | 0.0 | 13680 | 1542 | 0.11 | 0.19 |
| cpu | compact | plan | 4 | 14.0 | 16980 | 1542 | 0.09 | 0.15 |
| cpu | compact | plan | 4 | 0.0 | 13959 | 1542 | 0.11 | 0.18 |
| cpu | fixed | single | 1 | 14.0 | 485 | 1542 | 3.18 | 1.61 |
| cpu | fixed | single | 1 | 0.0 | 434 | 1542 | 3.55 | 1.66 |
| cpu | fixed | plan | 1 | 14.0 | 23975 | 1542 | 0.06 | 0.11 |
| cpu | fixed | plan | 1 | 0.0 | 21034 | 1542 | 0.07 | 0.13 |
| cpu | fixed | plan | 2 | 14.0 | 23887 | 1542 | 0.06 | 0.11 |
| cpu | fixed | plan | 2 | 0.0 | 20920 | 1542 | 0.07 | 0.13 |
| cpu | fixed | plan | 4 | 14.0 | 24033 | 1542 | 0.06 | 0.11 |
| cpu | fixed | plan | 4 | 0.0 | 21012 | 1542 | 0.07 | 0.13 |
| mps | fixed | single | 1 | 14.0 | 267 | 1542 | 5.77 | 1.84 |
| mps | fixed | single | 1 | 0.0 | 216 | 1542 | 7.12 | 1.90 |
| mps | fixed | plan | 1 | 14.0 | 12970 | 1542 | 0.12 | 0.20 |
| mps | fixed | plan | 1 | 0.0 | 10029 | 1542 | 0.15 | 0.25 |
| mps | fixed | plan | 2 | 14.0 | 11119 | 1542 | 0.14 | 0.23 |
| mps | fixed | plan | 2 | 0.0 | 8151 | 1542 | 0.19 | 0.30 |
| mps | fixed | plan | 4 | 14.0 | 10371 | 1542 | 0.15 | 0.24 |
| mps | fixed | plan | 4 | 0.0 | 7351 | 1542 | 0.21 | 0.32 |
<!-- /table:s3_model_three_5 -->

<!-- table:s3_model_three_sim_8 -->
**three_sim: ReconstructQ8.config, cpp rule False**: modelled quick-MC time per run vs the measured sequential quick MC (trial-rotation build c_rot us per orientation); whole_search = Amdahl estimate from the quick-MC share of the wall time

| device | variant | scheme | k | c_rot | modelled_ms | sequential_ms | ratio | whole_search |
|---|---|---|---|---|---|---|---|---|
| cpu | compact | single | 1 | 14.0 | 488 | 2269 | 4.65 | 1.87 |
| cpu | compact | single | 1 | 0.0 | 413 | 2269 | 5.49 | 1.94 |
| cpu | compact | plan | 1 | 14.0 | 33694 | 2269 | 0.07 | 0.11 |
| cpu | compact | plan | 1 | 0.0 | 27719 | 2269 | 0.08 | 0.13 |
| cpu | compact | plan | 2 | 14.0 | 33384 | 2269 | 0.07 | 0.11 |
| cpu | compact | plan | 2 | 0.0 | 27374 | 2269 | 0.08 | 0.13 |
| cpu | compact | plan | 4 | 14.0 | 33409 | 2269 | 0.07 | 0.11 |
| cpu | compact | plan | 4 | 0.0 | 27326 | 2269 | 0.08 | 0.13 |
| cpu | compact | memo | 4 | 14.0 | 6892 | 2269 | 0.33 | 0.45 |
| cpu | compact | memo | 4 | 0.0 | 5685 | 2269 | 0.40 | 0.53 |
| cpu | compact | memo | 16 | 14.0 | 7161 | 2269 | 0.32 | 0.44 |
| cpu | compact | memo | 16 | 0.0 | 5870 | 2269 | 0.39 | 0.52 |
| cpu | fixed | single | 1 | 14.0 | 708 | 2269 | 3.20 | 1.69 |
| cpu | fixed | single | 1 | 0.0 | 633 | 2269 | 3.58 | 1.74 |
| cpu | fixed | plan | 1 | 14.0 | 48444 | 2269 | 0.05 | 0.08 |
| cpu | fixed | plan | 1 | 0.0 | 42469 | 2269 | 0.05 | 0.09 |
| cpu | fixed | plan | 2 | 14.0 | 48006 | 2269 | 0.05 | 0.08 |
| cpu | fixed | plan | 2 | 0.0 | 41997 | 2269 | 0.05 | 0.09 |
| cpu | fixed | plan | 4 | 14.0 | 48069 | 2269 | 0.05 | 0.08 |
| cpu | fixed | plan | 4 | 0.0 | 41985 | 2269 | 0.05 | 0.09 |
| cpu | fixed | memo | 4 | 14.0 | 9912 | 2269 | 0.23 | 0.33 |
| cpu | fixed | memo | 4 | 0.0 | 8705 | 2269 | 0.26 | 0.37 |
| cpu | fixed | memo | 16 | 14.0 | 10308 | 2269 | 0.22 | 0.32 |
| cpu | fixed | memo | 16 | 0.0 | 9017 | 2269 | 0.25 | 0.36 |
| mps | fixed | single | 1 | 14.0 | 345 | 2269 | 6.58 | 2.01 |
| mps | fixed | single | 1 | 0.0 | 270 | 2269 | 8.41 | 2.09 |
| mps | fixed | plan | 1 | 14.0 | 23444 | 2269 | 0.10 | 0.15 |
| mps | fixed | plan | 1 | 0.0 | 17469 | 2269 | 0.13 | 0.20 |
| mps | fixed | plan | 2 | 14.0 | 21085 | 2269 | 0.11 | 0.17 |
| mps | fixed | plan | 2 | 0.0 | 15076 | 2269 | 0.15 | 0.23 |
| mps | fixed | plan | 4 | 14.0 | 20305 | 2269 | 0.11 | 0.18 |
| mps | fixed | plan | 4 | 0.0 | 14222 | 2269 | 0.16 | 0.24 |
| mps | fixed | memo | 4 | 14.0 | 5077 | 2269 | 0.45 | 0.58 |
| mps | fixed | memo | 4 | 0.0 | 3870 | 2269 | 0.59 | 0.71 |
| mps | fixed | memo | 16 | 14.0 | 4416 | 2269 | 0.51 | 0.64 |
| mps | fixed | memo | 16 | 0.0 | 3125 | 2269 | 0.73 | 0.82 |
<!-- /table:s3_model_three_sim_8 -->

<!-- table:s3_sim_e0 -->
**e0: ReconstructQ8.config (E0 clean), cpp rule False**, option A run for real (30 runs; sums over the levels of a run)

| k | rounds_per_run | counted | plan | memo | plan_over_counted | memo_over_counted | all_match |
|---|---|---|---|---|---|---|---|
| 1 | 90.2 | 2758 | 58177 | 15117 | 21.1 | 5.5 | True |
| 2 | 46.2 | 2758 | 59379 | 15387 | 21.5 | 5.6 | True |
| 4 | 24.5 | 2758 | 61793 | 15921 | 22.4 | 5.8 | True |
| 8 | 13.3 | 2758 | 66782 | 16979 | 24.2 | 6.2 | True |
| 16 | 8.0 | 2758 | 77004 | 19133 | 27.9 | 6.9 | True |
<!-- /table:s3_sim_e0 -->

<!-- table:s3_sim_three_sim -->
**three_sim: ReconstructQ8.config, cpp rule False**, option A run for real (1 runs; sums over the levels of a run)

| k | rounds_per_run | counted | plan | memo | plan_over_counted | memo_over_counted | all_match |
|---|---|---|---|---|---|---|---|
| 1 | 283.0 | 6420 | 427860 | 85896 | 66.6 | 13.4 | True |
| 4 | 73.0 | 6420 | 435590 | 87290 | 67.8 | 13.6 | True |
| 16 | 20.0 | 6420 | 468750 | 93288 | 73.0 | 14.5 | True |
<!-- /table:s3_sim_three_sim -->

<!-- table:s3_trace_e0_9 -->
**e0: ReconstructQ8.config (E0 clean), cpp rule False** (30 runs)

| level | runs | cand_median | cand_min | cand_max | failed | cand | rate | lo | hi | rounds_k1 | rounds_k2 | rounds_k4 | comp_over_counted_k1 | comp_over_counted_k2 | comp_over_counted_k4 | draws_ok | speedup_cpu_compact | speedup_cpu_fixed | speedup_mps_fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 30 | 168 | 138 | 229 | 1653 | 5118 | 0.323 | 0.310 | 0.336 | 55.6 | 28.1 | 14.4 | 25.2 | 25.6 | 26.5 | True | 5.9 | 3.4 | 8.9 |
| 1 | 30 | 44 | 32 | 63 | 766 | 1359 | 0.564 | 0.537 | 0.590 | 25.9 | 13.2 | 6.9 | 11.3 | 11.7 | 12.6 | True | 5.5 | 3.2 | 4.6 |
| 2 | 30 | 11 | 8 | 17 | 210 | 343 | 0.612 | 0.560 | 0.662 | 7.3 | 3.9 | 2.2 | 3.3 | 3.8 | 4.9 | True | 4.5 | 2.9 | 1.5 |
| 3 | 30 | 2 | 2 | 4 | 28 | 76 | 0.368 | 0.269 | 0.481 | 1.4 | 1.1 | 1.0 | 1.2 | 1.9 | 3.5 | True | 2.1 | 1.7 | 0.3 |
<!-- /table:s3_trace_e0_9 -->

<!-- table:s3_trace_three_0 -->
**three: Example2.Simulation.config, cpp rule False** (3 runs)

| level | runs | cand_median | cand_min | cand_max | failed | cand | rate | lo | hi | rounds_k1 | rounds_k2 | rounds_k4 | comp_over_counted_k1 | comp_over_counted_k2 | comp_over_counted_k4 | draws_ok | speedup_cpu_compact | speedup_cpu_fixed | speedup_mps_fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3 | 392 | 392 | 392 | 354 | 1176 | 0.301 | 0.275 | 0.328 | 119.0 | 59.7 | 30.3 | 49.9 | 50.3 | 51.1 | True | 1.8 | 1.3 | 3.9 |
| 1 | 3 | 101 | 99 | 108 | 164 | 308 | 0.532 | 0.477 | 0.587 | 55.3 | 28.0 | 14.3 | 25.5 | 25.9 | 26.7 | True | 1.9 | 1.2 | 3.8 |
| 2 | 3 | 25 | 24 | 29 | 52 | 78 | 0.667 | 0.556 | 0.761 | 17.7 | 9.3 | 5.0 | 7.4 | 7.9 | 8.8 | True | 1.9 | 1.2 | 3.5 |
| 3 | 3 | 6 | 6 | 7 | 11 | 19 | 0.579 | 0.363 | 0.769 | 3.7 | 2.0 | 1.3 | 2.4 | 2.7 | 3.9 | True | 1.9 | 1.2 | 1.9 |
| 4 | 3 | 1 | 1 | 1 | 0 | 3 | 0.000 | 0.000 | 0.562 | 1.0 | 1.0 | 1.0 | 1.0 | 1.8 | 3.5 | True | 1.4 | 1.0 | 0.4 |
| 5 | 3 | 1 | 1 | 1 | 1 | 3 | 0.333 | 0.061 | 0.792 | 1.0 | 1.0 | 1.0 | 1.0 | 1.8 | 3.5 | True | 1.4 | 1.0 | 0.4 |
<!-- /table:s3_trace_three_0 -->

<!-- table:s3_trace_three_2 -->
**three: ReconstructQ8.config, cpp rule False** (3 runs)

| level | runs | cand_median | cand_min | cand_max | failed | cand | rate | lo | hi | rounds_k1 | rounds_k2 | rounds_k4 | comp_over_counted_k1 | comp_over_counted_k2 | comp_over_counted_k4 | draws_ok | speedup_cpu_compact | speedup_cpu_fixed | speedup_mps_fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3 | 397 | 397 | 397 | 556 | 1191 | 0.467 | 0.439 | 0.495 | 186.0 | 93.3 | 47.0 | 78.5 | 78.9 | 79.7 | True | 5.5 | 3.6 | 10.6 |
| 1 | 3 | 105 | 104 | 106 | 213 | 315 | 0.676 | 0.623 | 0.725 | 71.3 | 36.0 | 18.3 | 29.2 | 29.6 | 30.5 | True | 5.4 | 3.5 | 7.3 |
| 2 | 3 | 27 | 26 | 27 | 51 | 80 | 0.637 | 0.528 | 0.734 | 17.0 | 8.7 | 4.7 | 7.1 | 7.5 | 8.4 | True | 5.2 | 3.4 | 3.4 |
| 3 | 3 | 6 | 6 | 6 | 7 | 18 | 0.389 | 0.203 | 0.614 | 3.0 | 1.7 | 1.0 | 2.0 | 2.5 | 3.5 | True | 3.7 | 2.7 | 0.9 |
<!-- /table:s3_trace_three_2 -->

<!-- table:s3_trace_three_4 -->
**three: ReconstructQ8.config, cpp rule True** (3 runs)

| level | runs | cand_median | cand_min | cand_max | failed | cand | rate | lo | hi | rounds_k1 | rounds_k2 | rounds_k4 | comp_over_counted_k1 | comp_over_counted_k2 | comp_over_counted_k4 | draws_ok | speedup_cpu_compact | speedup_cpu_fixed | speedup_mps_fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 3 | 273 | 273 | 273 | 390 | 819 | 0.476 | 0.442 | 0.510 | 130.3 | 65.7 | 33.0 | 57.6 | 58.0 | 58.8 | True | 5.3 | 3.4 | 9.6 |
| 1 | 3 | 69 | 68 | 70 | 148 | 207 | 0.715 | 0.650 | 0.772 | 50.0 | 25.3 | 13.0 | 22.5 | 22.9 | 23.8 | True | 5.3 | 3.4 | 5.8 |
| 2 | 3 | 17 | 17 | 17 | 33 | 51 | 0.647 | 0.510 | 0.764 | 11.0 | 6.0 | 3.3 | 5.1 | 5.5 | 6.5 | True | 4.7 | 3.2 | 2.4 |
| 3 | 3 | 4 | 4 | 4 | 4 | 12 | 0.333 | 0.138 | 0.609 | 1.7 | 1.0 | 1.0 | 1.3 | 1.8 | 3.5 | True | 3.1 | 2.4 | 0.5 |
<!-- /table:s3_trace_three_4 -->

<!-- table:s3_trace_three_sim_6 -->
**three_sim: ReconstructQ8.config, cpp rule False** (1 runs)

| level | runs | cand_median | cand_min | cand_max | failed | cand | rate | lo | hi | rounds_k1 | rounds_k2 | rounds_k4 | comp_over_counted_k1 | comp_over_counted_k2 | comp_over_counted_k4 | draws_ok | speedup_cpu_compact | speedup_cpu_fixed | speedup_mps_fixed |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | 1 | 397 | 397 | 397 | 189 | 397 | 0.476 | 0.427 | 0.525 | 190.0 | 95.0 | 48.0 | 81.4 | 81.7 | 82.5 | True | 5.5 | 3.6 | 10.6 |
| 1 | 1 | 105 | 105 | 105 | 74 | 105 | 0.705 | 0.612 | 0.784 | 74.0 | 37.0 | 19.0 | 29.9 | 30.3 | 31.2 | True | 5.4 | 3.5 | 7.3 |
| 2 | 1 | 27 | 27 | 27 | 17 | 27 | 0.630 | 0.442 | 0.785 | 17.0 | 9.0 | 5.0 | 7.5 | 7.9 | 8.8 | True | 5.2 | 3.4 | 3.4 |
| 3 | 1 | 6 | 6 | 6 | 1 | 6 | 0.167 | 0.030 | 0.564 | 2.0 | 1.0 | 1.0 | 1.4 | 1.8 | 3.5 | True | 3.7 | 2.7 | 0.9 |
<!-- /table:s3_trace_three_sim_6 -->
