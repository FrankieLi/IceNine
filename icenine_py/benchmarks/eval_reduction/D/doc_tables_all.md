<!-- table:evalD_step1_components -->
| detector | class | components | with a true partner | median pixels | p90 pixels | median sources | max sources |
|---|---|---|---|---|---|---|---|
| 0 | single-source | 10 | 0 | 10 | 11 | 1 | 1 |
| 0 | one grain, several voxels | 28411 | 24671 | 97 | 223 | 44 | 185 |
| 0 | several grains | 1770 | 1630 | 230 | 462 | 106 | 296 |
| 1 | single-source | 20 | 20 | 3 | 10 | 1 | 1 |
| 1 | one grain, several voxels | 26355 | 26355 | 103 | 235 | 44 | 185 |
| 1 | several grains | 823 | 823 | 249 | 472 | 103 | 296 |
<!-- /table:evalD_step1_components -->

<!-- table:evalD_step1_links -->
| class of pair | true links | peaks (voxel, reflection) |
|---|---|---|
| single-source | 0 | 0 |
| one grain, several voxels | 23502 | 1134225 |
| several grains | 4386 | 240828 |
<!-- /table:evalD_step1_links -->

<!-- table:evalD_step2_false_pairs -->
| cause | pairs |
|---|---|
| total | 3101 |
| both spots have no true partner | 0 |
| one spot has no true partner | 63 |
| both have a true partner, but another | 3038 |
<!-- /table:evalD_step2_false_pairs -->

<!-- table:evalD_step2_origin_error -->
| ray z-slope (abs u_z) | pairs | lab-frame error median (um) | lab-frame error p90 (um) | sample-frame error median (um) | sample-frame error p90 (um) | predicted sigma median (um) |
|---|---|---|---|---|---|---|
| 0 to 0.02 | 2127 | 30.8 | 92.8 | 30.8 | 92.8 | 302.0 |
| 0.02 to 0.05 | 5555 | 4.2 | 19.6 | 4.6 | 19.8 | 149.8 |
| 0.05 to 0.1 | 5821 | 1.6 | 4.3 | 2.4 | 5.1 | 61.7 |
| 0.1 to 1 | 9702 | 1.0 | 2.6 | 1.9 | 4.0 | 30.6 |
<!-- /table:evalD_step2_origin_error -->

<!-- table:evalD_step2_presets -->
| preset | predicted pairs | precision | link recall | peak recall | recall, one grain | recall, several grains |
|---|---|---|---|---|---|---|
| hull_tight | 54646 | 0.476 [0.472, 0.480] (26019/54646) | 0.933 [0.930, 0.936] (26019/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.884 [0.874, 0.893] (3878/4386) |
| hull_plan | 68984 | 0.377 [0.374, 0.381] (26034/68984) | 0.934 [0.931, 0.936] (26034/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.888 [0.878, 0.897] (3893/4386) |
| ring_tight | 27513 | 0.887 [0.883, 0.891] (24412/27513) | 0.875 [0.871, 0.879] (24412/27888) | 0.872 | 0.987 [0.986, 0.989] (23208/23502) | 0.275 [0.262, 0.288] (1204/4386) |
| ring_plan | 54706 | 0.497 [0.493, 0.501] (27196/54706) | 0.975 [0.973, 0.977] (27196/27888) | 0.977 | 0.998 [0.997, 0.998] (23453/23502) | 0.853 [0.843, 0.864] (3743/4386) |
| shape_tight | 24968 | 0.932 [0.929, 0.935] (23272/24968) | 0.834 [0.830, 0.839] (23272/27888) | 0.836 | 0.964 [0.962, 0.966] (22656/23502) | 0.140 [0.130, 0.151] (616/4386) |
| shape_plan | 43452 | 0.560 [0.555, 0.564] (24324/43452) | 0.872 [0.868, 0.876] (24324/27888) | 0.883 | 0.974 [0.972, 0.976] (22892/23502) | 0.326 [0.313, 0.341] (1432/4386) |
| ring_tight_c0.5 | 25121 | 0.945 [0.942, 0.947] (23728/25121) | 0.851 [0.847, 0.855] (23728/27888) | 0.845 | 0.980 [0.978, 0.982] (23033/23502) | 0.158 [0.148, 0.170] (695/4386) |
| ring_tight_c2 | 32754 | 0.778 [0.773, 0.782] (25481/32754) | 0.914 [0.910, 0.917] (25481/27888) | 0.914 | 0.992 [0.991, 0.993] (23313/23502) | 0.494 [0.480, 0.509] (2168/4386) |
| shape_tight_c0.5 | 23610 | 0.968 [0.966, 0.970] (22857/23610) | 0.820 [0.815, 0.824] (22857/27888) | 0.818 | 0.957 [0.954, 0.959] (22481/23502) | 0.086 [0.078, 0.094] (376/4386) |
| shape_tight_c2 | 27890 | 0.852 [0.848, 0.856] (23756/27890) | 0.852 [0.848, 0.856] (23756/27888) | 0.858 | 0.968 [0.966, 0.971] (22759/23502) | 0.227 [0.215, 0.240] (997/4386) |
| ring_tight_binary | 26111 | 0.880 [0.876, 0.884] (22969/26111) | 0.824 [0.819, 0.828] (22969/27888) | 0.815 | 0.929 [0.926, 0.932] (21835/23502) | 0.259 [0.246, 0.272] (1134/4386) |
| shape_tight_binary | 23634 | 0.927 [0.924, 0.931] (21919/23634) | 0.786 [0.781, 0.791] (21919/27888) | 0.782 | 0.908 [0.904, 0.912] (21344/23502) | 0.131 [0.121, 0.141] (575/4386) |
<!-- /table:evalD_step2_presets -->

<!-- table:evalD_step2_rejections -->
| first failing gate | true links |
|---|---|
| gate_ray | 0 |
| gate_hull | 48 |
| gate_sigma | 0 |
| gate_ring | 3405 |
| gate_size | 0 |
| passes all gates | 24435 |
<!-- /table:evalD_step2_rejections -->

<!-- table:evalD_step3_attach_spread300 -->
| attachment | true peaks found | attached pairs per voxel | false attached per voxel | share of attached that are false | recall, interior | recall, boundary |
|---|---|---|---|---|---|---|
| A r=2 um | 0.004 [0.003, 0.005] (71/16717) | 0.3 | 0.1 | 0.211 | 0.014 | 0.002 |
| A r=5 um | 0.028 [0.026, 0.031] (474/16717) | 2.0 | 0.4 | 0.199 | 0.081 | 0.018 |
| A r=10 um | 0.118 [0.113, 0.123] (1966/16717) | 8.3 | 1.7 | 0.207 | 0.227 | 0.095 |
| A r=15 um | 0.251 [0.245, 0.258] (4197/16717) | 18.4 | 4.4 | 0.239 | 0.405 | 0.219 |
| A r=20 um | 0.421 [0.414, 0.429] (7040/16717) | 33.7 | 10.2 | 0.303 | 0.604 | 0.384 |
| A r=30 um | 0.697 [0.690, 0.704] (11653/16717) | 77.8 | 38.9 | 0.501 | 0.776 | 0.681 |
| A r=50 um | 0.821 [0.815, 0.827] (13731/16717) | 217.3 | 171.6 | 0.789 | 0.826 | 0.820 |
| A r=100 um | 0.858 [0.853, 0.863] (14348/16717) | 827.8 | 779.9 | 0.942 | 0.850 | 0.860 |
| B footprint | 0.859 [0.854, 0.864] (14364/16717) | 422.2 | 374.3 | 0.887 | 0.849 | 0.861 |
<!-- /table:evalD_step3_attach_spread300 -->

<!-- table:evalD_step3_attach_voxels30 -->
| attachment | true peaks found | attached pairs per voxel | false attached per voxel | share of attached that are false | recall, interior | recall, boundary |
|---|---|---|---|---|---|---|
| A r=2 um | 0.004 [0.002, 0.008] (7/1714) | 0.3 | 0.1 | 0.222 | 0.008 | 0.002 |
| A r=5 um | 0.020 [0.014, 0.028] (34/1714) | 1.7 | 0.5 | 0.320 | 0.049 | 0.003 |
| A r=10 um | 0.127 [0.112, 0.143] (217/1714) | 9.3 | 2.0 | 0.219 | 0.320 | 0.015 |
| A r=15 um | 0.230 [0.211, 0.250] (394/1714) | 18.9 | 5.8 | 0.306 | 0.492 | 0.078 |
| A r=20 um | 0.370 [0.348, 0.394] (635/1714) | 34.5 | 13.4 | 0.387 | 0.545 | 0.270 |
| A r=30 um | 0.661 [0.638, 0.683] (1133/1714) | 79.8 | 42.1 | 0.527 | 0.667 | 0.657 |
| A r=50 um | 0.813 [0.794, 0.830] (1393/1714) | 208.8 | 162.3 | 0.778 | 0.804 | 0.818 |
| A r=100 um | 0.844 [0.826, 0.861] (1447/1714) | 783.6 | 735.4 | 0.938 | 0.839 | 0.847 |
| B footprint | 0.840 [0.822, 0.857] (1440/1714) | 374.7 | 326.7 | 0.872 | 0.836 | 0.843 |
<!-- /table:evalD_step3_attach_voxels30 -->

<!-- table:evalD_step3_separation -->
| wrong orientations | score | median best wrong / truth | max over voxels | voxels >= 0.5 truth | voxels >= truth |
|---|---|---|---|---|---|
| random axis, 1 degree and more | pixel overlap (C++ rule) | 0.110 | 0.173 | 0 | 0 |
| random axis, 1 degree and more | attributed, A r=15.7 um | 0.178 | 1.035 | 5 | 1 |
| random axis, 1 degree and more | attributed, A r=30 um | 0.138 | 0.983 | 2 | 0 |
| random axis, 1 degree and more | attributed, B footprint | 0.139 | 0.260 | 0 | 0 |
| random axis, 1 degree and more | attributed, oracle (true pairs) | 0.139 | 0.260 | 0 | 0 |
| random axis, 0.25 to 0.5 degree | pixel overlap (C++ rule) | 0.640 | 0.818 | 28 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, A r=15.7 um | 0.962 | 1.035 | 23 | 11 |
| random axis, 0.25 to 0.5 degree | attributed, A r=30 um | 0.849 | 1.526 | 30 | 3 |
| random axis, 0.25 to 0.5 degree | attributed, B footprint | 0.830 | 0.921 | 30 | 0 |
| random axis, 0.25 to 0.5 degree | attributed, oracle (true pairs) | 0.830 | 0.921 | 30 | 0 |
| rotations about <100>, <110>, <111> | pixel overlap (C++ rule) | 0.310 | 0.350 | 0 | 0 |
| rotations about <100>, <110>, <111> | attributed, A r=15.7 um | 0.316 | 1.000 | 4 | 1 |
| rotations about <100>, <110>, <111> | attributed, A r=30 um | 0.309 | 0.667 | 2 | 0 |
| rotations about <100>, <110>, <111> | attributed, B footprint | 0.310 | 0.435 | 0 | 0 |
| rotations about <100>, <110>, <111> | attributed, oracle (true pairs) | 0.307 | 0.435 | 0 | 0 |
<!-- /table:evalD_step3_separation -->

<!-- table:evalD_step3_sharpness -->
| misorientation (deg) | pixel overlap (C++ rule) | attributed, A r=15.7 um | attributed, A r=30 um | attributed, B footprint | attributed, oracle (true pairs) |
|---|---|---|---|---|---|
| 0.05 | 0.927 | 0.976 | 0.967 | 0.958 | 0.958 |
| 0.1 | 0.819 | 0.944 | 0.922 | 0.921 | 0.920 |
| 0.25 | 0.512 | 0.764 | 0.732 | 0.731 | 0.733 |
| 0.5 | 0.235 | 0.338 | 0.338 | 0.342 | 0.342 |
| 1 | 0.061 | 0.072 | 0.085 | 0.087 | 0.088 |
| 2 | 0.017 | 0.012 | 0.017 | 0.020 | 0.018 |
| 5 | 0.004 | 0.000 | 0.003 | 0.003 | 0.003 |
| 10 | 0.001 | 0.000 | 0.000 | 0.000 | 0.000 |
| 20 | 0.003 | 0.000 | 0.000 | 0.000 | 0.000 |
<!-- /table:evalD_step3_sharpness -->

<!-- table:evalD_step3_tolerance -->
| frame window | direction tolerance (deg) | B completeness at truth (median) | oracle purity at truth (median) | B half-score radius (deg) |
|---|---|---|---|---|
| 0 | 0.02 | 0.736 | 0.885 | 0.184 |
| 0 | 0.05 | 0.792 | 0.954 | 0.377 |
| 0 | 0.1 | 0.812 | 0.979 | 0.593 |
| 0 | 0.2 | 0.814 | 0.991 | 0.717 |
| 1 | 0.02 | 0.736 | 0.885 | 0.217 |
| 1 | 0.05 | 0.792 | 0.954 | 0.555 |
| 1 | 0.1 | 0.812 | 0.979 | 1.003 |
| 1 | 0.2 | 0.821 | 0.991 | 1.528 |
<!-- /table:evalD_step3_tolerance -->

<!-- table:evalD_step3_truth_scores -->
| score | median at the truth | minimum at the truth | half-score radius (deg) |
|---|---|---|---|
| pixel overlap (C++ rule) | 0.936 | 0.894 | 0.258 |
| attributed, A r=15.7 um | 0.078 | 0.000 | 0.384 |
| attributed, A r=30 um | 0.719 | 0.018 | 0.376 |
| attributed, B footprint | 0.792 | 0.696 | 0.377 |
| attributed, oracle (true pairs) | 0.789 | 0.696 | 0.378 |
<!-- /table:evalD_step3_truth_scores -->
