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
| total | 4241 |
| both spots have no true partner | 0 |
| one spot has no true partner | 92 |
| both have a true partner, but another | 4149 |
<!-- /table:evalD_step2_false_pairs -->

<!-- table:evalD_step2_origin_error -->
| ray z-slope (abs u_z) | pairs | lab-frame error median (um) | lab-frame error p90 (um) | sample-frame error median (um) | sample-frame error p90 (um) | predicted sigma median (um) |
|---|---|---|---|---|---|---|
| 0 to 0.02 | 2131 | 30.8 | 92.8 | 30.7 | 92.8 | 364.8 |
| 0.02 to 0.05 | 5569 | 4.2 | 19.7 | 4.7 | 19.8 | 180.7 |
| 0.05 to 0.1 | 5838 | 1.6 | 4.3 | 2.4 | 5.2 | 74.7 |
| 0.1 to 1 | 9714 | 1.0 | 2.6 | 1.9 | 4.0 | 37.0 |
<!-- /table:evalD_step2_origin_error -->

<!-- table:evalD_step2_presets -->
| preset | predicted pairs | precision | link recall | peak recall | recall, one grain | recall, several grains |
|---|---|---|---|---|---|---|
| hull_tight | 56175 | 0.463 [0.459, 0.467] (26023/56175) | 0.933 [0.930, 0.936] (26023/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.885 [0.875, 0.894] (3882/4386) |
| hull_plan | 69314 | 0.376 [0.372, 0.379] (26034/69314) | 0.934 [0.931, 0.936] (26034/27888) | 0.947 | 0.942 [0.939, 0.945] (22141/23502) | 0.888 [0.878, 0.897] (3893/4386) |
| ring_tight | 29029 | 0.854 [0.850, 0.858] (24788/29029) | 0.889 [0.885, 0.892] (24788/27888) | 0.887 | 0.990 [0.988, 0.991] (23256/23502) | 0.349 [0.335, 0.364] (1532/4386) |
| ring_plan | 56027 | 0.487 [0.483, 0.491] (27283/56027) | 0.978 [0.977, 0.980] (27283/27888) | 0.980 | 0.998 [0.998, 0.999] (23459/23502) | 0.872 [0.862, 0.881] (3824/4386) |
| shape_tight | 25798 | 0.909 [0.905, 0.913] (23452/25798) | 0.841 [0.837, 0.845] (23452/27888) | 0.844 | 0.966 [0.964, 0.968] (22705/23502) | 0.170 [0.159, 0.182] (747/4386) |
| shape_plan | 44653 | 0.545 [0.540, 0.550] (24338/44653) | 0.873 [0.869, 0.877] (24338/27888) | 0.884 | 0.974 [0.972, 0.976] (22897/23502) | 0.329 [0.315, 0.343] (1441/4386) |
| ring_tight_c0.5 | 25910 | 0.927 [0.923, 0.930] (24007/25910) | 0.861 [0.857, 0.865] (24007/27888) | 0.856 | 0.985 [0.983, 0.986] (23139/23502) | 0.198 [0.186, 0.210] (868/4386) |
| ring_tight_c2 | 35907 | 0.723 [0.718, 0.727] (25954/35907) | 0.931 [0.928, 0.934] (25954/27888) | 0.932 | 0.993 [0.992, 0.994] (23344/23502) | 0.595 [0.580, 0.610] (2610/4386) |
| shape_tight_c0.5 | 24090 | 0.957 [0.955, 0.960] (23058/24090) | 0.827 [0.822, 0.831] (23058/27888) | 0.826 | 0.961 [0.959, 0.964] (22588/23502) | 0.107 [0.098, 0.117] (470/4386) |
| shape_tight_c2 | 29689 | 0.807 [0.802, 0.811] (23948/29689) | 0.859 [0.855, 0.863] (23948/27888) | 0.868 | 0.970 [0.967, 0.972] (22788/23502) | 0.264 [0.252, 0.278] (1160/4386) |
| ring_tight_binary | 28200 | 0.849 [0.845, 0.853] (23942/28200) | 0.859 [0.854, 0.863] (23942/27888) | 0.852 | 0.957 [0.954, 0.959] (22488/23502) | 0.332 [0.318, 0.346] (1454/4386) |
| shape_tight_binary | 25027 | 0.906 [0.903, 0.910] (22684/25027) | 0.813 [0.809, 0.818] (22684/27888) | 0.812 | 0.935 [0.931, 0.938] (21968/23502) | 0.163 [0.153, 0.174] (716/4386) |
<!-- /table:evalD_step2_presets -->

<!-- table:evalD_step2_rejections -->
| first failing gate | true links |
|---|---|
| gate_ray | 0 |
| gate_hull | 36 |
| gate_sigma | 0 |
| gate_ring | 3036 |
| gate_size | 0 |
| passes all gates | 24816 |
<!-- /table:evalD_step2_rejections -->

<!-- table:evalD_step3_attach_spread300 -->
| attachment | true peaks found | attached pairs per voxel | false attached per voxel | share of attached that are false | recall, interior | recall, boundary |
|---|---|---|---|---|---|---|
| A r=2 um | 0.004 [0.003, 0.005] (71/16717) | 0.3 | 0.1 | 0.253 | 0.014 | 0.002 |
| A r=5 um | 0.029 [0.026, 0.031] (477/16717) | 2.1 | 0.5 | 0.234 | 0.081 | 0.018 |
| A r=10 um | 0.118 [0.113, 0.123] (1973/16717) | 8.6 | 2.0 | 0.237 | 0.228 | 0.095 |
| A r=15 um | 0.252 [0.245, 0.258] (4208/16717) | 19.2 | 5.2 | 0.269 | 0.406 | 0.220 |
| A r=20 um | 0.423 [0.416, 0.431] (7072/16717) | 35.1 | 11.5 | 0.328 | 0.607 | 0.385 |
| A r=30 um | 0.702 [0.695, 0.709] (11732/16717) | 81.0 | 41.8 | 0.517 | 0.783 | 0.685 |
| A r=50 um | 0.830 [0.824, 0.836] (13875/16717) | 225.9 | 179.7 | 0.795 | 0.839 | 0.828 |
| A r=100 um | 0.872 [0.867, 0.877] (14572/16717) | 862.4 | 813.8 | 0.944 | 0.866 | 0.873 |
| B footprint | 0.870 [0.865, 0.875] (14545/16717) | 434.4 | 386.0 | 0.888 | 0.864 | 0.871 |
<!-- /table:evalD_step3_attach_spread300 -->

<!-- table:evalD_step3_attach_voxels30 -->
| attachment | true peaks found | attached pairs per voxel | false attached per voxel | share of attached that are false | recall, interior | recall, boundary |
|---|---|---|---|---|---|---|
| A r=2 um | 0.004 [0.002, 0.008] (7/1714) | 0.3 | 0.1 | 0.222 | 0.008 | 0.002 |
| A r=5 um | 0.020 [0.014, 0.028] (34/1714) | 1.7 | 0.5 | 0.320 | 0.049 | 0.003 |
| A r=10 um | 0.127 [0.112, 0.144] (218/1714) | 9.6 | 2.3 | 0.240 | 0.320 | 0.016 |
| A r=15 um | 0.231 [0.212, 0.252] (396/1714) | 19.7 | 6.5 | 0.330 | 0.494 | 0.079 |
| A r=20 um | 0.372 [0.349, 0.395] (637/1714) | 35.8 | 14.6 | 0.407 | 0.546 | 0.271 |
| A r=30 um | 0.663 [0.640, 0.685] (1136/1714) | 82.8 | 44.9 | 0.542 | 0.670 | 0.658 |
| A r=50 um | 0.821 [0.802, 0.838] (1407/1714) | 218.1 | 171.2 | 0.785 | 0.817 | 0.823 |
| A r=100 um | 0.863 [0.846, 0.878] (1479/1714) | 817.5 | 768.2 | 0.940 | 0.858 | 0.866 |
| B footprint | 0.851 [0.833, 0.867] (1458/1714) | 386.8 | 338.2 | 0.874 | 0.852 | 0.850 |
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
| 0 | 0.02 | 0.736 | 0.865 | 0.184 |
| 0 | 0.05 | 0.792 | 0.925 | 0.377 |
| 0 | 0.1 | 0.817 | 0.956 | 0.594 |
| 0 | 0.2 | 0.825 | 0.980 | 0.719 |
| 1 | 0.02 | 0.736 | 0.865 | 0.217 |
| 1 | 0.05 | 0.792 | 0.925 | 0.555 |
| 1 | 0.1 | 0.817 | 0.956 | 1.003 |
| 1 | 0.2 | 0.825 | 0.980 | 1.526 |
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
