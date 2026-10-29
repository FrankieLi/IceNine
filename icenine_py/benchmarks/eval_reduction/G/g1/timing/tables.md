<!-- table:evalG1_cache_build -->
| timing | SmoothCost constructor (s) | lit-pixel trees (s) | observed rays (s) | build total (s) |
|---|---|---|---|---|
| single-worker | 0.002 | 4.6 | 12.1 | 16.7 |
<!-- /table:evalG1_cache_build -->

<!-- table:evalG1_cache_timing -->
| arm (SmoothCostRank) | timing | key off, 15901 / 17034 (s) | 15901: cold cache / warm cache (s) | 17034: fresh build / warm cache (s) | saved per voxel by the cache (s), A / B | extra time over key off, as before (fresh build), B | extra time over key off, warm cache, A / B | warm equals fresh (B), cold equals warm (A) |
|---|---|---|---|---|---|---|---|---|
| E3+E4,tau=0.25 | single-worker | 112.1 / 119.2 | 157.3 / 139.7 | 161.6 / 144.4 | 17.6 / 17.2 | +35.5% | +24.6% / +21.1% | yes / yes |
| E1+E2+E3+E4 | single-worker | 112.1 / 119.2 | 257.1 / 239.3 | 263.9 / 246.8 | 17.8 / 17.1 | +121.3% | +113.5% / +106.9% | yes / yes |
<!-- /table:evalG1_cache_timing -->
