# C++ IceNine elastic-strain / lattice-deformation survey (AR6 G1c, idea I31)

Read-only survey (2026-11-16) by a subagent through `git grep` / `git show` on the refs below,
plus the non-git checkout `~/Research/IceNine.others/IceNine.Strain`. Placed by the main session
from the text the agent returned (its tools were read-only).

## Summary
1. Strain exists only in the forward image simulator; reconstruction and fitting never touch it.
2. Per-voxel `oDeformation` M (3x3, symmetric on read) acts in the crystal frame: g = O M g0.
3. .mic: 19 columns, cols 13-18 = M00 M11 M22 M01 M12 M02 (identity if fewer columns).
4. A strain-optimizer config (StrainOpt.config) is parsed and serialized but never used.
5. All five branches are effectively identical for strain; IceNine.Strain is an older,
   less complete draft. No strain-resolution results or notes were found anywhere.

## Source 1: git refs in the IceNine repository
References are origin/master unless marked GF (= origin/GrainFitting).

### Parameterization
- `XDM++/libXDM/Voxel.h:56`: `SMatrix3x3 oDeformation; // This isn't exactly strain. More like a
  metric modifier`. Defaults to identity (Voxel.h:65).
- `Src/ForwardSimulation.cpp:324` (GF:328): "oDeformation is M, where M * A_0 = A, A is the
  'new' B matrix". So M multiplies the reciprocal basis: B' = M B0.
- There is no F, U, eps or lattice-parameter parameterization. The unit cell (a, b, c, alpha,
  beta, gamma) is fixed per phase: `XDM++/libXDM/CrystalStructure.cpp:217-237` (a along x, b in
  the xy plane) and `:261-282` (a* = 2 pi (b x c)/V). Units are 1/Angstrom with the 2 pi
  included; k = E[keV] * keV_over_hbar_c_in_ang (`Src/Simulation.cpp:73`).
- Physical reading (the agent's inference): if the direct lattice deforms as a' = F a (crystal
  frame), then B' = F^-T B0, so M = F^-T, about I - eps for small symmetric strain. The comment at
  ForwardSimulation.cpp:342 ("g'' = (I + S) g ... we've taken the inverse of (I+S)") says the
  stored matrix is already inverted, but does not state the sign convention for S. Nothing in the
  code fixes the sign, so AR6 defines it and tests it with a pure dilation.

### Forward model
- `Src/ForwardSimulation.cpp:294-370` `SimulatePeaksWithStrain_Debug` (GF:298):
  - :327 `g = oDeformation * g0` (crystal frame), then `g.Transform(oOrientMatrix)` (sample
    frame; O is the active ZXZ Euler matrix, `MicIO.cpp:265`).
  - :349 `GetScatteringOmegas(..., g, g.GetLength())`: Bragg angle and omega come from the
    strained |g| (`Simulation.cpp:73-90`: sin(theta) = |g|/2k).
  - The diffracted ray is a mirror reflection about the normalized, omega-rotated g
    (`Simulation.tmpl.cpp` near :129 GetReflectedRayDir), so 2-theta and eta follow the strained g.
  - Bug: :351 computes `fSinTheta` for the acceptance functor from the unstrained `fMag`. This
    affects only the eta / sin(2-theta) intensity filter, not peak positions.
- Reached only through `SimulateDetectorImage` when `bStrainEnabled` (`Src/Driver.cpp:133`,
  GF:156) and `ForwardSimulation.cpp:104,170`. The unstrained `SimulatePeaks` (:232-245) has only
  a placeholder comment.
- The reconstruction cost path ignores M: `Src/DiscreteSearch.h:277` "Need to change this to get
  to strain compatible - need a template parameter"; DiffractionCore.h uses unstrained vectors.

### Optimizer and config (stubs only)
- Main-config keys (`Src/ConfigFile.cpp:159,166,167`; handlers :663-686; GF:161,172,173,694-718):
  `EnableStrain 0|1`; `StrainOptConfigFilename <file>` (sets bStrainEnabled and parses the file);
  `LazyStrain` (pushes FIT_STRAIN_LBFS). ConfigFile.cpp:735 requires EnableStrain to agree with
  StrainOptConfigFilename being given.
- `Src/XDMServer.cpp:278-281` `RunStrainLazyBFS`: `RUNTIME_ASSERT(0, "Strain LazyBFS
  disabled")`; the body is commented out (it would have written `.strain` output).
- `Src/OptimizationConfig.h:62-131`, `.cpp:48-342`: `StrainConfigFile` with keys
  `StrainSearchType` [Global|Local|Mix], `StrainType` [Diag|Sym|Asym], `OptimizationMethod`
  [MonteCarlo|ConjugateGradient], `StrainStepSize` (3x3 matrix over 3 lines), `MaxStrainStep`
  (int), `ConvergenceMethod` [Quality fQ|HitRatio fH|Differential fTol]. The struct is serialized
  (ConfigFile.cpp near :900, :995) but nothing reads StrainSteps.
- Example `ConfigFiles/StrainOpt.config` ("Not all options are currently available"): Local,
  MaxStrainStep 100, StrainStepSize diag(0.001, 0.001, 0.001), Quality 1.0, MonteCarlo, Diag.
  These suggest a per-voxel MC with 1e-3 diagonal steps, orientation and geometry fixed (an
  inference; no implementation exists).
- `XDM++/libXDM/MicGrid.h:530`: Nye-tensor code assumes elastic strain is negligible.

### .mic I/O
- Read, `XDM++/libXDM/MicIO.cpp:274-288` (also MicIO.h:331, TriangularMic.h:397): only when a
  line has exactly 19 tokens; col 10 = cost, col 11 = overlap, col 12 = fitting time (skipped on
  read); cols 13-18 = M00 M11 M22 M01 M12 M02, mirrored to make M symmetric; otherwise identity.
- Write, MicIO.cpp:355-368: cost, overlap, time, then M00 M11 M22 M10 M12 M20 (same as the read
  order because M is symmetric).
- Problem: no precision is set on the stream, so the default is 6 significant digits. Diagonal
  entries near 1 are quantized to 1e-5, which is 10% of a 1e-4 strain.
- Data trail: GF `ConfigFiles/OneGrain.Config:64-68` points to StrainDataCornell simulated Ti
  (BeforeStrain10.mic, AfterStrain.mic). The AfterStrain.mic copy in XDM++ has only 10 columns,
  i.e. no deformation stored.

### Branch differences (git diff --stat on the strain files)
- master vs GF: ConfigFile.cpp only drops the DEBUG_*SearchType keyword guards;
  OptimizationConfig uses std::tie instead of boost::tie. Same strain logic.
- IceNine_v3.0.x: same as master.
- fix_pmc_opt: same logic, reformatted (strain simulation at ForwardSimulation.cpp:304-370).
- archive/local-snapshot-2025-11-28: adds a debug print and exit(0) in unstrained SimulatePeaks.
- None is more complete; GF and master are the cleanest to read.

## Source 2: `~/Research/IceNine.others/IceNine.Strain` (non-git, files dated 2017)
- Older draft. `XDM++.strain/libXDM/MicIO.cpp:261-270`: a 16-column format, cols 10-15 =
  S00 S11 S22 S01 S12 S02, stored as M = I + S ("strain tensor (S + I)"). Does not match the
  19-column layout above.
- `Src/ForwardSimulation.cpp:160-186`: strain is only a comment and is never applied. :186 has
  `fSinTheta = fMag / fWavenumber` (the factor 2 is missing).
- `Src/ConfigFile.cpp:124,594`: `EnableStrain` only; no strain-optimizer config.
- No notes on strain in TODO.txt or WishList.txt.
- `~/Research/IceNine.others/XDM++`: no strain math, only MicGrid.h:507 (the same Nye note) and
  AfterStrain.mic example inputs.

## Source 3: Python port, `icenine_py/icenine/mic_file.py`
- `Voxel.deformation` is an optional (3, 3) float32 (:63, :77, :89-94).
- `_parse_deformation_tensor` (:183-203) builds a symmetric matrix from [m00 m11 m22 m01 m12 m02].
- Parsed only when there are exactly 19 tokens (:253-256), matching C++.
- Write (:375-395) always emits 19 columns: time 0.0, deformation identity if None, order
  D00 D11 D22 D10 D12 D20, format 13.6E (7 significant digits, about 1e-6 near 1). Better than
  the C++ 6-digit default, but the float32 storage limits it to about 6e-8 near 1.
- Nothing else in `icenine/` uses `deformation`; the forward model ignores it.

## What to reuse for AR6
- Convention to adopt (matches the C++ data path and I14):
  g0 = B0 [h k l] (crystal frame, 2 pi / Angstrom); g = O M g0;
  M = F^-T, about (I + eps)^-1, about I - eps; sin(theta) = |g| / (2k); omega from the Bragg
  solution using |g|; the outgoing ray mirrors k_in about g/|g|. Here eps is the crystal-frame
  symmetric strain; the sample-frame strain is O eps O^T.
- Pin the sign with a pure dilation: eps = d I must give |g| = |g0| / (1 + d) and a smaller
  2-theta.
- .mic: keep the 19-column layout and store M (not eps) in cols 13-18 for C++ compatibility.
  Write at least %.9E for strain work; convert eps = M^-1 - I at the I/O boundary.
- Config keys to mirror for opt-in strain fitting: EnableStrain, StrainType (Diag | Sym = 6
  components | Asym), StrainSearchType (Local | Global), StrainStepSize (3x3, start 1e-3 and
  refine toward 1e-4), MaxStrainStep, ConvergenceMethod (Quality | HitRatio | Differential).
- Fix when porting: compute sin(2-theta) for the acceptance and intensity filter from the
  strained |g| (C++ uses the unstrained magnitude).
- There is no strain optimizer, cost function or resolution benchmark to port. Use HEXRD (I30)
  for the fitting method and these files only for the convention and I/O.
