# k4DTF test samples

Inputs and tools for testing the decay-tree fitter on FCC-ee fast-simulation
samples (Pythia8 + EvtGen + Delphes IDEA, EDM4hep output):

| Step | What | Where |
|---|---|---|
| 1 | Generate signal samples | `generation/generate.sh` |
| 2 | Build candidates, attach MC truth, run k4DTF → per-event ntuple | `analysis/stage1_dtf_inputs.py` (FCCAnalyses) |
| 3 | Re-run k4DTF on that ntuple → one row per candidate with pulls | `k4dtf-validate` (`tools/`) |

Decay modes (each with its charge conjugate):

| Mode | Descriptor | Nodes (pre-order) |
|---|---|---|
| `Bd2DstPi` | `B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+` | `B Dst D0 D0_K D0_pi Dst_pi B_pi` |
| `Bu2D0K_D02KPi` | `B- -> (D0 -> K- pi+) K-` | `B D0 D0_K D0_pi B_K` |
| `Bu2D0K_D02KSPiPi` | `B- -> (D0 -> (KS0 -> pi+ pi-) pi+) K-` | `B D0 KS KS_pip KS_pim D0_pip D0_pim B_K` |
| `Bu2D0Pi_D02KPiPi0` | `B- -> (D0 -> K- pi+ (pi0 -> gamma gamma)) pi-` | `B D0 D0_K D0_pi D0_pi0 pi0_g1 pi0_g2 B_pi` |
| `Bd2D0Pi0_D02KK` | `B0 -> (D~0 -> K- K+) (pi0 -> gamma gamma) ` | `B D0 D0_Km D0_Kp D0_pi0 pi0_g1 pi0_g2` |

## Quick start

```bash
source /cvmfs/sw.hsf.org/key4hep/setup.sh

# k4DTF (from the repository root)
cmake -B build -DCMAKE_INSTALL_PREFIX=$PWD/install && cmake --build build -j8 && cmake --install build
k4_local_repo

# 1. 10k events, seed 1 -> samples/p8_ee_Zbb_ecm91_EvtGen_Bd2DstPi_seed1.root
examples/dtf_test_samples/generation/generate.sh Bd2DstPi 10000 1 samples

# 2. candidates + truth + in-situ fit (FCCAnalyses)
fccanalysis run examples/dtf_test_samples/analysis/stage1_dtf_inputs.py -j 4 \
    -i samples/p8_ee_Zbb_ecm91_EvtGen_Bd2DstPi_seed1.root \
    -o dtf_inputs_Bd2DstPi.root -- --mode Bd2DstPi

# 3. stand-alone fit and validation tree (re-run after every change to the fitter)
k4dtf-validate --input dtf_inputs_Bd2DstPi.root --output dtf_validation_Bd2DstPi.root \
    --descriptor "B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+" \
    --labels B,Dst,D0,D0_K,D0_pi,Dst_pi,B_pi --mass-constraint 421 --pv
```

Every generated event contains one signal B; the other b hadron decays
generically, so the samples also contain combinatorial background candidates.
Fit columns are `-999` for fits that did not converge.

**Signal only.** For high-multiplicity modes, building and fitting every
candidate is slow. Add `--signal-only` after `--` in step 2 to fit only the
truth-matched candidates. Candidates are then built only from reconstructed
particles whose MC particle comes from a decay of the head particle (any
charge, so mixed B0 count), which removes most of the combinatorics up front.
Only candidates where every node is matched are kept. The ntuple then has no
background candidates (`nCand` counts the signal ones), so use it for
resolutions and pulls, not efficiencies. The primary-vertex fit still uses all
tracks.

## Stage-1 ntuple (`events` tree)

One entry per event with at least one candidate (`--keep-empty` keeps all).
Per-candidate columns are vectors of length `nCand`. `<L>` is a node label
from the table above; `-999` marks missing values.

**Event**

| Column | Type | Content |
|---|---|---|
| `nCand` | int | number of candidates |
| `Bz` | float | magnetic field [T] (`magFieldBz` from the file) |
| `PV_reco` | `edm4hep::VertexData` | FCCAnalyses primary-vertex fit (`VertexFitterSimple`, beam-spot constraint 4.5 µm × 20 nm × 300 µm) |
| `PV_nTracks` | int | tracks in that fit |
| `PV_true_{x,y,z}` | double | true primary vertex [mm] |

**Per candidate**

| Column | Type | Content |
|---|---|---|
| `cand_isCC` | int | built for the charge-conjugate descriptor |
| `cand_isSignal` | int | every node matched to the generated decay chain |
| `cand_nMatched` | int | number of matched nodes |
| `cand_nLeavesInPV` | int | candidate tracks used in `PV_reco` (PV biased if > 0) |
| `<L>_mc` | `edm4hep::MCParticleData` | matched MC particle of every node (PDG 0 if unmatched) |
| `<leaf>_recoIdx` | int | index in `ReconstructedParticles` |
| `<leaf>_reco` | `edm4hep::ReconstructedParticleData` | with the descriptor's mass hypothesis and PDG; 4-momentum covariance propagated from the track |
| `<leaf>_track` | `edm4hep::TrackState` | reconstructed helix with covariance (fitter input) |
| `<leaf>_trueTrack` | `edm4hep::TrackState` | exact helix of the matched MC particle, same conventions, zero covariance |
| `<photon>_cluster` | `edm4hep::ClusterData` | photon leaves: cluster position, energy and errors as given to the fitter (see Photons below) |
| `<node>_M` | float | invariant mass before the fit [GeV] |
| `<node>_trueProdVtx_{x,y,z}`, `<node>_trueDecayVtx_{x,y,z}` | double | true vertices [mm] |
| `<node>_trueCtau` | double | true proper decay length [mm] |
| `deltaM` | float | `Dst_M - D0_M` (Bd2DstPi only) |

**Per candidate and fit configuration** (`DTF`: no constraints; `DTFPV`: mass and primary-vertex constraints)

| Column | Content |
|---|---|
| `<fit>_status`, `_chi2`, `_ndf`, `_nIter`, `_prob` | fit summary; status is the `k4DTF::ErrCode` bit pattern, 0 = success |
| `<fit>_<L>_M`, `_MErr`, `_p{x,y,z}` | fitted mass, its uncertainty, momentum of every node |
| `<fit>_<node>_vtx_{x,y,z}`, `_vtxErr_{x,y,z}` | fitted decay vertex of composite nodes [mm] |
| `<fit>_<node>_ctau`, `_ctauErr`, `_L`, `_LErr` | proper decay length and decay length [mm]; composites with a lifetime only (not the D*) |

## Validation tree (`k4dtf-validate`)

One row per candidate (`--signal-only` keeps matched ones): fit summary,
`<L>_{true,fit}{M,P}`, `<L>_pullP`, decay-vertex and `ctau` values and pulls for
composite nodes, primary-vertex pulls (`pv_pull_{x,y,z}`), and for every leaf
the pulls of the input track parameters against the MC helix
(`<leaf>_pull{D0,Phi,Omega,Z0,TanLambda}`),or for photons the pull of the
cluster energy (`<photon>_pullE`). These pulls do not depend on the
fit and should be standard normal already; they are the first check that the
inputs and conventions below are right.

## Conventions

- **Units**: mm, GeV, T. Masses of leaves and constraints are those of EvtGen's
  `evt.pdl`, which the samples are generated with (`k4DTF::ParticleTable`).
- **Track parameters** (`edm4hep::TrackState` as written by k4SimDelphes, see
  `k4DTF::DTFTrack`): `(D0, phi, omega, Z0, tanLambda)` at reference point
  (0,0,0), `phi` is the momentum azimuth at the point of closest approach,
  `omega = -2C` with Delphes' half-curvature `C`, so `sign(omega) = sign(q·Bz)`
  and `pT = 0.2998e-3 · Bz / |omega|`; the PCA is `(-D0 sin phi, D0 cos phi, Z0)`.
  Covariance: first 15 elements of `covMatrix`, lower triangle in that
  parameter order; the time row is zero.
- **Reconstructed particles** from Delphes have the pion mass, PDG ±211 and an
  empty 4-momentum covariance. `<leaf>_reco` and the fitter inputs use the
  descriptor's mass hypothesis and the covariance propagated from the track.
- **Photons** are the neutral reconstructed particles whose cluster is in
  `EFlowPhoton`. The fitter needs the cluster position and the errors of
  position and energy. k4SimDelphes writes them from v00-09 (Key4hep nightlies
  from September 2026); older files have the energy only. What is missing is
  taken from `CALORIMETER` in the stage-1 script (IDEA card defaults): the
  position where the photon direction meets the calorimeter front face, the
  position error from the tower size, and the ECAL energy resolution. Delphes
  places a photon uniformly within its 0.02 × 0.02 tower, so the angular
  error with respect to the true impact point is 0.02/√6 = 8.2 mrad.
  k4SimDelphes v00-09 writes 0.02/√12, half the variance. In the toy
  (B- -> D0(K- pi+ pi0) pi-) that gives ⟨chi2⟩ = 7.3 for ndf 6 and
  pi0/D0 momentum pulls of 1.25; `--photon-errors model` uses the
  `CALORIMETER` errors also for files that have their own.  
- **PID** is MC truth (the PDG of the linked MC particle), i.e. perfect.
- **Truth matching**: a composite node is matched if all its daughters are, they
  share one MC mother with the right PDG code, and that mother has no other
  daughters except photons. The true decay vertex of a node is the production
  vertex of its daughters (k4SimDelphes does not fill `endpoint`); for the head,
  the production vertex is taken before any B0 oscillation in the MC record.
- **Primary vertex**: the Pythia card smears the luminous region with the same
  4.5 µm × 20 nm × 300 µm used as beam-spot constraint in the PV fit; keep the
  two consistent (`BEAMSPOT_UM` in the stage-1 script).

## Caveats

- `DelphesPythia8EvtGen_EDM4HEP_k4Interface` seeds EvtGen with a fixed value
  (1000), so jobs with different Pythia seeds can share EvtGen random sequences.
  Prefer one long job per mode, or check for repeated signal decay times when
  merging jobs.
- Delphes track parameters come from a parametrised covariance smearing of the
  true helix, not from a track fit; pulls are Gaussian by construction, so
  these samples test the fitter, not the tracking.
- The primary-vertex fit uses all tracks; signal tracks can enter it
  (`cand_nLeavesInPV`).

## Adding a decay mode

1. Write an EvtGen user decay file with a signal alias (see `generation/*.dec`)
   and add a case to `generate.sh` (signal PDG and alias).
2. Add an entry to `MODES` in `analysis/stage1_dtf_inputs.py`: descriptor, one
   label per node in pre-order, mass windows, fit configurations.
3. Pass the same descriptor and labels to `k4dtf-validate`.
