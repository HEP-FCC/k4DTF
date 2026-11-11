# k4DTF — Key4hep Decay Tree Fitter

A kinematic decay-tree fitter for the [Key4hep](https://key4hep.github.io/key4hep-doc/) software stack.
The algorithm is exposed as a standalone C++ library, a [Gaudi](https://gaudi.web.cern.ch/) algorithm
via [k4FWCore](https://github.com/key4hep/k4FWCore), and a set of
[RDataFrame](https://root.cern/doc/master/classROOT_1_1RDataFrame.html)-compatible free functions
for use inside [FCCAnalyses](https://github.com/HEP-FCC/FCCAnalyses).

---

## Table of contents

- [Overview](#overview)
- [Repository structure](#repository-structure)
- [Dependencies](#dependencies)
- [Development](#development)
  - [1. Source the Key4hep nightlies](#1-source-the-key4hep-nightlies)
  - [2. Build and install](#2-build-and-install)
  - [3. Activate the local build](#3-activate-the-local-build-with-k4_local_repo)
  - [4. Run the tests](#4-run-the-tests)
- [Usage](#usage)
  - [Core C++ library](#core-c-library)
  - [Gaudi / k4FWCore](#gaudi--k4fwcore)
  - [FCCAnalyses / RDataFrame](#fccanalyses--rdataframe)
- [Algorithm](#algorithm)
- [References](#references)
- [Contributing](#contributing)
- [Acknowledgements](#acknowledgements)

---

## Overview

The decay-tree fitter (DTF) performs a least-squares kinematic fit of a full
decay chain.  Starting from the measured four-momenta and their covariance
matrices of the final-state particles, it enforces a set of constraints
(four-momentum conservation at each decay node, optional mass constraints,
optional vertex constraints) using Lagrange multipliers, following the
approach Hulsbergen 2024 (arXiv https://arxiv.org/abs/physics/0503191, NIMA Nucl.Instrum.Meth. A552 (2005) 566-575.)

The fit output includes:

- Fitted four-momenta, vertex positions and updated covariances for all particles in the decay
  tree.
- A chi-square value and number of degrees of freedom.

---

## Repository structure

```
k4DTF/
├── CMakeLists.txt                        # Top-level build
├── cmake/k4DTFConfig.cmake.in            # Makes k4DTF find_package()-able
│
├── k4DTF/                                # Core library
│   ├── include/k4DTF/
│   │   ├── DecayTreeFitter.h             # Fitter
│   │   ├── FitResult.h                   # chi2, status, per-node results
│   │   ├── ErrCode.h, ChiSquare.h
│   │   ├── core/DecayNodeInput.h         # Fitter input: tree of nodes
│   │   ├── core/ParticleMeasurement.h    # Track (helix), photon, composite, ... measurements
│   │   ├── core/Helix.h                  # Helix conventions and helpers
│   │   ├── DecayDescriptor.h             # Descriptor parser (nested decays, pre-order nodes)
│   │   ├── ParticleTable.h               # Particle names, PDG codes, masses, resonances
│   │   ├── DTFParticle.h                 # Particle + track with ROOT/EDM4hep conversions
│   │   └── TreeInput.h                   # Descriptor + particles -> DecayNodeInput
│   └── src/
│       ├── DecayTreeFitter.cpp           # Iteration, result extraction
│       └── core/                         # Nodes, FitParams, KalmanCalculator, ...
│
├── k4DTFGaudi/                           # Gaudi algorithm (k4FWCore::MultiTransformer)
├── k4DTFFCCAnalyses/                     # FCCAnalyses / RDataFrame
│   ├── include/k4DTFFCCAnalyses/
│   │   ├── DTFInput.h                    # Candidate building and MC truth matching
│   │   ├── DTFFunctions.h                # DTFRunner, fit* extractors, DTF_* functions
│   │   └── LinkDef.h
│   └── src/
│
├── tools/k4dtf_validate.cpp              # k4dtf-validate: fit + pulls on a stage-1 ntuple
├── examples/
│   ├── analysis_B0_to_Kpi.py
│   └── dtf_test_samples/                 # Generation, FCCAnalyses stage 1, docs
│
└── test/
    ├── test_core.cpp                     # Parser, table, conversions, fitter API
    ├── test_toy.cpp                      # Helix conventions; toy-MC pulls and chi2
    ├── test_input.cpp                    # Candidate building, truth matching, DTFRunner
    ├── ToyDecays.h, SyntheticEvent.h
    └── test_fccanalyses_pipeline.py
```

---

## Dependencies

| Package | Required by | Notes |
|---|---|---|
| ROOT (MathCore, Matrix, Physics, ROOTVecOps) | Core, FCCAnalyses wrapper | >= 6.28 |
| EDM4HEP | Core, Gaudi wrapper, FCCAnalyses wrapper | |
| podio | Gaudi wrapper, FCCAnalyses wrapper | |
| Gaudi | Gaudi wrapper | Required |
| k4FWCore | Gaudi wrapper | Required |
| Catch2 >= 3 | Tests | Optional |

All dependencies are provided by the Key4hep software stack — no manual
installation is required.

---

## Development

### 1. Source the Key4hep nightlies

The nightlies ship a complete, self-consistent set of all dependencies.
Source the stack **before** every build session:

```bash
# Latest nightly (rolling)
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh

# Or pin to a specific date for reproducibility
source /cvmfs/sw-nightlies.hsf.org/key4hep/releases/2026-02-26/x86_64-almalinux9-gcc14.2.0-opt/key4hep-stack/*/setup.sh
```

> **Note:** CVMFS must be mounted (`/cvmfs/sw-nightlies.hsf.org`).
> On a machine without CVMFS, use a Key4hep container or Spack instead.

### 2. Build and install

```bash
git clone https://github.com/HEP-FCC/k4DTF.git
cd k4DTF
cmake -B build -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j$(nproc)
cmake --install build
```

To disable the test suite:

```bash
cmake -B build -DCMAKE_INSTALL_PREFIX=$PWD/install -DBUILD_TESTING=OFF
```

All three sub-packages (`k4DTF`, `k4DTFGaudi`, `k4DTFFCCAnalyses`) are always
built. Gaudi and k4FWCore are required dependencies.

### 3. Activate the local build with `k4_local_repo`

After installing, run `k4_local_repo` from the repository root to prepend
the local `install/` tree to all relevant environment variables
(`PATH`, `LD_LIBRARY_PATH`, `CMAKE_PREFIX_PATH`, `ROOT_INCLUDE_PATH`, …):

```bash
# From the repository root (expects install/ to exist there)
k4_local_repo
```

The shell function is provided by the Key4hep setup script sourced in step 1.
It reads the `install/` directory relative to the current working directory,
so always call it from the repository root.

To verify the library is visible:

```bash
root -l -e 'gSystem->Load("libk4DTFFCCAnalyses"); std::cout << "OK\n";'
```

### 4. Run the tests

Run the full test suite from the build directory:

```bash
ctest --test-dir build --output-on-failure
```

Run only the C++ unit tests (requires Catch2):

```bash
ctest --test-dir build -R test_k4DTF --output-on-failure
```

Run only the FCCAnalyses RDataFrame pipeline test:

```bash
ctest --test-dir build -R test_fccanalyses_pipeline --output-on-failure
```

The pipeline test exercises the full `DTF_chi2` / `DTF_prob` / `DTF_mass` /
`DTF_fittedParticles` chain with synthetic in-memory data — no input file
required. It can also be run directly with Python:

```bash
LD_LIBRARY_PATH=$PWD/build/k4DTFFCCAnalyses:$PWD/build/k4DTF:$LD_LIBRARY_PATH \
    python test/test_fccanalyses_pipeline.py
```

To run with verbose CTest output and parallel jobs:

```bash
ctest --test-dir build --output-on-failure -j$(nproc) -V
```

---

## Usage

### Core C++ library

The fitter takes the decay tree as a `k4DTF::DecayNodeInput` (one node per
particle, measurements on the leaves) and an optional interaction point.
From a decay descriptor and the measured final-state particles,
`k4DTF::makeDecayNodeInput` builds it:

```cpp
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/TreeInput.h"

auto desc = k4DTF::DecayDescriptor::fromString("B- -> (D0 -> K- pi+) K-");

// one DTFParticle per leaf, in depth-first order: K-, pi+, K-
// charged particles with their track (helix), e.g.
//   k4DTF::DTFParticle::fromEDM4hep(recoParticle, trackState, Bz)
std::vector<k4DTF::DTFParticle> leaves = {kaon1, pion, kaon2};

auto input = k4DTF::makeDecayNodeInput(desc, leaves, /*mass constraints*/ {421});
k4DTF::InteractionPointInput pv = k4DTF::makeInteractionPoint(primaryVertex);

const k4DTF::DecayTreeFitter dtf;
k4DTF::FitResult r = dtf.fit(input, pv);

if (r.converged()) {
  std::cout << "chi2 / ndf = " << r.chi2 << " / " << r.ndf << ", prob " << r.prob() << "\n";
  const auto& D0 = r.nodes[1];  // nodes in descriptor pre-order: B-, D0, K-, pi+, K-
  std::cout << "D0 ctau = " << D0.ctau << " +- " << D0.ctauErr << " mm\n";
}
```

- Descriptors may nest decays in parentheses and use LHCb, EvtGen or ROOT
  particle names (`k4DTF::ParticleTable`).
- Composites that decay strongly or electromagnetically (D*, J/psi, ...) are
  resonances: no vertex or flight length of their own. Others get a decay
  vertex and a flight length from their production vertex.
- Tracks are helix measurements (`k4DTF::HelixMeasurement`, conventions in
  `k4DTF/core/Helix.h`, which are those of k4SimDelphes `TrackState`s); their
  mass hypothesis is the nominal mass of the descriptor's particle.
- `FitResult::nodes` has, for every node: 4-momentum and covariance, mass and
  error, decay vertex and covariance, flight length and c*tau with errors.

### Gaudi / k4FWCore

`DTFAlgorithm` is a `k4FWCore::MultiTransformer` that reads an
`edm4hep::ReconstructedParticleCollection` with the final-state particles of
one candidate (in the depth-first leaf order of the descriptor; charged
particles with their tracks) and writes every fitted node of the tree (head
first, with daughter relations) and one `Vertex` per decay vertex.

```python
from Configurables import DTFAlgorithm

dtf = DTFAlgorithm("BDTF",
    InputParticles   = ["BDaughters"],
    OutputParticles  = ["BFitted"],
    OutputVertices   = ["BFittedVertices"],
    DecayDescriptor  = "B- -> (D0 -> K- pi+) K-",
    MassConstraints  = {"421": 0.0},   # PDG code (str) -> mass [GeV], <= 0: nominal
    Bz               = 2.0,
)
```

### FCCAnalyses / RDataFrame

[`examples/dtf_test_samples`](examples/dtf_test_samples/README.md) shows the
full workflow: `k4DTF::CandidateBuilder` builds candidates from a descriptor,
`k4DTF::TruthMatcher` matches every node to MC truth, `k4DTF::DTFRunner` fits
them with their tracks and the primary vertex, and `k4DTF::fitNode*`
functions extract per-node columns.

The earlier one-call functions `DTF_chi2`, `DTF_prob`, `DTF_mass`,
`DTF_fittedParticles` and `DTF_chi2_single` remain; they fit momenta only
(no tracks), from the 4-momentum covariance of the particles, which Delphes
leaves empty.

---

## Algorithm

The implementation follows the BaBar/LHCb DecayTreeFitter (Hulsbergen 2005):
all parameters of the tree (vertex positions, flight lengths, momenta) are
fitted at once, by applying the constraints one after the other as Kalman
filter steps, iterated until the chi2 is stable.

| Node | Parameters | Constraints |
|---|---|---|
| Interaction point (optional) | position | measured position |
| Composite with lifetime | decay vertex, flight length, (px, py, pz, E) | 4-momentum conservation; flight direction along p from the production vertex; optional mass |
| Resonance | (px, py, pz, E) | 4-momentum conservation; optional mass |
| Track | (px, py, pz) | helix at the mother's vertex = measured helix (5) |

Each iteration linearises all constraints at the result of the previous one,
starting from an inflated diagonal of the previous covariance as weak prior.
The covariance update uses the Joseph form for numerical stability.

`test/test_toy.cpp` validates the fit with a toy MC (B- -> D0 K-,
B0 -> D*- pi+): pulls of momenta, vertices and c*tau are unit Gaussians and
<chi2> = ndf.

---

## References

### Fit algorithm

- R. Moser & M. Roussarie, *"Mathematical methods for B-meson oscillation analyses"*,
  Nucl. Instrum. Meth. A **384** (1997) 491–505.
  The canonical reference for the Lagrange-multiplier decay-tree fitter formalism
  used here (constraint matrix, Newton update, chi-square expression).

### Software frameworks and data model

- [Key4hep](https://key4hep.github.io/key4hep-doc/) — umbrella software stack
  for future collider experiments.

- [EDM4hep](https://github.com/key4hep/EDM4hep) — common event data model;
  defines `ReconstructedParticle`, `Vertex`, and the 10-element upper-triangle
  covariance convention used by `DTFParticle`.

- [k4FWCore](https://github.com/key4hep/k4FWCore) — Gaudi-based framework layer
  providing the `MultiTransformer` functional algorithm pattern used by
  `DTFAlgorithm`.

- [FCCAnalyses](https://github.com/HEP-FCC/FCCAnalyses) — RDataFrame-based
  analysis framework for FCC studies; defines the `Analysis` class structure
  and `RVec`-based free-function conventions used in `DTFFunctions`.

---

## License

This project is licensed under the [Apache License, Version 2.0](LICENSE),
consistent with the rest of the Key4hep software ecosystem (k4FWCore, EDM4hep,
podio, FCCAnalyses).

Source files should carry the following SPDX header:

```cpp
// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0
```

## Contributing

Contributions are welcome. Please open an issue or pull request on GitHub.
The project follows the Key4hep coding conventions and uses `clang-format`
(LLVM style) for formatting.

## Acknowledgements

The initial skeleton of this project was developed with the assistance of
[Claude](https://www.anthropic.com/claude) (Anthropic).
