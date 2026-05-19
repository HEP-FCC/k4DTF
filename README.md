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
- [Implementing the algorithm](#implementing-the-algorithm)
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
approach of Moser & Roussarie (*NIM* A 384, 1997).

The fit output includes:

- Fitted four-momenta and updated covariances for all particles in the decay
  tree.
- A chi-square value and number of degrees of freedom.
- A fitted production vertex (when vertex constraints are enabled).

---

## Repository structure

```
k4DTF/
├── CMakeLists.txt                        # Top-level build
├── cmake/
│   └── k4DTFConfig.cmake.in             # Makes k4DTF find_package()-able
│
├── k4DTF/                               # Core library
│   ├── include/k4DTF/
│   │   ├── DTFParticle.h                # Internal particle (p4 + covariance + EDM4hep conversions)
│   │   ├── DecayDescriptor.h            # Decay topology parser
│   │   ├── FitResult.h                  # Fit output structure
│   │   └── DecayTreeFitter.h            # Main fitter class
│   └── src/
│       ├── DTFParticle.cpp              # fromEDM4hep / toEDM4hep implementations
│       ├── DecayDescriptor.cpp
│       └── DecayTreeFitter.cpp          # ← implement the fit here
│
├── k4DTFGaudi/                          # Gaudi plugin (requires k4FWCore)
│   ├── include/k4DTFGaudi/
│   │   └── DTFAlgorithm.h               # k4FWCore::MultiTransformer wrapper
│   └── src/
│       └── DTFAlgorithm.cpp
│
├── k4DTFFCCAnalyses/                    # FCCAnalyses / RDataFrame wrapper
│   ├── include/k4DTFFCCAnalyses/
│   │   ├── DTFFunctions.h               # RVec-returning free functions
│   │   └── LinkDef.h                    # ROOT dictionary directives
│   └── src/
│       └── DTFFunctions.cpp
│
├── examples/
│   └── analysis_B0_to_Kpi.py           # FCCAnalyses analysis example
│
└── test/
    ├── CMakeLists.txt                   # Catch2 + CTest
    ├── test_core.cpp                    # C++ unit tests
    └── test_fccanalyses_pipeline.py     # RDataFrame pipeline test
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

```cpp
#include "k4DTF/DecayDescriptor.h"
#include "k4DTF/DecayTreeFitter.h"
#include "k4DTF/DTFParticle.h"

// Build the decay descriptor
auto desc = k4DTF::DecayDescriptor::fromString("B0 -> K+ pi-");

// Configure the fitter
k4DTF::DecayTreeFitter dtf(std::move(desc));
dtf.addMassConstraint(511, 5.27963)   // B0 mass in GeV
   .setMaxIterations(100)
   .setConvergenceCriterion(1e-6);

// Provide measured final-state particles (order must match the descriptor)
std::vector<k4DTF::DTFParticle> daughters = { kaon, pion };

k4DTF::FitResult result = dtf.fit(daughters);

if (result.converged()) {
    std::cout << "chi2 / ndf = " << result.chi2 << " / " << result.ndf << "\n";
    std::cout << "fit prob   = " << result.prob() << "\n";
}
```

`DTFParticle` holds a `ROOT::Math::PxPyPzEVector` and a 10-element
upper-triangle covariance over `(px, py, pz, E)`, matching the layout of
`edm4hep::ReconstructedParticle::covMatrix`.

### Gaudi / k4FWCore

`DTFAlgorithm` is a `k4FWCore::MultiTransformer` that reads an
`edm4hep::ReconstructedParticleCollection` of final-state daughters and writes
back a fitted `ReconstructedParticleCollection` and a `VertexCollection`.

Python steering (Gaudi option file):

```python
from Configurables import DTFAlgorithm

dtf = DTFAlgorithm("B0DTF",
    InputParticles   = ["B0Daughters"],
    OutputParticles  = ["B0FittedDaughters"],
    OutputVertices   = ["B0FittedVertices"],
    DecayDescriptor  = "B0 -> K+ pi-",
    MassConstraints  = {"511": 5.27963},   # PDG code (str) -> mass [GeV]
    VertexConstraint = False,
    MaxIterations    = 100,
    ConvergenceEps   = 1e-6,
)

ApplicationMgr().TopAlg += [dtf]
```

The input collection must contain the final-state particles in the same order
as the leaves of the decay descriptor.

### FCCAnalyses / RDataFrame

Load the library in your analysis script and use the free functions as
`RDataFrame::Define` expressions:

```python
import ROOT
ROOT.gSystem.Load("libk4DTFFCCAnalyses")

# Make k4DTF functions and edm4hep types visible to the ROOT JIT compiler
ROOT.gInterpreter.Declare('#include "k4DTFFCCAnalyses/DTFFunctions.h"')
ROOT.gInterpreter.Declare('#include "edm4hep/ReconstructedParticleData.h"')

df = (
    ROOT.RDataFrame("events", "input.root")
    .Define("B0_DTF_chi2",
            "k4DTF::DTF_chi2(B0_daughters, \"B0 -> K+ pi-\", 511, 5.27963f)")
    .Define("B0_DTF_prob",
            "k4DTF::DTF_prob(B0_daughters, \"B0 -> K+ pi-\", 511, 5.27963f)")
    .Define("B0_DTF_mass",
            "k4DTF::DTF_mass(B0_daughters, \"B0 -> K+ pi-\", 511, 5.27963f)")
)
```

Or with a lambda to avoid repeating the descriptor string:

```python
df = df.Define("B0_DTF_chi2",
    """[](const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& v){
        return k4DTF::DTF_chi2(v, "B0 -> K+ pi-", 511, 5.27963f);
    }""",
    ["B0_daughters"])
```

Available functions (all in namespace `k4DTF`):

| Function | Return type | Description |
|---|---|---|
| `DTF_chi2` | `RVec<float>` | chi-square of the fit |
| `DTF_prob` | `RVec<float>` | chi-square probability (`TMath::Prob`) |
| `DTF_mass` | `RVec<float>` | fitted invariant mass of the mother |
| `DTF_fittedParticles` | `RVec<ReconstructedParticleData>` | all fitted particles |
| `DTF_chi2_single` | `float` | convenience wrapper for a single candidate |

All functions return `-1` (or an empty vector) when the fit does not converge.

---

## Algorithm

The fitter uses Lagrange multipliers to solve the constrained least-squares
problem

```
min (x - x0)^T C0^{-1} (x - x0)   subject to  h(x) = 0
```

where `x` is the state vector of four-momenta, `x0` the measured values,
`C0` their covariance, and `h(x)` the constraint vector.  The Newton update
at each iteration is:

```
lambda  = (H C H^T)^{-1} (H (x - x0) + h)
x_new   = x0 - C H^T lambda
chi2    = lambda^T (H C H^T) lambda
```

Convergence is declared when `|delta chi2| < eps` (default `1e-6`).

Supported constraints:

| Constraint | Status |
|---|---|
| Four-momentum conservation at each decay node | to implement in `buildConstraints()` |
| Mass constraint on any node | to implement in `buildConstraints()` |
| Vertex constraint (production vertex) | to implement — requires track parameters |

The two main TODOs in `k4DTF/src/DecayTreeFitter.cpp` are
`buildConstraints()` (populate `h` and the Jacobian `H`) and `iterate()`
(the Newton step shown above).

---

## Implementing the algorithm

The repository is a skeleton.  All interfaces, build system, tests, and
framework wrappers are in place; the actual math in
`k4DTF/src/DecayTreeFitter.cpp` is left as stubs.  The four functions to fill
in, in order, are:

### 1. `initialiseState()`

Build the state vector `state.x` (and its copy `state.x0`) from the measured
four-momenta of the input `DTFParticle` array, and fill `state.C` (and
`state.C0`) from the corresponding covariance matrices.

Each particle contributes four parameters `(px, py, pz, E)`, so for `n`
daughters the state vector has `4n` elements.  The 10-element upper-triangle
covariance of each particle maps to a `4×4` block on the diagonal of the full
covariance matrix (the cross-particle covariances are zero at input).

The EDM4hep upper-triangle index convention for `(px, py, pz, E)` is:

```
cov[0]  = σ²(px,px)
cov[1]  = σ(px,py)   cov[2]  = σ²(py,py)
cov[3]  = σ(px,pz)   cov[4]  = σ(py,pz)   cov[5]  = σ²(pz,pz)
cov[6]  = σ(px,E)    cov[7]  = σ(py,E)    cov[8]  = σ(pz,E)    cov[9]  = σ²(E,E)
```

### 2. `buildConstraints()`

Populate the constraint vector `state.h` and the Jacobian `state.H`
(`∂h/∂x`) for the current value of `state.x`.

**Four-momentum conservation** at a decay node with daughters `i` and `j`
gives four constraints:

```
h_μ = Σ_daughters p_μ^daughter - p_μ^mother = 0
```

Because the mother four-momentum is the sum of the daughters, this reduces to
a linear constraint — `H` has `+1` entries at the daughter columns and `-1`
(or `0`) at the mother column if the mother is also a free parameter.

**Mass constraint** on a node (PDG code `pdg`, target mass `m`) adds one
constraint:

```
h = (Σ p_μ)^2 - m^2 = E² - |p|² - m² = 0
```

The Jacobian row for this constraint is `∂h/∂p_μ = 2(Σ p_μ)` evaluated at
the current `state.x`.

Vertex constraints (production vertex position) are optional and require
track helix parameters; see the `setVertexConstraint()` method.

### 3. `iterate()`

Apply one Newton step of the Lagrange-multiplier update (Moser & Roussarie
eq. 14):

```
W      = H C H^T          (constraint covariance, size n_c × n_c)
lambda = W^{-1} (h + H (x - x0))
x_new  = x0 - C H^T lambda
chi2   = lambda^T W lambda
```

`W` is inverted with `TMatrixD::Invert()`.  After each step, call
`buildConstraints()` again with the updated `state.x` and check
`|Δchi2| < m_convergenceEps` for convergence.

### 4. `extractResult()`

Unpack `state.x` back into `FitResult::fittedParticles`: create one
`DTFParticle` per input daughter, copy the fitted `(px, py, pz, E)` from the
state vector, and copy the updated covariance from the diagonal block of
`state.C`.

Set `result.chi2`, `result.ndf` (number of constraints minus number of free
parameters), `result.nIter`, and `result.status = FitResult::Status::Success`.

### Where to look

| File | What to change |
|---|---|
| `k4DTF/src/DecayTreeFitter.cpp` | The four stub functions above |
| `k4DTF/include/k4DTF/FitResult.h` | `prob()` stub — replace `return 0.0` with `TMath::Prob(chi2, ndf)` |
| `test/test_core.cpp` | Tighten the fitter tests once results are non-trivial |
| `test/test_fccanalyses_pipeline.py` | Replace the loose `is not None` assertions with numerical checks |

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
