# SPDX-FileCopyrightText: Copyright (c) 2026 CERN
# SPDX-License-Identifier: Apache-2.0

"""
CTest script: verify the k4DTFFCCAnalyses RDataFrame pipeline runs end-to-end.

Creates a one-event in-memory RDataFrame with synthetic K+ and pi- candidates,
applies k4DTF::DTF_chi2_single / DTF_mass / DTF_fittedParticles, and checks
that the pipeline executes without error.

Exit code 0 = pass, non-zero = fail (CTest convention).
"""

import sys
import ROOT

# ---------------------------------------------------------------------------
# Load the library (LD_LIBRARY_PATH must include the build tree — CTest sets
# this via set_tests_properties ENVIRONMENT)
# ---------------------------------------------------------------------------
if ROOT.gSystem.Load("libk4DTFFCCAnalyses") < 0:
    print("ERROR: could not load libk4DTFFCCAnalyses", file=sys.stderr)
    sys.exit(1)

# Make the k4DTF namespace and edm4hep types visible to the JIT compiler
# used by RDataFrame string-based Define expressions.
ROOT.gInterpreter.Declare('#include "k4DTFFCCAnalyses/DTFFunctions.h"')
ROOT.gInterpreter.Declare('#include "edm4hep/ReconstructedParticleData.h"')

# ---------------------------------------------------------------------------
# C++ helper: build a two-daughter vector (K+ and pi-) on-shell
# ---------------------------------------------------------------------------
ROOT.gInterpreter.Declare(r"""
#include <cmath>
#include "edm4hep/ReconstructedParticleData.h"

ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> makeBDaughters() {
    ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData> v(2);

    // K+  (PDG 321, mass 493.7 MeV, pz = +2.5 GeV)
    const double mK = 0.4937, pzK = 2.5;
    v[0].PDG      = 321;
    v[0].charge   = 1.f;
    v[0].momentum = {0.f, 0.f, static_cast<float>(pzK)};
    v[0].energy   = static_cast<float>(std::sqrt(mK*mK + pzK*pzK));
    v[0].mass     = static_cast<float>(mK);

    // pi- (PDG -211, mass 139.6 MeV, pz = -2.5 GeV)
    const double mPi = 0.1396, pzPi = -2.5;
    v[1].PDG      = -211;
    v[1].charge   = -1.f;
    v[1].momentum = {0.f, 0.f, static_cast<float>(pzPi)};
    v[1].energy   = static_cast<float>(std::sqrt(mPi*mPi + pzPi*pzPi));
    v[1].mass     = static_cast<float>(mPi);

    // momentum covariance (lower triangle px, py, pz, E): sigma = 50 MeV
    for (auto& p : v)
        for (int i : {0, 2, 5, 9}) p.covMatrix[i] = 0.05f * 0.05f;

    return v;
}
""")

# ---------------------------------------------------------------------------
# Build a one-event RDataFrame and exercise the DTF functions
# ---------------------------------------------------------------------------
df = (
    ROOT.RDataFrame(1)
    .Define("B0_daughters", "makeBDaughters()")
    .Define("B0_DTF_chi2",
            'k4DTF::DTF_chi2_single(B0_daughters, "B0 -> K+ pi-", 511, 5.27963f)')
    .Define("B0_DTF_prob",
            'k4DTF::DTF_prob(B0_daughters, "B0 -> K+ pi-", 511, 5.27963f)')
    .Define("B0_DTF_mass",
            'k4DTF::DTF_mass(B0_daughters, "B0 -> K+ pi-", 511, 5.27963f)')
    .Define("B0_DTF_daughters",
            'k4DTF::DTF_fittedParticles(B0_daughters, "B0 -> K+ pi-", 511, 5.27963f)')
)

# Trigger the event loop
chi2     = df.Take["float"]("B0_DTF_chi2").GetValue()[0]
prob     = df.Take["ROOT::VecOps::RVec<float>"]("B0_DTF_prob").GetValue()[0]
mass     = df.Take["ROOT::VecOps::RVec<float>"]("B0_DTF_mass").GetValue()[0]
nfitted  = df.Take["ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>"]("B0_DTF_daughters").GetValue()[0]

print(f"DTF chi2          = {chi2}")
print(f"DTF prob          = {prob}")
print(f"DTF fitted mass   = {mass}")
print(f"DTF n fitted ptcl = {len(nfitted)}")

# ---------------------------------------------------------------------------
# Assertions
# ---------------------------------------------------------------------------
# Momentum-only fit with a B0 mass constraint: the invariant mass of the
# inputs is 5.05 GeV, so the constraint pulls it to 5.27963 with chi2 > 0.
assert chi2 > 0, f"Expected chi2 > 0, got {chi2}"
assert 0 <= prob[0] <= 1, f"Expected prob in [0,1], got {prob[0]}"
assert abs(mass[0] - 5.27963) < 1e-4, f"Fitted mass is not the constraint: {mass[0]}"
assert len(nfitted) == 3, f"Expected 3 fitted particles (B0, K+, pi-), got {len(nfitted)}"

print("PASS")
sys.exit(0)
