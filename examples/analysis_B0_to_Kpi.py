# SPDX-FileCopyrightText: Copyright (c) 2026 CERN
# SPDX-License-Identifier: Apache-2.0

"""
FCCAnalyses analysis: B0 -> K+ pi- with decay-tree fit

Selects B0 -> K+ pi- candidates from an EDM4hep input file, applies the
k4DTF decay-tree fitter with a B0 mass constraint, and writes the fitted
quantities to an output ROOT file.

Run with:
    fccanalysis run examples/analysis_B0_to_Kpi.py
"""

import ROOT

# ---------------------------------------------------------------------------
# Load k4DTF — adjust the path if k4_local_repo has not been called
# ---------------------------------------------------------------------------
ROOT.gSystem.Load("libk4DTFFCCAnalyses")

ROOT.gInterpreter.Declare('#include "edm4hep/ReconstructedParticleData.h"')

# ---------------------------------------------------------------------------
# Physics constants
# ---------------------------------------------------------------------------
B0_PDG   = 511
B0_MASS  = 5.27963   # GeV
KAON_PDG = 321       # K+
PION_PDG = -211      # pi-

# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------
class Analysis:

    def analysers(df):
        df2 = (
            df

            # ------------------------------------------------------------------
            # 1. Select kaon and pion candidates from the reconstructed
            #    particle collection.  In a real analysis these would come from
            #    a PID algorithm; here we select by PDG code as a placeholder.
            # ------------------------------------------------------------------
            .Define("Kaons",
                    f"ROOT::VecOps::Filter(ReconstructedParticles,"
                    f" [](const edm4hep::ReconstructedParticleData& p)"
                    f"{{ return p.PDG == {KAON_PDG}; }})")
            .Define("Pions",
                    f"ROOT::VecOps::Filter(ReconstructedParticles,"
                    f" [](const edm4hep::ReconstructedParticleData& p)"
                    f"{{ return p.PDG == {PION_PDG}; }})")

            # ------------------------------------------------------------------
            # 2. Form B0 candidates: concatenate [K+, pi-] into a flat
            #    daughter vector matching the order in the decay descriptor.
            # ------------------------------------------------------------------
            .Define("B0_daughters",
                    "ROOT::VecOps::Concatenate(Kaons, Pions)")

            .Filter("Kaons.size() >= 1 && Pions.size() >= 1",
                    "At least one K+ and one pi-")

            # ------------------------------------------------------------------
            # 3. Unconstrained invariant mass of the B0 candidate
            # ------------------------------------------------------------------
            .Define("B0_mass_raw",
                    """[](const ROOT::VecOps::RVec<edm4hep::ReconstructedParticleData>& v) {
                        ROOT::Math::PxPyPzEVector p4;
                        for (const auto& d : v)
                            p4 += ROOT::Math::PxPyPzEVector(
                                d.momentum.x, d.momentum.y, d.momentum.z, d.energy);
                        return static_cast<float>(p4.M());
                    }(B0_daughters)""")

            # ------------------------------------------------------------------
            # 4. Decay-tree fit with B0 mass constraint
            # ------------------------------------------------------------------
            .Define("B0_DTF_chi2",
                    f"k4DTF::DTF_chi2_single(B0_daughters,"
                    f' "B0 -> K+ pi-", {B0_PDG}, {B0_MASS}f)')

            .Define("B0_DTF_prob",
                    f"k4DTF::DTF_prob(B0_daughters,"
                    f' "B0 -> K+ pi-", {B0_PDG}, {B0_MASS}f)[0]')

            .Define("B0_DTF_mass",
                    f"k4DTF::DTF_mass(B0_daughters,"
                    f' "B0 -> K+ pi-", {B0_PDG}, {B0_MASS}f)[0]')

            .Define("B0_DTF_daughters",
                    f"k4DTF::DTF_fittedParticles(B0_daughters,"
                    f' "B0 -> K+ pi-", {B0_PDG}, {B0_MASS}f)')

            # Fitted kaon and pion momenta magnitude
            .Define("Kaon_DTF_p",
                    """B0_DTF_daughters.size() > 0
                        ? std::sqrt(std::pow(B0_DTF_daughters[0].momentum.x, 2)
                                  + std::pow(B0_DTF_daughters[0].momentum.y, 2)
                                  + std::pow(B0_DTF_daughters[0].momentum.z, 2))
                        : -999.f""")
            .Define("Pion_DTF_p",
                    """B0_DTF_daughters.size() > 1
                        ? std::sqrt(std::pow(B0_DTF_daughters[1].momentum.x, 2)
                                  + std::pow(B0_DTF_daughters[1].momentum.y, 2)
                                  + std::pow(B0_DTF_daughters[1].momentum.z, 2))
                        : -999.f""")

            # Quality selection: require fit converged (chi2 > 0) and good probability
            .Filter("B0_DTF_chi2 > 0",       "DTF converged")
            .Filter("B0_DTF_prob > 0.001",   "DTF prob > 0.1%")
        )
        return df2

    def output():
        return [
            "B0_mass_raw",
            "B0_DTF_chi2",
            "B0_DTF_prob",
            "B0_DTF_mass",
            "Kaon_DTF_p",
            "Pion_DTF_p",
        ]
