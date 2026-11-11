// SPDX-FileCopyrightText: Copyright (c) 2026 CERN
// SPDX-License-Identifier: Apache-2.0

#ifdef __CLING__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// Functions are selected one by one rather than the whole namespace, to keep
// in-memory classes such as DTFParticle (std::optional member) out of the
// I/O dictionary. Analysis scripts declare the headers, so classes are
// available to the interpreter without dictionary entries.
#pragma link C++ function k4DTF::DTF_chi2;
#pragma link C++ function k4DTF::DTF_prob;
#pragma link C++ function k4DTF::DTF_mass;
#pragma link C++ function k4DTF::DTF_fittedParticles;
#pragma link C++ function k4DTF::DTF_chi2_single;

// candidate building, truth matching and candidate fits (DTFInput.h, DTFFunctions.h)
#pragma link C++ function k4DTF::recoToMC;
#pragma link C++ function k4DTF::truthPDG;
#pragma link C++ function k4DTF::trackStatePerReco;
#pragma link C++ function k4DTF::firstOr;
#pragma link C++ function k4DTF::signalPID;
#pragma link C++ function k4DTF::nCandidates;
#pragma link C++ function k4DTF::selectCandidates;
#pragma link C++ function k4DTF::leafRecoIndex;
#pragma link C++ function k4DTF::leafParticles;
#pragma link C++ function k4DTF::leafTrackStates;
#pragma link C++ function k4DTF::nodeMass;
#pragma link C++ function k4DTF::nLeavesInPV;
#pragma link C++ function k4DTF::nodeMCParticle;
#pragma link C++ function k4DTF::nodeTrueProductionVertex;
#pragma link C++ function k4DTF::nodeTrueDecayVertex;
#pragma link C++ function k4DTF::nodeTrueCtau;
#pragma link C++ function k4DTF::leafTrueTrackState;
#pragma link C++ function k4DTF::truePrimaryVertex;
#pragma link C++ function k4DTF::fitStatus;
#pragma link C++ function k4DTF::fitChi2;
#pragma link C++ function k4DTF::fitNdf;
#pragma link C++ function k4DTF::fitNIter;
#pragma link C++ function k4DTF::fitProb;
#pragma link C++ function k4DTF::fitNodeMass;
#pragma link C++ function k4DTF::fitNodeMassErr;
#pragma link C++ function k4DTF::fitNodeP;
#pragma link C++ function k4DTF::fitNodeMomentum;
#pragma link C++ function k4DTF::fitNodeDecayVertex;
#pragma link C++ function k4DTF::fitNodeDecayVertexErr;
#pragma link C++ function k4DTF::fitNodeCtau;
#pragma link C++ function k4DTF::fitNodeCtauErr;
#pragma link C++ function k4DTF::fitNodeDecayLength;
#pragma link C++ function k4DTF::fitNodeDecayLengthErr;

#endif
