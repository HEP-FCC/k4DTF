! Pythia8 card for k4DTF test samples: e+e- -> Z -> b bbar at the Z pole.
! Based on FCC-config FCCee/Generator/Pythia8/p8_ee_Zbb_ecm91_EvtGen_Bu2D0Pi.cmd,
! with luminous-region smearing added so that the primary vertex is not
! exactly at the origin. generate.sh appends Main:numberOfEvents and
! Random:seed to a per-job copy of this file.

! 1) Settings used in the main program.
Random:setSeed = on
Main:timesAllowErrors = 5          ! how many aborts before run stops

! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showChangedParticleData = off ! list changed particle data
Next:numberCount = 100             ! print message every n events
Next:numberShowInfo = 1            ! print event information n times
Next:numberShowProcess = 1         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

! 3) Beam parameter settings.
Beams:idA = 11
Beams:idB = -11
Beams:eCM = 91.188                 ! CM energy of collision

! 4) Luminous region [mm]. These match the beam-spot constraint used by the
!    FCCAnalyses primary-vertex fit in the analysis script
!    (4.5 um x 20 nm x 300 um); keep the two consistent if you change them.
Beams:allowVertexSpread = on
Beams:sigmaVertexX = 4.5e-3
Beams:sigmaVertexY = 2.0e-5
Beams:sigmaVertexZ = 0.30

! 5) Hard process: Z -> b bbar
WeakSingleBoson:ffbar2gmZ = on
23:onMode = off
23:onIfAny = 5
