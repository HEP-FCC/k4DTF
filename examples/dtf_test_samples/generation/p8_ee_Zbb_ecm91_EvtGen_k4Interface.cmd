! Pythia8 card for e+e- -> Z -> bbbar at 91.188 GeV, for use with
!   DelphesPythia8EvtGen_EDM4HEP_k4Interface  (11-argument, signal-forcing mode)

! 1) Settings used in the main program.
Main:timesAllowErrors = 5          ! how many aborts before run stops

! 2) Settings related to output in init(), next() and stat().
Init:showChangedSettings = on      ! list changed settings
Init:showChangedParticleData = off ! list changed particle data
Next:numberCount = 100             ! print message every n events
Next:numberShowInfo = 1            ! print event information n times
Next:numberShowProcess = 1         ! print process record n times
Next:numberShowEvent = 0           ! print event record n times

! 3) Beam parameter settings.
Beams:idA = 11                     ! first beam, e-
Beams:idB = -11                    ! second beam, e+
Beams:allowMomentumSpread = off

! Vertex smearing (FCC-ee Z-pole luminous region):
Beams:allowVertexSpread = on
Beams:sigmaVertexX = 5.96e-3
Beams:sigmaVertexY = 23.8E-6
Beams:sigmaVertexZ = 0.397
Beams:sigmaTime = 10.89            ! 36.3 ps

! 4) Hard process: Z -> bbbar at Ecm = 91.188 GeV
Beams:eCM = 91.188
WeakSingleBoson:ffbar2gmZ = on
23:onMode = off
23:onIfAny = 5
PartonLevel:ISR = on               ! initial-state radiation
PartonLevel:FSR = on               ! final-state radiation

! 5) Stop at parton level: hadronization is done by the EvtGen interface
HadronLevel:all = off

! 6) Leave all b hadrons undecayed in Pythia: EvtGen decays them
511:onMode = off
521:onMode = off
531:onMode = off
541:onMode = off
5122:onMode = off
5132:onMode = off
5142:onMode = off
5232:onMode = off
5242:onMode = off
5332:onMode = off
5342:onMode = off
5412:onMode = off
5414:onMode = off
5422:onMode = off
5424:onMode = off
5432:onMode = off
5434:onMode = off
5442:onMode = off
5444:onMode = off
5512:onMode = off
5514:onMode = off
5522:onMode = off
5524:onMode = off
5532:onMode = off
5534:onMode = off
5542:onMode = off
5544:onMode = off
