# SPDX-FileCopyrightText: Copyright (c) 2026 CERN
# SPDX-License-Identifier: Apache-2.0
'''
FCCAnalyses stage 1 for k4DTF tests: build decay-tree candidates, attach MC
truth for every node of the tree, and run k4DTF on each candidate.

The output has one entry per event (events with at least one candidate) and
RVec branches with one element per candidate; see ../README.md for the full
list. All fitter inputs are kept (track states with covariances, Bz, primary
vertex), so the fit can be re-run on this ntuple with tools/k4dtf-validate
without going back to the EDM4hep file.

Run (after building k4DTF and calling k4_local_repo in its source tree):

    fccanalysis run stage1_dtf_inputs.py -j 4 \\
        -i p8_ee_Zbb_ecm91_EvtGen_Bd2DstPi_seed1.root \\
        -o dtf_inputs_Bd2DstPi.root -- --mode Bd2DstPi

Extra options after "--": --mode, --max-candidates, --no-cc, --keep-empty,
--signal-only, --photon-errors

With --signal-only, candidates are built only from reconstructed particles
that come from a decay of the head particle (by MC truth), and only fully
truth-matched candidates are kept and fitted. This is much faster for
high-multiplicity modes, but gives no background candidates: use it for
resolution and pull studies, not for efficiencies or background.
'''

from argparse import ArgumentParser

import ROOT

# ---------------------------------------------------------------------------
# Decay modes
# ---------------------------------------------------------------------------
# labels: one name per node of the descriptor in depth-first pre-order (node 0
#         is the head); they become branch-name prefixes.
# mass_windows: |PDG| -> (min, max) [GeV] on the invariant mass of a node,
#         applied while building candidates (leaf masses from the descriptor).
# fits:   name -> (PDG codes to mass-constrain, primary-vertex constraint);
#         constraints apply to the charge conjugates too.
# derived: extra per-candidate columns.
MODES = {
    'Bd2DstPi': {
        'descriptor': 'B0 -> (D*(2010)- -> (D~0 -> K+ pi-) pi-) pi+',
        'labels': ['B', 'Dst', 'D0', 'D0_K', 'D0_pi', 'Dst_pi', 'B_pi'],
        'mass_windows': {511: (4.9, 5.7), 413: (1.95, 2.07), 421: (1.80, 1.93)},
        'fits': {'DTF': ([], False), 'DTFPV': ([421], True)},
        'derived': {'deltaM': 'Dst_M - D0_M'},
    },
    'Bu2D0K_D02KPi': {
        'descriptor': 'B- -> (D0 -> K- pi+) K-',
        'labels': ['B', 'D0', 'D0_K', 'D0_pi', 'B_K'],
        'mass_windows': {521: (4.9, 5.7), 421: (1.80, 1.93)},
        'fits': {'DTF': ([], False), 'DTFPV': ([421], True)},
        'derived': {},
    },
    'Bu2D0K_D02KSPiPi': {
        'descriptor': 'B- -> (D0 -> (KS0 -> pi+ pi-) pi+ pi-) K-',
        'labels': ['B', 'D0', 'KS', 'KS_pip', 'KS_pim', 'D0_pip', 'D0_pim', 'B_K'],
        'mass_windows': {521: (4.9, 5.7), 421: (1.80, 1.93), 310 : (0.400,0.600) },
        'fits': {'DTF': ([], False), 'DTFPV': ([421,310], True)},
        'derived': {},
    },
    'Bu2D0Pi_D02KPiPi0': {
        'descriptor': 'B- -> (D0 -> K- pi+ (pi0 -> gamma gamma)) pi-',
        'labels': ['B', 'D0', 'D0_K', 'D0_pi', 'D0_pi0', 'pi0_g1', 'pi0_g2', 'B_pi'],
        'mass_windows': {521: (4.8, 5.8), 421: (1.70, 2.00), 111: (0.09, 0.18)},
        'fits': {'DTF': ([], False), 'DTFPV': ([421, 111], True)},
        'derived': {},
    },
    'Bd2D0Pi0_D02KK': {
        'descriptor': 'B0 -> (D~0 -> K- K+) (pi0 -> gamma gamma)',
        'labels': ['B', 'D0', 'D0_Km', 'D0_Kp', 'B_pi0', 'pi0_g1', 'pi0_g2'],
        'mass_windows': {511: (4.9, 5.7), 421: (1.70, 2.00), 111: (0.09, 0.18)},
        'fits': {'DTF': ([], False), 'DTFPV': ([511,421, 111], True)},
        'derived': {},
    },    
}

# Photon resolution for clusters written without errors (k4SimDelphes before
# v00-09) or without a position: k4DTF::CalorimeterResolution. The defaults
# match the Delphes IDEA card (crystal ECAL, 0.02 x 0.02 towers with
# SmearTowerCenter, front face at R = 2250 mm, |z| = 2500 mm); change them
# here for another card.
CALORIMETER = {
    'towerSize': 0.02,
    'smearTowerCenter': True,
    'stochastic': 0.03,
    'noise': 0.002,
    'constant': 0.005,
    'radius': 2250.,
    'halfLength': 2500.,
 }

# Luminous region [um] for the primary-vertex fit; must match the Pythia card.
BEAMSPOT_UM = (4.5, 20e-3, 300.)

# Column names in EDM4hep v1 files from k4SimDelphes, with older spellings
# as fallbacks. The first name found in the input is used.
COLUMNS = {
    'k4dtf_linkFrom': ['_RecoMCLink_from.index', '_MCRecoAssociations_rec.index',
                       'MCRecoAssociations#0.index'],
    'k4dtf_linkTo': ['_RecoMCLink_to.index', '_MCRecoAssociations_sim.index',
                     'MCRecoAssociations#1.index'],
    'k4dtf_mcParents': ['_Particle_parents.index', 'Particle#0.index'],
    'k4dtf_mcDaughters': ['_Particle_daughters.index', 'Particle#1.index'],
    'k4dtf_rpTracks': ['_ReconstructedParticles_tracks.index', 'ReconstructedParticles#1.index'],
    'k4dtf_trackStates': ['_EFlowTrack_trackStates', 'EFlowTrack_1'],
}
# Only needed for modes with photons.
PHOTON_COLUMNS = {
    'k4dtf_rpClusters': ['_ReconstructedParticles_clusters.index', 'ReconstructedParticles#0.index'],
    'k4dtf_photonClusters': ['EFlowPhoton'],
}

def load_k4dtf():
    '''Load the k4DTF FCCAnalyses library and its headers (once).'''
    if ROOT.gSystem.Load('libk4DTFFCCAnalyses') < 0:
        raise RuntimeError('Cannot load libk4DTFFCCAnalyses: build k4DTF and '
                           'run k4_local_repo in its source tree first.')
    ROOT.gInterpreter.Declare('#include "k4DTFFCCAnalyses/DTFFunctions.h"')


def cpp_windows(windows):
    '''{511: (4.9, 5.7)} -> C++ initialiser of CandidateBuilder::MassWindows.'''
    return '{' + ', '.join(f'{{{pdg}, {{{lo}, {hi}}}}}' for pdg, (lo, hi) in windows.items()) + '}'


def declare_config(mode, cfg, max_candidates, include_cc):
    '''Declare the candidate builder, truth matcher and fitters as C++ globals.

    They are constructed once (not per event) and used from the RDataFrame
    graph; their call operators are const, so this is safe with several threads.
    '''
    ns = f'k4dtf_{mode}'
    fits = '\n'.join(
        f'const k4DTF::DTFRunner fit_{name}(builder, {{{", ".join(map(str, pdgs))}}}, '
        f'{"true" if pv else "false"});'
        for name, (pdgs, pv) in cfg['fits'].items())
    calo = '\n'.join(f'  c.{k} = {str(v).lower() if isinstance(v, bool) else v};'
                     for k, v in CALORIMETER.items())    
    ROOT.gInterpreter.Declare(f'''
#ifndef K4DTF_CFG_{mode.upper()}
#define K4DTF_CFG_{mode.upper()}
namespace {ns} {{
const k4DTF::CalorimeterResolution calo = [] {{
  k4DTF::CalorimeterResolution c;
{calo}
  return c;
}}();
const k4DTF::CandidateBuilder builder("{cfg["descriptor"]}", {cpp_windows(cfg["mass_windows"])},
                                      {"true" if include_cc else "false"}, {max_candidates});
const k4DTF::TruthMatcher matcher(builder);
{fits}
}}
#endif
''')
    return ns


def describe_nodes(cfg):
    '''Return [(label, node index, leaf index or -1)] in pre-order.'''
    desc = ROOT.k4DTF.DecayDescriptor.fromString(cfg['descriptor'])
    flat = desc.flatten()
    if flat.size() != len(cfg['labels']):
        raise ValueError(f'{len(cfg["labels"])} labels given for {flat.size()} nodes of '
                         f'"{desc.toString()}"')
    return [(label, i, flat[i].leafIndex) for i, label in enumerate(cfg['labels'])]


class Analysis():
    '''Decay-tree candidates, MC truth and k4DTF fit results.'''

    def __init__(self, cmdline_args):
        parser = ArgumentParser(description='k4DTF test-sample options',
                                usage='Provided after "--"')
        parser.add_argument('--mode', default='Bd2DstPi', choices=sorted(MODES),
                            help='Decay mode (default: Bd2DstPi)')
        parser.add_argument('--max-candidates', default=10000, type=int,
                            help='Maximum number of candidates per event and mode')
        parser.add_argument('--no-cc', action='store_true',
                            help='Do not build the charge-conjugate mode')
        parser.add_argument('--keep-empty', action='store_true',
                            help='Keep events without candidates')
        parser.add_argument('--photon-errors', default='file', choices=['file', 'model'],
                            help='Photon position/energy errors from the file when it has them (default), '
                                 'or always from CALORIMETER')        
        parser.add_argument('--signal-only', action='store_true',
                            help='Build candidates only from particles of the signal decay (MC truth) '
                                 'and keep only truth-matched candidates; much faster, no background')
        self.ana_args, _ = parser.parse_known_args(cmdline_args.get('remaining', []))

        self.cfg = MODES[self.ana_args.mode]
        load_k4dtf()
        self.ns = declare_config(self.ana_args.mode, self.cfg, self.ana_args.max_candidates,
                                 not self.ana_args.no_cc)
        self.nodes = describe_nodes(self.cfg)
        flat = ROOT.k4DTF.DecayDescriptor.fromString(self.cfg['descriptor']).flatten()
        self.pdgs = [abs(int(flat[n].pdg)) for n in range(flat.size())]
        self.head_pdg = int(flat[0].pdg)
        self.has_photons = any(self.pdgs[n] == 22 for _, n, leaf in self.nodes if leaf >= 0)        


    def _plan(self, bz_expr='2.0f'):
        '''Ordered list of (column, expression, written to output).

        A column named FILTER is a Filter on its expression. Shared by
        analyzers() and output(), since FCCAnalyses asks for the output list
        before building the graph.
        '''
        ns, cfg = self.ns, self.cfg
        bsx, bsy, bsz = BEAMSPOT_UM
        plan = [
            # per reconstructed particle
            ('Bz', bz_expr, True),
            ('RP_mcIdx', 'k4DTF::recoToMC(ReconstructedParticles, k4dtf_linkFrom, k4dtf_linkTo)', False),
            ('RP_pid', 'k4DTF::truthPDG(RP_mcIdx, Particle)', False),
            ('RP_ts', 'k4DTF::trackStatePerReco(ReconstructedParticles, k4dtf_rpTracks, EFlowTrack, '
                      'k4dtf_trackStates)', False),
        ]
        if self.has_photons:
            # photon clusters with position and errors, filled in from the
            # calorimeter model where the file has none
            plan.append(('RP_cluster', 'k4DTF::photonClusterPerReco(ReconstructedParticles, k4dtf_rpClusters, '
                                       f'k4dtf_photonClusters, {ns}::calo, '
                                       f'{"true" if self.ana_args.photon_errors == "model" else "false"})', False))        
        # candidates and truth matching
        matcher = f'{ns}::matcher({{}}, RP_mcIdx, Particle, k4dtf_mcParents, k4dtf_mcDaughters)'
        if self.ana_args.signal_only:
            # Only particles from a decay of the head particle enter the
            # combinatorics; of their combinations only the truth-matched
            # ones are kept (and fitted).
            plan += [
                ('RP_pidSignal', f'k4DTF::signalPID(RP_pid, RP_mcIdx, Particle, k4dtf_mcParents, '
                                 f'{self.head_pdg})', False),
                ('cands_built', f'{ns}::builder(ReconstructedParticles, RP_pidSignal)', False),
                ('truth_built', matcher.format('cands_built'), False),
                ('cands', 'k4DTF::selectCandidates(cands_built, truth_built.isSignal)', False),
                ('truth', 'k4DTF::selectCandidates(truth_built, truth_built.isSignal)', False),
                ('nCand', 'k4DTF::nCandidates(cands)', True),
            ]
            if not self.ana_args.keep_empty:
                plan.append(('FILTER', 'nCand > 0', False))
        else:
            plan.append(('cands', f'{ns}::builder(ReconstructedParticles, RP_pid)', False))
            plan.append(('nCand', 'k4DTF::nCandidates(cands)', True))
            if not self.ana_args.keep_empty:
                plan.append(('FILTER', 'nCand > 0', False))
            plan.append(('truth', matcher.format('cands'), False))
        plan += [
            ('cand_isCC', 'cands.isCC', True),
            ('cand_isSignal', 'truth.isSignal', True),
            ('cand_nMatched', 'truth.nMatched', True),
            # primary vertex (FCCAnalyses fit with beam-spot constraint)
            ('PV_tracks', 'FCCAnalyses::VertexFitterSimple::get_PrimaryTracks('
                          f'k4dtf_trackStates, true, {bsx}, {bsy}, {bsz}, 0., 0., 0.)', False),
            ('PV_object', 'FCCAnalyses::VertexFitterSimple::VertexFitter_Tk('
                          f'1, PV_tracks, true, {bsx}, {bsy}, {bsz}, 0., 0., 0.)', False),
            ('PV_reco', 'FCCAnalyses::VertexingUtils::get_VertexData(PV_object)', True),
            ('PV_nTracks', 'PV_object.ntracks', True),
        ]
        plan += [(f'PV_true_{a}', f'k4DTF::truePrimaryVertex(Particle, {k})', True)
                 for k, a in enumerate('xyz')]
        plan.append(('cand_nLeavesInPV', 'k4DTF::nLeavesInPV(cands, RP_ts, PV_tracks)', True))

        # leaves and composite nodes
        for label, node, leaf in self.nodes:
            plan.append((f'{label}_mc', f'k4DTF::nodeMCParticle(truth, Particle, {node})', True))
            if leaf >= 0:
                plan += [
                    (f'{label}_recoIdx', f'k4DTF::leafRecoIndex(cands, {leaf})', True),
                    (f'{label}_reco', f'k4DTF::leafParticles(cands, ReconstructedParticles, RP_ts, '
                                      f'{ns}::builder, {leaf}, Bz)', True),
                    (f'{label}_track', f'k4DTF::leafTrackStates(cands, RP_ts, {leaf})', True),
                    (f'{label}_trueTrack', f'k4DTF::leafTrueTrackState(truth, Particle, {ns}::builder, '
                                           f'{leaf}, Bz)', True),
                ]
                if self.has_photons and self.pdgs[node] == 22:
                    plan.append((f'{label}_cluster', f'k4DTF::leafClusters(cands, RP_cluster, {leaf})', True))                
            else:
                plan += [
                    (f'{label}_M', f'k4DTF::nodeMass(cands, ReconstructedParticles, {ns}::builder, {node})',
                     True),
                    (f'{label}_trueCtau', f'k4DTF::nodeTrueCtau(truth, Particle, k4dtf_mcDaughters, {node})',
                     True),
                ]
                for k, a in enumerate('xyz'):
                    plan += [
                        (f'{label}_trueProdVtx_{a}',
                         f'k4DTF::nodeTrueProductionVertex(truth, Particle, {node}, {k})', True),
                        (f'{label}_trueDecayVtx_{a}',
                         f'k4DTF::nodeTrueDecayVertex(truth, Particle, k4dtf_mcDaughters, {node}, {k})', True),
                    ]
        plan += [(name, expr, True) for name, expr in cfg['derived'].items()]

        # decay-tree fits
        for fit in cfg['fits']:
            fits = f'{fit}_fits'
            clusters = 'RP_cluster, ' if self.has_photons else ''
            plan.append((fits, f'{ns}::fit_{fit}(cands, ReconstructedParticles, RP_ts, {clusters}Bz, PV_reco)',
                         False))
            plan += [(f'{fit}_{col}', f'k4DTF::{func}({fits})', True)
                     for col, func in [('status', 'fitStatus'), ('chi2', 'fitChi2'), ('ndf', 'fitNdf'),
                                       ('nIter', 'fitNIter'), ('prob', 'fitProb')]]
            for label, node, leaf in self.nodes:
                p = f'{fit}_{label}'
                plan += [(f'{p}_M', f'k4DTF::fitNodeMass({fits}, {node})', True),
                         (f'{p}_MErr', f'k4DTF::fitNodeMassErr({fits}, {node})', True)]
                plan += [(f'{p}_p{a}', f'k4DTF::fitNodeMomentum({fits}, {node}, {k})', True)
                         for k, a in enumerate('xyz')]
                if leaf < 0:
                    for k, a in enumerate('xyz'):
                        plan += [(f'{p}_vtx_{a}', f'k4DTF::fitNodeDecayVertex({fits}, {node}, {k})', True),
                                 (f'{p}_vtxErr_{a}', f'k4DTF::fitNodeDecayVertexErr({fits}, {node}, {k})',
                                  True)]
                    plan += [(f'{p}_{col}', f'k4DTF::{func}({fits}, {node})', True)
                             for col, func in [('ctau', 'fitNodeCtau'), ('ctauErr', 'fitNodeCtauErr'),
                                               ('L', 'fitNodeDecayLength'), ('LErr', 'fitNodeDecayLengthErr')]]
        return plan

    def analyzers(self, dframe):
        '''Define candidates, truth and fit columns.'''
        d = dframe

        # input columns: EDM4hep v1 names, older ones as fallback
        names = set(str(c) for c in d.GetColumnNames())

        def has(col):
            return col in names or d.HasColumn(col)

        columns = dict(COLUMNS, **(PHOTON_COLUMNS if self.has_photons else {}))
        for alias, choices in columns.items():            
            found = next((c for c in choices if has(c)), None)
            if found is None:
                raise RuntimeError(f'None of {choices} found in the input; check the file with podio-dump.')
            d = d.Alias(alias, found)
        bz = 'k4DTF::firstOr(magFieldBz, 2.0f)' if has('magFieldBz') else '2.0f'

        for name, expr, _ in self._plan(bz):
            d = d.Filter(expr, expr) if name == 'FILTER' else d.Define(name, expr)
        return d

    def output(self):
        '''Columns written to the output file.'''
        return [name for name, _, keep in self._plan() if keep]
