#!/usr/bin/env bash
# Generate k4DTF test samples with Pythia8 + EvtGen + Delphes (IDEA), EDM4hep output.
# Follows the FCC tutorial "Sample generation: Pythia8 + EvtGen".
#
# Usage:
#   source /cvmfs/sw.hsf.org/key4hep/setup.sh
#   ./generate.sh [MODE] [NEVENTS] [SEED] [OUTDIR]
#
#   MODE     Bd2DstPi (default)  B0 -> (D*- -> (D~0 -> K+ pi-) pi-) pi+          (+ c.c.)
#            Bu2D0K_D02KPi       B- -> (D0 -> K- pi+) K-                         (+ c.c.)
#            Bu2D0K_D02KSPiPi    B- -> (D0 -> KS pi+ pi-) K-                     (+ c.c.)
#            Bu2D0Pi_D02KPiPi0   B- -> (D0 -> K- pi+ (pi0 -> gamma gamma)) pi-   (+ c.c.)
#            Bd2D0Pi0_D02KK      B0 -> (D~0 -> K- K+) (pi0 -> gamma gamma)       (+ c.c.)
#   NEVENTS  number of events (default 10000)
#   SEED     Pythia8 random seed (default 1)
#   OUTDIR   output directory (default ./samples)
#
# Every event contains one signal B decaying through the chain above; the other
# b hadron decays generically.
set -euo pipefail

MODE=${1:-Bd2DstPi}
NEVENTS=${2:-10000}
SEED=${3:-1}
OUTDIR=${4:-samples}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "${MODE}" in
  Bd2DstPi) DECFILE=Bd2DstPi_D0Kpi.dec; SIGPDG=511; SIGALIAS=Bd_SIGNAL ;;
  Bu2D0K_D02KPi)   DECFILE=Bu2D0K_D0Kpi.dec;   SIGPDG=521; SIGALIAS=Bu_SIGNAL ;;
  Bu2D0K_D02KSPiPi)   DECFILE=Bu2D0K_D0KSpipi.dec;   SIGPDG=521; SIGALIAS=Bu_SIGNAL ;;
  Bu2D0Pi_D02KPiPi0)   DECFILE=Bu2D0Pi_D0Kpipi0.dec;   SIGPDG=521; SIGALIAS=Bu_SIGNAL ;;
  Bd2D0Pi0_D02KK) DECFILE=Bd2D0Pi0_D0KK.dec; SIGPDG=511; SIGALIAS=Bd_SIGNAL ;;  
  *) echo "Unknown MODE '${MODE}' (see comments in this file for options)" >&2; exit 1 ;;
esac

EXE=DelphesPythia8EvtGen_EDM4HEP_k4Interface
if ! command -v "${EXE}" > /dev/null; then
  echo "${EXE} not found: source the Key4hep stack first" >&2
  exit 1
fi

mkdir -p "${OUTDIR}"
WORK="${OUTDIR}/work_${MODE}_${SEED}"
mkdir -p "${WORK}"

# Detector card, EDM4hep output steering and EvtGen tables from FCC-config.
# The Key4hep stack ships FCC-config and exports $FCCCONFIG; fall back to GitHub.
fetch() {  # fetch <path inside FCC-config/FCCee>
  local rel=$1 name
  name=$(basename "${rel}")
  if [[ -n "${FCCCONFIG:-}" && -f "${FCCCONFIG}/FCCee/${rel}" ]]; then
    cp "${FCCCONFIG}/FCCee/${rel}" "${WORK}/${name}"
  else
    wget -q -O "${WORK}/${name}" \
      "https://raw.githubusercontent.com/HEP-FCC/FCC-config/main/FCCee/${rel}"
  fi
}
#fetch Delphes/card_IDEA.tcl
cp "${HERE}/card_IDEA.tcl" "${WORK}"
fetch Delphes/edm4hep_IDEA.tcl
fetch Generator/EvtGen/DECAY.DEC
fetch Generator/EvtGen/evt.pdl
cp "${HERE}/${DECFILE}" "${WORK}/"

CARD="${WORK}/p8_ee_Zbb_ecm91_EvtGen_${MODE}.cmd"
#cp "${HERE}/p8_ee_Zbb_ecm91_EvtGen_DTF.cmd" "${CARD}"
cp "${HERE}/p8_ee_Zbb_ecm91_EvtGen_k4Interface.cmd" "${CARD}"
{
  echo ""
  echo "! added by generate.sh"
  echo "Main:numberOfEvents = ${NEVENTS}"
  echo "Random:seed = ${SEED}"
} >> "${CARD}"

OUTFILE="$(cd "${OUTDIR}" && pwd)/p8_ee_Zbb_ecm91_EvtGen_${MODE}_seed${SEED}.root"

# Argument order: detector card, EDM4hep steering, Pythia card, output,
# DECAY.DEC, evt.pdl, user decay file, signal PDG, signal alias, force flag.
(
  cd "${WORK}"
  "${EXE}" card_IDEA.tcl edm4hep_IDEA.tcl "$(basename "${CARD}")" "${OUTFILE}" \
    DECAY.DEC evt.pdl "${DECFILE}" "${SIGPDG}" "${SIGALIAS}" 1
)

echo "Wrote ${OUTFILE}"
echo "Check the collections with:  podio-dump ${OUTFILE}"
