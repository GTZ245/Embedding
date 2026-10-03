#!/usr/bin/env bash
#source /opt/phenix/core/bin/phenix_setup.csh

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <index>"
  exit 1
fi

idx_raw="$1"

if ! [[ "$idx_raw" =~ ^[0-9]+$ ]]; then
  echo "Error: index must be an integer"
  exit 1
fi

idx=$((10#$idx_raw))
tag=$(printf "%05d" "$idx")

mcinput="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/more_condor/pack/dsts/dst_out_pi0_${tag}.root"
oscarsinput="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/more_condor/pack/oscars/${tag}.oscar.particles.dat"
#mcinput="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/error_files/pack/dsts/dst_out_pi0_${tag}.root"
#oscarsinput="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/error_files/pack/oscars/${tag}.oscar.particles.dat"
ntananame="ana_${tag}.root"
ntname="/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/work/output_kek/kek_${tag}.root"

echo "Running tag = ${tag}"
echo "${mcinput}"
echo "${oscarsinput}"
echo "${ntananame}"
echo "${ntname}"

root -l -b -q "Fun4All_embedeval_svx.C( \
  900, \
  \"/phenix/plhf/tongzhouguo/yuri_embed/real/work/CNTmerge_MB-0000372586-0011.root\", \
  \"${mcinput}\", \
  \"/phenix/hhj/lebedev/chi_c/simulation/pairobj/pairobject_chisigembed_251500_100.root\", \
  111, \
  \"${ntananame}\", \
  \"/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt\", \
  \"${oscarsinput}\", \
  372586, \
  \"${ntname}\" \
)"

