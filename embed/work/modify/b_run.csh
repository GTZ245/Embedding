#!/bin/csh -f

# HTCondor invokes this script with one process index. Re-enter the same
# script inside the SL7 container so no additional wrapper scripts are needed.
set inside_sl7 = 0
if ($#argv == 2) then
  if ("$1" == "--inside-sl7") set inside_sl7 = 1
endif

if ($inside_sl7 == 0) then
  if ($#argv != 1) then
    echo "Usage: $0 <index>"
    exit 1
  endif

  echo "$1" | grep -Eq '^[0-9]+$'
  if ($status != 0) then
    echo "Error: index must be a non-negative integer"
    exit 1
  endif

  source /opt/phenix/core/bin/phenix_alma9_setup.csh -n
  sl7cmd ./b_run.csh --inside-sl7 "$1"
  exit $status
endif

set idx = "$2"

if (-f /opt/phenix/core/bin/phenix_setup.csh) then
  source /opt/phenix/core/bin/phenix_setup.csh -n
endif

if (-f /opt/phenix/core/bin/phenix_alma9_setup.csh) then
  source /opt/phenix/core/bin/phenix_alma9_setup.csh -n
endif

setenv MYLOG /gpfs/mnt/gpfs02/phenix/plhf/plhf1/$USER/log
setenv MYINSTALL /direct/phenix+u/$USER/install
setenv MYCVS /direct/phenix+u/$USER/cvsdir
setenv PATH .:$HOME/bin:$MYINSTALL/bin:$PATH
setenv LD_LIBRARY_PATH .:$MYINSTALL/lib:$LD_LIBRARY_PATH
setenv TSEARCHPATH $MYINSTALL

set tag = `printf "%05d" "$idx"`
set nevents = 900
@ real_skip = $idx * $nevents

set mcinput = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/more_condor/pack/dsts/dst_out_pi0_${tag}.root"
set oscarsinput = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/more_condor/pack/oscars/${tag}.oscar.particles.dat"
set ntananame = "ana_${tag}.root"
set ntname = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/work/output_kek/kek_${tag}.root"

echo "Running tag = ${tag}"
echo "Events = ${nevents}"
echo "Real-data skip = ${real_skip}"
echo "MC input = ${mcinput}"
echo "OSCAR input = ${oscarsinput}"
echo "Analysis output = ${ntananame}"
echo "Evaluation output = ${ntname}"

set root_call = "Fun4All_embedeval_svx.C(${nevents},\"/phenix/plhf/tongzhouguo/yuri_embed/real/work/CNTmerge_MB-0000372586-0011.root\",\"${mcinput}\",\"/phenix/hhj/lebedev/chi_c/simulation/pairobj/pairobject_chisigembed_251500_100.root\",111,\"${ntananame}\",\"/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt\",\"${oscarsinput}\",372586,\"${ntname}\",${real_skip})"

root -l -b -q "$root_call"
set root_status = $status

if ($root_status != 0) then
  echo "ROOT failed with status ${root_status}"
endif

exit $root_status
