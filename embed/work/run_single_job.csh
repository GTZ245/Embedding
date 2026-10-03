#!/usr/local/bin/tcsh -f

#Do the embedding

################################################################
# This is based on the original logic by Takashi.
#
# The idea is to set up the environment in a uniform way to facilitate
# accessing data among the various stages of the embedding process.
#
# For that reason the script embed/setup/embedding_setup template
# needs to be used (after modification appropriate for the user)

setenv prompt 1
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source $i
end
source $HOME/.login
unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv ROOTSYS

setenv DATADIR /gpfs/mnt/gpfs02/phenix/plhf/plhf1/berd/sim/output/embedding

source /opt/phenix/bin/phenix_setup.csh new
setenv LD_LIBRARY_PATH `pwd`/../lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/hhj/hhj1/hachiya/15.08/source/svx_cent_ana/install/lib:$LD_LIBRARY_PATH

##################################

if( $#argv != 1) then
	echo "run_single_job.csh requires job number as an argument to run:"
exit -1
endif

echo "Setting internal variables..."

set outputdir = "$DATADIR/Embed"
if ( ! -d $outputdir) then
	mkdir -p $outputdir
endif

set jobnum     = $1
set evtnum    = 0
set runnum = 231429
#set runnum = 240121
set particle = "akaon"
set centrality = "60-93"
#set inputreal = "$DATADIR/CNT/CNT_MB_$centrality-0000$runnum-$runnum_part.root"
set inputreal = "$DATADIR/CNT/MB_$centrality-$runnum.root"
set inputsim  = "$DATADIR/pisaToDST/${particle}_${jobnum}.root"
set outdst    = "$outputdir/outdst_${particle}_${jobnum}_${centrality}.root"
set outntuple = "$outputdir/${particle}_${jobnum}_${centrality}.root"

set scriptdir = `pwd`

set tmpdir      = "/home/tmp/${USER}_job_$jobnum-$centrality"

set inreal = `basename $inputreal`
set insim  = `basename $inputsim`

echo "runnum       $runnum     "
echo "jobnum       $jobnum     "
echo "evtnum       $evtnum     "
echo "inputreal    $inputreal  "
echo "inputsim     $inputsim   "
echo "outdst       $outdst     "
echo "outntuple    $outntuple  "
echo "scriptdir    $scriptdir  "
echo "tmpdir       $tmpdir     "

#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir


##pisaToDSTLinker with new library
/opt/phenix/core/bin/LuxorLinker.pl -1 $runnum

yes y | cp ${scriptdir}/embed.C .
yes y | cp ${scriptdir}/run_embed.C .
yes y | cp ${scriptdir}/Fun4All_embedeval.C .
yes y | cp ${scriptdir}/embed_IOManager.C    .
yes y | cp ${scriptdir}/svxPISA.par         .
yes y | cp ${scriptdir}/svx_threshold.dat   .

echo "pwd"
pwd
ls

echo ".x Fun4All_embedeval.C($evtnum, "'"'$inputreal'"'", "'"'$inputsim'"'", "'"'$outdst'"'", "'"'$outntuple'"'", $runnum);" >  cmd.input
#echo ".x run_embed.C($evtnum, "'"'$inputreal'"'", "'"'$inputsim'"'", "'"'$outdst'"'");" >  cmd.input
#echo ".x embed.C($evtnum, "'"'$inputreal'"'", "'"'$inputsim'"'", "'"'$outdst'"'", "'"'$outntuple'"'", $runnum);" >  cmd.input
echo ".q" >> cmd.input

##run root
root -b < cmd.input

#remove tmp dir
cd $HOME
#rm -fr $tmpdir
echo "removed $tmpdir"
