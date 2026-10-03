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

# --Maxim Potekhin /mxp/--

#if (! $?EMBEDDING_HOME) then       
#echo "Environment not set, exiting..."
#exit -1   
#endif

setenv prompt 1
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
source $i
end
source $HOME/.login
unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv ROOTSYS

set shift = 2000

setenv DATADIR /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embed/work

source /opt/phenix/bin/phenix_setup.csh new
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embedreco/install/lib:$LD_LIBRARY_PATH
setenv LD_LIBRARY_PATH /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embed/svx_cent_ana/install/lib:$LD_LIBRARY_PATH

echo $LD_LIBRARY_PATH

setenv sourcedir /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/embed/svx_cent_ana/build
if( $1 == 0) then
cd $sourcedir
make -j4
make install
endif
##################################

if( $#argv != 7) then
echo "reco_rawhit.csh num"
echo "   num       = job number"
echo "   evtnum    = Nevent"
echo "   inputreal = inputfile"
echo "   inputsim  = inputfile"
echo "   outdst    = outdst"
echo "   outntuple = outdst"
echo "   outntana  = outdst"
exit -1
endif

echo "Setting internal variables..."

set runnum    = "409471"
set jobno     = $1
set evtnum    = $2
#set inputsim  = $3
set inputreal = $4
set outdst    = $5
set outntuple = $6

set INPUT = `expr $shift + $jobno`
echo $INPUT

set outntana  = "$INPUT"_$7

set DIR = `printf "%05d" $INPUT`
set inputsim = /phenix/plhf/mitran/Simul/Dileptons/output_single/dst/dst_out_single_$DIR.root
set inputoscar = /phenix/plhf/mitran/Simul/Dileptons/output_single/single/$DIR.oscar.particles.dat
set inputvtx = /gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/output/vertexes.txt

set scriptdir = "$DATADIR"

set outdstdir   = $DATADIR/output
set outntdir    = $DATADIR/output
set outntanadir = $DATADIR/output
set tmpdir      = "/home/tmp/${USER}_job_$INPUT"

set inreal = `basename $inputreal`
set insim  = `basename $inputsim`

set runnum = `echo $inreal | awk -F'-' '{printf "%d", $2}'`



echo "runnum       $runnum     "
echo "jobno        $jobno      "
echo "evtnum       $evtnum     "
echo "inputreal    $inputreal  "
echo "inputsim     $inputsim   "
echo "outdst       $outdst     "
echo "outntuple    $outntuple  "
echo "outntana     $outntana   "
echo "inputvtx     $inputvtx   "
echo "inputoscar   $inputoscar "
echo "scriptdir    $scriptdir  "
echo "outdstdir    $outdstdir  "
echo "outntdir     $outntdir   "
echo "outntanadir  $outntanadir"
echo "tmpdir       $tmpdir     "

# exit


#move to wrk directory
if( ! -d $tmpdir ) then
mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir


##pisaToDSTLinker with new library
/opt/phenix/core/bin/LuxorLinker.pl -1 $runnum

cp ${scriptdir}/Fun4All_embedeval.C .
cp ${scriptdir}/Fun4All_embedeval_svx.C .
cp ${scriptdir}/embed_IOManager.C    .
cp ${scriptdir}/svxPISA.par         .
cp ${scriptdir}/svx_threshold.dat   .

#copy input file
cp $inputreal .
cp $inputsim .

echo "pwd"
pwd
ls -ltr

echo "yolo"
echo ".x Fun4All_embedeval_svx.C($evtnum, "'"'$insim'"'", "'"'$inreal'"'", "'"'$outdst'"'", "'"'$outntuple'"'", "'"'$outntana'"'", "'"'$inputvtx'"'", "'"'$inputoscar'"'", $runnum);" >  cmd.input
echo ".q" >> cmd.input

##run root
root -b < cmd.input #>& $HOME/root.log
#ls -ltr
#move
#mv -f $outdst    $outdstdir
mv -f $outntuple $outntdir
mv -f $outntana  $outntanadir
mv -f my-$outntana $outntanadir

#remove tmp dir
cd $HOME
rm -fr $tmpdir
echo "removed $tmpdir"

