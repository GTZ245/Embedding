#!/usr/local/bin/tcsh -f

setenv HOME /phenix/u/mitran
setenv prompt 1
source /etc/csh.login
foreach i (/etc/profile.d/*.csh)
  source $i
end
source $HOME/.login
unsetenv OFFLINE_MAIN
unsetenv ONLINE_MAIN
unsetenv ROOTSYS

source /opt/phenix/core/bin/phenix_setup.csh new
setenv LD_LIBRARY_PATH /phenix/plhf/mitran/Simul/Dileptons/real/dstmerge/install/lib:/phenix/plhf/mitran/Simul/Dileptons/embedreco/install/lib/:$LD_LIBRARY_PATH
setenv TSEARCHPATH /phenix/hhj/hachiya/15.08/embed/embed

##################################

if( $#argv != 3) then
  echo "reco_rawhit.csh num"
  echo "   num = job number"
  echo "   evtnum = Nevent"
  echo "   infile = inputfile"
  exit -1
endif

set sourcedir = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/dstmerge/build/"
cd $sourcedir
make -j4
make install

set jobno     = $1
set evtnum    = $2
set inputfile = $3

set scriptdir = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work"
set outputdir = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/work/outputfull"
set tmpdir    = "/gpfs/mnt/gpfs02/phenix/plhf/plhf1/mitran/Simul/Dileptons/real/output/000_$jobno"

set infile = `basename $inputfile`
set runnum = `echo $infile | awk -F'-' '{printf "%d", $2}'`
set seqnum = `basename $inputfile .PRDFF| awk -F'-' '{printf "%s", $3}'`

echo all parameters:
echo $jobno
echo $evtnum
echo $inputfile
echo $runnum
echo $seqnum
echo $scriptdir
echo $tmpdir


#move to wrk directory
if( ! -d $tmpdir ) then
  mkdir -p $tmpdir
endif
echo "cd $tmpdir"
cd       $tmpdir


##pisaToDSTLinker with new library
/afs/rhic/phenix/software/calibration/data/LuxorLinker.pl -1 $runnum
/opt/phenix/core/bin/LuxorLinker.pl -1 $runnum

cp ${scriptdir}/Fun4All_CA_merge.C .
cp ${scriptdir}/OutputManager.C    .
cp ${scriptdir}/QA.C               .
cp ${scriptdir}/TrigSelect.C.run14auau200  TrigSelect.C
cp ${scriptdir}/rawdatacheck.C     .

#copy input file
#cp $inputfile .
#${HOME}/copy_prdf.pl $infile


#generate input file
echo ".x Fun4All_CA_merge.C($evtnum, "'"'$inputfile'"'");" >  cmd.input
echo ".q"                                               >> cmd.input


##run root
root -b < cmd.input

#move
ls
echo "mv -f CNTmerge_MB_*.root $outputdir"
mv -f CNTmerge_MB*.root $outputdir
mv vertexes.txt $outputdir/vertexes_MB-0000$runnum-$seqnum.txt

#remove tmp dir
cd $HOME
rm -fr $tmpdir
echo "removed $tmpdir"

