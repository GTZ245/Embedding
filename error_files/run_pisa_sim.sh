#!/bin/sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/mnt/gpfs02/phenix/plhf/plhf3/nnovitzk/mazsi_Test/ccnt/source/emc-evaluation/build/.libs/:/direct/phenix+u/tongzhouguo/install
export INSTALL=/direct/phenix+u/tongzhouguo/install
export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH
export TSEARCHPATH=/direct/phenix+u/tongzhouguo/install
export DCACHE_DOOR=phnxcore03.rcf.bnl.gov:1094
export GSEARCHPATH ${GSEARCHPATH}:DCACHE

echo $LD_LIBRARY_PATH

chmod g+rx ${_CONDOR_SCRATCH_DIR}
cd ${_CONDOR_SCRATCH_DIR}
#cd /gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/error_files/
echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="

INPUT=$(( $1 + $2 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`
mkdir -p $DIR
pushd $DIR
echo start the HELIOS simulation with PISA in directory $DIR

root -l -b << EOF
gSystem->Load("/phenix/plhf/roli/HELIOS/simulation/WriteEvent_C.so");
.x /gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/error_files/pack/WriteROOT2Oscar10K.C("/phenix/plhf/roli/HELIOS/simulation/singlePi0.root","oscar.particles.dat",$INPUT,$3)
EOF

NEVT=1000

echo "==============================================="
echo "============= START PISA NOW =================="
echo "==============================================="
mkdir pisa
pushd pisa

cp ../oscar.particles.dat /phenix/plhf/tongzhouguo/error_files/pack/oscars/$DIR.oscar.particles.dat
mv ../oscar.particles.dat .
#cp /phenix/plhf/tongzhouguo/error_files/pisa_wo_vtx/* .
cp /phenix/plhf/tongzhouguo/pisa_svx/* .

sed -i '4s/.*/ptrig '$NEVT'/' glogon.kumac

pisa<pisa.input

echo "==============================================="
echo "================ PISA TO DST =================="
echo "==============================================="

#root -b -q -l pisaToDST.C
root -b -q -l pisaToDST_VTX.C




cp dst_out.root /phenix/plhf/tongzhouguo/error_files/pack/dsts/dst_out_pi0_$DIR.root
mv dst_out.root ../dst_out_pi0_$DIR.root
popd  

#SimFileName=`echo dst_out_pi0_$DIR.root`
#ln -s /phenix/plhf/tongzhouguo/error_files/run_TTreemaker.C .  
#output=/phenix/plhf/tongzhouguo/error_files/trees/sim_trees_$DIR.root
#echo $SimFileName
#root -l -b -q 'run_TTreemaker.C("'$SimFileName'", "'$output'", 0, 1)'

rm -rf pisa
popd
rm -r $DIR
