#!/bin/sh


export INSTALL=/direct/phenix+u/tongzhouguo/install
export LD_LIBRARY_PATH=$INSTALL/lib:$LD_LIBRARY_PATH

echo "==============================================="
echo "============= START SIMULATION  ==============="
echo "==============================================="

INPUT=$(( $1 + $2 ))
echo $INPUT

DIR=`printf "%05d" $INPUT`

SYSTEM=0;

pushd output
ln -s /gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/yuri_embed/embed/analysis/emb/lookup_3D_one_phi.root .
ln -s $PWD/../get_Ntag_uncorr_C.so .
file=/gpfs/mnt/gpfs02/phenix/plhf/plhf1/tongzhouguo/test_tree5/trees/sim_trees_$DIR.root
output=$PWD/output_$DIR.root

root -l -b << EOF
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisEvent");
gSystem->Load("/direct/phenix+u/tongzhouguo/install/lib/libDileptonAnalysisReco");
gSystem->Load("/direct/phenix+u/roli/scratch/HELIOS/WriteEvent_C.so");
gSystem->Load("get_Ntag_uncorr_C.so");
get_Ntag_uncorr("$file","$output",$SYSTEM,$INPUT,$3)
EOF

