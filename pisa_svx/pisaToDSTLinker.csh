#! /bin/csh -f

echo ""
echo "pisaToDSTLinker setup for run15 simulations"

# configuration files
ln -sf /afs/rhic.bnl.gov/phenix/software/simulation/run15/pisaToDSTLinker.csh .

# rich cabling
ln -sf /afs/rhic.bnl.gov/phenix/software/simulation/head/crk_cabling_vrdc.txt .

# Dch calibrations for MC  xxx
ln -sf /afs/rhic.bnl.gov/phenix/software/calibration/sim00/DchCalibration.Mc .
ln -sf /afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run15/DchEfficiency.Real .
ln -sf /afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/run15/DchAlwaysDead.dat .

# central arm integrated field map xxx
ln -sf /afs/rhic.bnl.gov/phenix/software/calibration/run2004/fieldIntegral++.dat.run04 fieldIntegral.dat

# 3D magnetic field map
ln -sf /afs/rhic.bnl.gov/phenix/PHENIX_LIB/simulation/Sim3D++.root .

# muon arm material map
ln -sf /afs/rhic.bnl.gov/phenix/software/calibration/run15/pisafile.dat.cZ .

# root macros
cp /afs/rhic.bnl.gov/phenix/software/simulation/run15/pisaToDST_IOManager.C .
cp /afs/rhic.bnl.gov/phenix/software/simulation/run15/Fun4Muons_Pisa.C .


# input file - we generate PISAEvent.root locally

echo ""
echo "setup script completed"
echo "Type root -b -q Fun4Muons_Pisa.C"
echo "for 1000 dimuons."
echo "You might need to change "
echo "  the pisaToDST macro; "
echo "  the link the magnetic field maps;"
echo "  and the input file."
echo ""
echo "Please report any broken link by email to phenix-off-l@lists.bnl.gov"
