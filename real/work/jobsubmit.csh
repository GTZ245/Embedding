#!/usr/local/bin/tcsh -f

if ( $#argv < 2) then
  echo "Usage : jobsubmit.csh [script] [arg]"
  exit 0
endif

if( ! -f $1 ) then
  echo "No mergelist : $1"
  exit 0
endif

###################

source /opt/phenix/bin/phenix_setup.csh

set submitdir   = $PWD
set jobbasedir  = "$submitdir/jobs"
set script      = $1
set jobno       = $2
set nevent      = $3
set inputfile   = $4

echo $script
echo $jobno    
echo $nevent   
echo $inputfile
set arglist = ($jobno $nevent $inputfile)

##################
# make jobfile

set jobdir = "$jobbasedir/job_$jobno"
if( ! -d $jobdir ) then
  mkdir $jobdir
endif

set jobfile = "$jobdir/job_${jobno}.job"
set outfile = "$jobdir/job_${jobno}.out"
set errfile = "$jobdir/job_${jobno}.err"
set logfile = "$jobdir/job_${jobno}.log"

echo "Universe        = vanilla"                   >  $jobfile
echo "Executable      = $script"                   >> $jobfile
echo "Arguments       = "\"$arglist\"              >> $jobfile
echo "GetEnv          = True"                      >> $jobfile
echo "+Job_Type       = "\"cas\"                   >> $jobfile
#echo "+Job_Type       = "\"VtxProduction\"            >> $jobfile
echo "+Experiment     = "\"phenix\"                   >> $jobfile
echo "Initialdir      = $submitdir"                >> $jobfile
echo "Requirements    = TotalDisk > 20000000"      >> $jobfile
echo "Output          = $outfile"                  >> $jobfile
echo "Error           = $errfile"                  >> $jobfile
echo "Log             = $logfile"                  >> $jobfile
#echo "Priority        = 9950"                      >> $jobfile
echo "Queue"                                       >> $jobfile


#/opt/condor/bin/condor_submit $jobfile
condor_submit $jobfile

