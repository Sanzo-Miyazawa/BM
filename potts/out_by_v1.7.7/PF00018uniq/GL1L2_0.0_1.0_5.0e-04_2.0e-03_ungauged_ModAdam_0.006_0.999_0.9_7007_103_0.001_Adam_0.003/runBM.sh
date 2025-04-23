#! /usr/bin/sh

echo "Please copy this script together with one of RunBM.scala, RunMC_A.scala or RunMC_B.scala into OUTDIR/...."
echo "And then modify the Run*.scala for your purpose."
echo "At last, modify this script, and run it in the OUTDIR/..., into which all output and necessary files are located."
#exit

## 

OUTDIR=`dirname $0`
LIBDIR=$OUTDIR      # directory which includes jar Files.

PIPIJ="$OUTDIR/PiaPijab.out.gz"
REP="$OUTDIR/representativeMSA.fasta.gz"

PRGM=$OUTDIR/RunBM.scala
LOG=$OUTDIR/RunBM.log

## 

JAR=$LIBDIR/org-sanzo-potts_3-1.7.7.jar
FATJAR=$LIBDIR/org-sanzo-potts-assembly-1.7.7.jar

export SCALA="scala3 -classpath $FATJAR"
export SCALA="scala3"

$PRGM $OUTDIR $PIPIJ $REP 1> $LOG 2>&1  < /dev/null &

