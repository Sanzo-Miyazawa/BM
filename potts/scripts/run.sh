#! /usr/bin/sh

echo "Please copy this script together with one of RunBM.scala, RunMC_A.scala or RunMC_B.scala into OUTDIR/...."
echo "And then modify the Run*.scala for your purpose."
echo "At last, modify this script, and run it in the OUTDIR/..., into which all output and necessary files are located."
#exit

## 
OUTDIR=`dirname $0`
LIBDIR=$OUTDIR      # directory which includes jar Files.

JAR=$LIBDIR/org-sanzo-potts_3-1.7.7.jar
FATJAR=$LIBDIR/org-sanzo-potts-assembly-1.7.7.jar

export SCALA="scala3"
export SCALA="scala3 -classpath $FATJAR"
## 

PIPIJ="$OUTDIR/PiaPijab.out.gz"
REP="$OUTDIR/representativeMSA.fasta.gz"

case $0 in
   *runBM.sh)

	PRGM=$OUTDIR/RunBM.scala
	LOG=$OUTDIR/RunBM.log

	echo $PRGM $OUTDIR $PIPIJ $REP
	$PRGM $OUTDIR $PIPIJ $REP 1> $LOG 2>&1  < /dev/null &
	;;
   *runMC.sh)
	INT="$OUTDIR/log_Interactions_hJ_6187.txt.gz"
	PRGM=$OUTDIR/RunMC*.scala
	LOG=`dirname $PRGM`/`basename $PRGM .scala`.log

	echo $PRGM $OUTDIR $PIPIJ $REP $INT
	$PRGM $OUTDIR $PIPIJ $REP $INT 1> $LOG 2>&1  < /dev/null &
	;;
    *)
	exit
	;;
esac

