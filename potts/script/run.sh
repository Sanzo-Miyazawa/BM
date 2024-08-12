#! /bin/sh

echo "Please define PRGMDIR as the root of this program."
echo "Do $ cp -p $PRGMDIR/script/RunBM.scala  $OUTDIR/
echo "And then modify RunBM.scala if necessary."
exit

PRGMDIR=".../potts" 
PROTEIN="PF00018uniq"

OUTDIR=$PRGMDIR/out/$PROTEIN
# cp -p "$PRGMDIR/script/RunBM.scala"  $OUTDIR/
PRGM=$OUTDIR/RunBM.scala
LOG=$OUTDIR/RunBM.log

cd $PRGMDIR

JAR=$PRGMDIR/jarFiles/org-sanzo-potts_3-1.7.6.jar
FATJAR=$PRGMDIR/jarFiles/org-sanzo-potts-assembly-1.7.6.jar

#export SCALA="scala3"
export SCALA="scala3 -classpath $FATJAR"

PIPIJ="$PRGMDIR/data/PF00018uniq/PiaPijab.out.gz"
REP="$PRGMDIR/data/PF00018uniq/representativeMSA.fasta.gz"

$PRGM $OUTDIR $PIPIJ $REP >& $LOG < /dev/null &

