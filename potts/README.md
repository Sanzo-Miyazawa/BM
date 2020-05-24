
## A program written in Scala for a Boltzman machine.


#### To build a jar file

In the following, it is assumed that the version of scala to build is 2.13. 
If not, change "ThisBuild / scalaVersion" in build.sbt.

>bash$ PRGMDIR=          # the root directory of this program

>bash$ cd $PRGMDIR

>bash$ sbt

>sbt> compile

>sbt> packageBin

>sbt> packageDoc

>sbt> assembly           # to make a fat jar

>sbt> exit

>bash$ JAR=$PRGMDIR/./target/scala-2.13/org-sanzo-potts_2.13-0.1.0.jar

>bash$ FATJAR=$PRGMDIR/./target/scala-2.13/org-sanzo-potts-assembly-0.1.0.jar

Two scripts, RunBM.scala and RunMC.scala, in the $PRGMDIR/scripts are provided as examples for using this program.
RunBM.scala includes parameter definitions that are appropriate to PF00595.
Thus if it is used for other proteins, they must be changed.
Please read the scripts to understand how to use this program.

#### To run a scala script on your system.

In the following it is assumed that the **scala version installed in your system is equal to the version of
scala with which the jar file was built.**

>bash$ export SCALA="scala -nc -cp $FATJAR"

>bash$ OUTDIR=...        # the directory into which output files are created.

>bash$ PIPIJ="$PRGMDIR/data/PF00595uniqSeq/PiaPijab.out.gz"

>bash$ $PRGMDIR/script/RunBM.scala $OUTDIR $PIPIJ >& $OUTDIR/RunBM.log < /dev/null &

#### Reference: 

 1. 

 2.

#### For any question, send a mail to sanzo.miyazawa@gmail.com.    

2020-04-06 Sanzo Miyazawa



