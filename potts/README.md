
## A program written in Scala for a Boltzman machine.


#### To build a jar file

Jar files are available in the subdirectory named "jarFiles".
If you want to make them, pleas refer to the following.

In the following, it is assumed that the version of scala to build is 3.3.3. 
If not, change "ThisBuild / scalaVersion" in build.sbt.

>bash$ PRGMDIR=          # the root directory of this program

>bash$ cd $PRGMDIR

>bash$ sbt

>sbt> compile

>sbt> packageBin

>sbt> packageDoc

>sbt> assembly           # to make a fat jar

>sbt> exit

>bash$ JAR=$PRGMDIR/./target/scala-3.2.0/org-sanzo-potts_3-0.3.0.jar

>bash$ FATJAR=$PRGMDIR/./target/scala-3.2.0/org-sanzo-potts-assembly-0.3.0.jar

Two scripts, RunBM.scala, RunMC_A.scala and RunMC_B.scala, in the $PRGMDIR/script are provided
as examples for using this program.  RunBM.scala includes parameter definitions
that are appropriate to PF00018.  Thus if it is used for other proteins,
they must be changed.

Please read the scripts to understand how to use this program.

#### To run the scala scripts on your system, see the following or run.sh that is found in the $PRGMDIR/script.

In the following it is assumed that the **scala version installed in your system is equal to the version 3 of
scala with which the jar file was built.**

>bash$ export SCALA="scala3 -classpath $FATJAR"	# please replace "scala3" with a command to run scala3.

>bash$ OUTDIR=		# the directory into which output files are created.

>bash$ PIPIJ="$PRGMDIR/data/PF00018uniq/PiaPijab.out.gz"

>bash$ REP="$PRGMDIR/data/PF00018uniq/representativeMSA.fasta.gz"

>bash$ $PRGMDIR/script/RunBM.scala $OUTDIR $PIPIJ $REP >& $OUTDIR/RunBM.log < /dev/null &

#### Reference: 

 1. arXiv:1909.05006 [q-bio.BM], 2019

 2. IEEE/ACM Transactions on Computational Biology and Bioinformatics, 2020 
    DOI: 10.1109/TCBB.2020.2993232

 3. https://gitlab.com/sanzo.miyazawa/BM/

#### For any question, send a mail to sanzo.miyazawa@gmail.com.    

2024-08-01 Sanzo Miyazawa



