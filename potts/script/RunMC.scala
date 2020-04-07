#! /bin/sh
# "$0" is passed as the first argument.
if [ "$SCALA" = "" ] ; then
        exec scala -savecompiled "$0" "$0" "$@"
else
        exec $SCALA -savecompiled "$0" "$0" "$@"
fi

!#

import org.sanzo.potts.MCMC
import org.sanzo.potts.MCMC.printPiPij

object RunMC {
  import scala.collection.immutable.Seq

  import org.biojava.nbio.core.util.InputStreamProvider
  import org.sanzo.sequence.CasePreservingProteinSequenceArray
  import org.sanzo.biojava.SequenceIO
  import org.sanzo.sequence.CompoundOrderMap
  import org.sanzo.sequence.SequenceArray.pairIndex

  import org.sanzo.potts.MCMC
  import org.sanzo.potts.BM

  import breeze.linalg.DenseVector
  import breeze.linalg.DenseMatrix

  def processFiles(args: Array[String], 
	outputFiles: collection.mutable.HashMap[String, java.io.File]) = {

    val isp = new InputStreamProvider()
    val file = new java.io.File( args(0) )
    val fin = isp.getInputStream(file)

    val bufSource = scala.io.Source.fromInputStream(fin)
    val lines = bufSource.getLines()

    val (fia, fijab, stateOrderString) = MCMC.readPiPij(lines)

    val file2 = new java.io.File( args(1) )
    val fin2 = isp.getInputStream(file2)

    val bufSource2 = scala.io.Source.fromInputStream(fin2)
    val lines2 = bufSource2.getLines()

    val (hia, jijab, stateOrderString2) = MCMC.readhJ(lines2)
    require( stateOrderString == stateOrderString2)

    val initialInteractions = MCMC.Interactions(hia, jijab)
	
  /*
    val (fia, fijab, stateOrderString) =
    {
		import breeze.linalg.DenseVector
		import breeze.linalg.DenseMatrix
                val nUnits = 10
                val nStatsOfUnit = 21
		val fia = (new Array[DenseVector[Double]](nUnits)).map(
			pi => DenseVector.fill(nStatsOfUnit){1.0/nStatsOfUnit} )
                val fijab = new Array[DenseMatrix[Double]]( (nUnits * (nUnits - 1)) / 2)
                (0 until fia.size).foreach ( i => {
                  (0 until i).foreach ( j => {
                        val (ij, a, b) = pairIndex(i, j, 0, 1)
                        if ( a == 0 )
                                fijab(ij) = fia(i) * fia(j).t
                        else
                                fijab(ij) = fia(j) * fia(i).t
                  } )
                } )
    		val stateOrderString = "ARNDCQEGHILKMFPSTWYV-"
		(fia, fijab, stateOrderString)
    }
   */

    val nUnits = fia.size
    val nStatesOfUnit = fia(0).size
    val nIndependentMC = 20
    //val neff = 19473.9 // PF00153
    val neff = 4748.76   // PF00595uniqSeq

    val nSamples = ( neff / nIndependentMC ).ceil.toInt

    val mcmc = new MCMC( outputFiles, stateOrderString,
		initialInteractions, proposedPia = fia ) 

        val (initialStates, nExtendedIterations, independentSamplings) =
          mcmc.runMC(
		initialConfigurations = (new Array[Array[Byte]](nIndependentMC)),
		nInitialIterationsPerUnit = 100,
		everyNIterationsPerUnit = 10,
		nSamples = nSamples,
		//initialT = 1.2,
		initialT = 2.0,
		annealingRate = 0.99,
                maxExtendedIterations = 10,
		kernel = mcmc.kernelMultiBlockMH
		)

 	val fileIndpendentMCsamplings = new java.io.File(outputFiles("outputDir"), s"MC_samples.fasta")
	val outSamples =  new java.io.PrintStream(fileIndpendentMCsamplings)
        MCMC.printIndependentMCsamplings(outSamples, stateOrderString,
                independentSamplings )
    
  }

  def main(argsP: Array[String]): Unit = {
    //val outfile = new File("/tmp/output.txt")
    //val fout = new FileOutputStream(outfile)

      val args = argsP.drop(1)
      val outputFiles = collection.mutable.HashMap( 
		("outputDir", new java.io.File(args(0) ) ) )
      args(1) match {
	case _ => {
    		//val fout = scala.sys.process.stdout	//System.out 
    		processFiles(args.drop(1), outputFiles)
		}
      }
   }

}

