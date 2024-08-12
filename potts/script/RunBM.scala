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

object RunBM {
  import scala.collection.immutable.Seq
  import scala.collection.immutable.HashMap

  import org.biojava.nbio.core.util.InputStreamProvider

  import org.sanzo.potts.MCMC
  import org.sanzo.potts.BM

  import breeze.linalg.DenseVector
  import breeze.linalg.DenseMatrix

  def processFiles(args: Array[String], 
        outputFiles: HashMap[String, java.io.File]) = {

    val isp = new InputStreamProvider()
    val file = new java.io.File( args(0) )
    val fin = isp.getInputStream(file)

    val bufSource = scala.io.Source.fromInputStream(fin)
    val lines = bufSource.getLines()

    val (fia, fijab, stateOrderString) = MCMC.readPiPij(lines)

    val initialInteractions: Option[MCMC.Interactions] = None           //null
        
  /*
    import org.sanzo.biojava.SequenceIO
    import org.sanzo.sequence.CasePreservingProteinSequenceArray
    import org.sanzo.sequence.CompoundOrderMap
    import org.sanzo.sequence.SequenceArray.pairIndex

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
    val nMC = neff      // 5000.0

    val maxNSamples = ( nMC / nIndependentMC ).ceil.toInt
    val (initialNSamples: Int, incrementSamples: Double ) = {
        val init = (5000.0 / nIndependentMC).ceil.toInt
        val incr = (maxNSamples - init) / 1000.0
        if ( incr >= 0.0 ) ( init, incr ) else ( maxNSamples, 0.0 )
    }

    val lambdaPhi = 0.001 / neff
    val lambdaPhij = math.pow(10.0, -1.0)

    val minLearningStepsForBestKL = 1000
    val minLearningSteps = 1200
    val maxNoLearnings = 100

    val mBoltzman = new BM(outputFiles, stateOrderString,
                        fia, fijab, 
                        regTerm = "GL1L2", propL1h = 0.0, propL1J = 1.0,
                        //regTerm = "L1L2", propL1h = 0.0, propL1J = 0.9,
                        //regTerm = "L1L2", propL1h = 0.9, propL1J = 0.9,
                        //regTerm = "L2", propL1h = 0.0, propL1J = 0.0,
                        lambdaPhi = lambdaPhi, lambdaPhij = lambdaPhij,
                        //betaV = 0.9, betaM = 0.0, eps = 0.0 )         // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.95, betaM = 0.0, eps = 0.0 )        // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.99, betaM = 0.0, eps = 0.0 )        // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-8 )    // for ModRAdam (= ModRAdamMax) 
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-8 )    // for RAdam
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-7 )    // for ModAdamSum
                        betaV = 0.999, betaM = 0.9, eps = 1.0e-8 )    // for ModAdam (= ModAdamMax) 
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-8 )    // for Adam
                        //betaV = 0.5, betaM = 1.2, eps = 1.0e-8 )      // for RPROP-LR; betaV=rateDecrease , betaM=rateIncrease

     val results = mBoltzman.runBM(
                optionInitialInteractions = initialInteractions,
                initialConfigurations = (new Array[Array[Byte]](nIndependentMC)),
                nInitialIterationsPerUnit = 100,
                everyNIterationsPerUnit = 10,

                nSamples = 0,                   // if nSamples == 0, initialNSamples, incrementSamples, and maxNSamples are used.
                initialNSamples = initialNSamples,
                incrementSamples = incrementSamples,
                maxNSamples = maxNSamples,
                mcmcKernel = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

                initialT = 1.2,
                annealingRate = 0.99,
                maxExtendedIterations = 10,

                gauge = "phi_zeroSum",  // "phi_zeroSum" for L2 and L2GL1, "ungauged" for all

                gradientDescentMethod = "ModAdam",

              //learningRate = 1.0E-1,          // for NAG
              //learningRate = 1.0E-2,          // for RPROP-LR
              //        minLearningRate = 1.0E-5,       // for RPROP-LR
              //        maxLearningRate = 10.0,         // for RPROP-LR
              //        rateDecrease = 0.5,             // for RPROP-LR
              //        rateIncrease = 1.2,             // for RPROP-LR
              //learningRate = 1.0E-3,          // for Adam
              //learningRate = 1.0E-3,          // for RAdam
              //learningRate = 1.0E-1,                  // for ModAdamSum
                learningRate = 1.0E-2, //5.0E-2,        // for ModAdam (= ModAdamMax)
              //learningRate = 1.0E-2,                  // for ModRAdamMax

                epsForKL = 1.0e-5,
                nBestKLs = 1,

                minLearningStepsForBestKL = minLearningStepsForBestKL,
                minLearningSteps = minLearningSteps,
                maxNoLearnings = maxNoLearnings,
                logInterval = 10
                )
    
  }

  def main(argsP: Array[String]): Unit = {
    //val outfile = new File("/tmp/output.txt")
    //val fout = new FileOutputStream(outfile)

      scala.sys.process.stdout.print("#")
      argsP.foreach( i => scala.sys.process.stdout.print(" " + i) )
      scala.sys.process.stdout.print("\n")

      val args = argsP.drop(1)
      val ioFiles = HashMap( 
                ("outputDir", new java.io.File(args(0) ) ),
                ("jobControlParamFile", new java.io.File(new java.io.File(args(0)),"jobControlParam.dat" ) )
                 )
      args(1) match {
        case _ => {
                //val fout = scala.sys.process.stdout   //System.out 
                processFiles(args.drop(1), ioFiles)
                }
      }
   }

}

