#! /bin/sh
# "$0" is passed as the first argument.
if [ "$SCALA" = "" ] ; then
        exec scala -savecompiled "$0" "$0" "$@"
else
        exec $SCALA -savecompiled "$0" "$0" "$@"
fi

!#


object RunBM {
  import java.io.PrintStream

  import scala.collection.immutable.Seq
  import scala.collection.immutable.ArraySeq
  import scala.collection.immutable.Vector
  import scala.collection.immutable.HashMap

  import org.biojava.nbio.core.util.InputStreamProvider

  import org.sanzo.potts.MCMC
  import org.sanzo.potts.MCMC.readFasta
  import org.sanzo.potts.MCMC.printPiPij
  import org.sanzo.potts.BM
  import org.sanzo.potts.LearningRate

  import breeze.linalg.DenseVector
  import breeze.linalg.DenseMatrix

  def processFiles(args: Array[String], 
        ioFiles: HashMap[String, java.io.File]) = {

    val isp = new InputStreamProvider()

    val (effectiveNumberOfSamples, fia, fijab, stateOrderString) =
      if ( args.size < 1 ) {
        require( args.size > 0)
        (0.0, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]](), "")
      } else {
        val file = new java.io.File( args(0) )
        val fin = isp.getInputStream(file)

        val bufSource = scala.io.Source.fromInputStream(fin)
        val lines = bufSource.getLines()

        val (effectiveNumberOfSamples, fia, fijab, stateOrderString) = MCMC.readPiPij(lines)
        bufSource.close()
        (effectiveNumberOfSamples, fia, fijab, stateOrderString)
      }

  //val miniBatch = false   // for the v0.3 type of paramaters
    val miniBatch = true    // =true is better for minibatch of  >v1.0

    val (initialConfIDs, initialConfigurations: ArraySeq[IArray[Byte]], nSamples, maxNSamples, initialNSamples, incrementSamples: Double,
         nInitialIterationsPerUnit, everyNIterationsPerUnit ) =
      if ( ! miniBatch ) { 
        // obsolete: for the v0.3 type of paramaters

        val nIndependentMC = 20
        val nMC = effectiveNumberOfSamples      // 5000.0

        val nSamples = 0   // for the v0.3 type of paramaters
        val maxNSamples = ( nMC / nIndependentMC ).ceil.toInt
        val (initialNSamples: Int, incrementSamples: Double ) = {
          val init = (5000.0 / nIndependentMC).ceil.toInt
          val incr = (maxNSamples - init) / 1000.0
          if ( incr >= 0.0 )
            ( init, incr ) 
          else
            ( maxNSamples, 0.0 )
        }
        val nInitialIterationsPerUnit = 100
        val everyNIterationsPerUnit = 10

        (ArraySeq[String](), ArraySeq.fill(nIndependentMC)(IArray[Byte]()), nSamples,
         maxNSamples, initialNSamples, incrementSamples, nInitialIterationsPerUnit, everyNIterationsPerUnit)
      } else {
        val fastaFile = new java.io.File( args(1) )
        val fastaIn = isp.getInputStream(fastaFile)
        val fastaLines = scala.io.Source.fromInputStream(fastaIn).getLines()

        val (idArraySeq, seqArraySeq, numArraySeq) =
            readFasta(fastaLines, stateOrderString = stateOrderString, caseSensitive = false)

        val initialConfigurations = numArraySeq

        val initialConfIDs = idArraySeq.map( id => id.split("[ \t]")(0) )

        fastaIn.close()

        val nSamples = 1   // the number of samples per MC for minibatch of the v1.0
        val (initialNSamples, incrementSamples, maxNSamples) = (nSamples, 0.0, nSamples)

        val nInitialIterationsPerUnit = 0
        val everyNIterationsPerUnit = 10
        
        (initialConfIDs, initialConfigurations , nSamples, maxNSamples, initialNSamples, incrementSamples, nInitialIterationsPerUnit, everyNIterationsPerUnit)
     }

    val initialInteractions: Option[MCMC.Interactions] =
      if ( args.size < 3 ) {
        None           //null for the v0.3 type of paramaters
      } else {
        val file = new java.io.File( args(2) )
        val fin = isp.getInputStream(file)
        val bufSource = scala.io.Source.fromInputStream(fin)
        val lines = bufSource.getLines()

      //val (hia, jijab, orderString) = MCMC.readhJ(lines)
        val (interactions, orderString) = MCMC.readInteractions(lines)
        bufSource.close()
        assert(stateOrderString == orderString)
      //Option(MCMC.Interactions(hia, jijab))
        Option(interactions)
      }

    // for PF00018uniq
    val lambdaPhij = 2.0e-03
    val lambdaPhi  = 4.0e-04

    val (nEff, nRep) = (effectiveNumberOfSamples, initialConfigurations.size)

    val fullBatchSize = initialConfigurations.size
    val miniBatchSize = { import scala.util.Sorting

                          val minSize = 100
                          if ( fullBatchSize <= minSize ) {
                              fullBatchSize
                          } else {
                              val r = Vector.range(minSize, (minSize * 1.1).toInt )
                              val t = r.map{i =>  (i, (fullBatchsize % i) / i.toDouble) }
                              t.fold(t(0)){ (x, y) => if (x._2 <= y._2) x else y }._1
                          }
                        }

    val miniBSize = BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
    val nMBsInFB = BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )
    val stepsPerEpoch = nMBsInFB

    val mBoltzman = new BM(ioFiles, stateOrderString,
                        effectiveNConfigs = nEff, fia, fijab, 
                        regTerm = "GL1L2", propL1h = 0.0, propL1J = 1.0 )

                        //regTerm = "L1L2", propL1h = 0.0, propL1J = 0.9,
                        //regTerm = "L1L2", propL1h = 0.9, propL1J = 0.9,
                        //regTerm = "L2", propL1h = 0.0, propL1J = 0.0,
                        //lambdaPhi = lambdaPhi, lambdaPhij = lambdaPhij,
                        //betaV = 0.9, betaM = 0.0, eps = 0.0 )         // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.95, betaM = 0.0, eps = 0.0 )        // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.99, betaM = 0.0, eps = 0.0 )        // for NAG:  learningRate / (1 - betaV) should be constant.
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-6 )    // for ModRAdam (= ModRAdamMax) 
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-6 )    // for RAdam
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-6 )    // for ModAdamSum
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-6 )      // for ModAdam (= ModAdamMax) 
                        //betaV = 0.999, betaM = 0.9, eps = 1.0e-6 )    // for Adam
                        //betaV = 0.9, betaM = 0.9, eps = 1.0e-5 )      // for Adadelta
                        //betaV = 0.5, betaM = 1.2, eps = 1.0e-6 )      // for RPROP-LR; betaV=rateDecrease , betaM=rateIncrease

    val gradientDescentMethod = "ModAdam"  // + Adam by jobControlParam.dat    // NAG, Adadelta, RPROP-LR, RAdam, Adam, ModRAdam, ModAdam

              //maxLR = 1.0E-1          // for NAG
              //maxLR = 1.0E-2          // for RPROP-LR
              //        minLearningRate = 1.0E-5       // for RPROP-LR
              //        maxLearningRate = 10.0         // for RPROP-LR
              //        rateDecrease = 0.5             // for RPROP-LR
              //        rateIncrease = 1.2             // for RPROP-LR
              //maxLR = 1.0E-3          // for Adam
              //maxLR = 1.0E-3          // for RAdam
              //maxLR = 1.0E-1                  // for ModAdamSum
              //maxLR = 1.0E-2 //5.0E-2         // for ModAdam (= ModAdamMax)
              //maxLR = 1.0E-2                  // for ModRAdamMax

    val maxLR = 0.6E-02    // 1.0E-2  // 5.0E-2,        // for ModAdam (= ModAdamMax)

    val nEpochs = (1800.0/stepsPerEpoch).ceil.toInt    //(2200.0/stepsPerEpoch).ceil.toInt
    val minLearningSteps = nEpochs * stepsPerEpoch
    val minLearningStepsForBestKL = (minLearningSteps * 0.95).toInt
    val maxNoLearnings =  (minLearningSteps * 0.05).toInt

    val learningRate = new LearningRate(warmupEpochs = 100.0,                //minLearningSteps * 0.05,
                     maxLR = maxLR, maxLREpochs = minLearningSteps - 100.0,  //minLearningSteps * 0.95, //scala.Int.MaxValue.toDouble , // minLearningSteps.toDouble,
                     coolingMethod = "(1+at)^b", constPerEpoch = 0.1, power = -0.5,
                   //coolingMethod = "1/(1+at)", epochsUntilHalf = 10.0,
                   //coolingMethod = "rate", coolRate = 0.99, minLR = maxLR * 0.01,
                     stepsPerEpoch = 1.0 )

  // The second SGD method following the first SGD method defined in the RunBM(...).
    val minLearningStepsFor2thSGD = minLearningSteps + ( (1000.0/stepsPerEpoch).ceil.toInt ) * stepsPerEpoch
    val outJobControlParamFile = new PrintStream(ioFiles("jobControlParamFile") )
    outJobControlParamFile.print(s"minLearningSteps\t${minLearningStepsFor2thSGD}\n")
    outJobControlParamFile.print("gradientDescentMethod\tAdam\n")
    outJobControlParamFile.print("betaV\t0.999\n")
    outJobControlParamFile.print("betaM\t0.9\n")
    outJobControlParamFile.print("maxLR\t0.3E-2\n")
    outJobControlParamFile.close()

    val logInterval = ( (math.max(100.0, stepsPerEpoch) / stepsPerEpoch + 0.5).floor * stepsPerEpoch).toInt

    val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) * nMBsInFB  + 0.5 ).floor.toInt
  //val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) + 0.5 ).floor.toInt * nMBsInFB

    val results = mBoltzman.runBM(
                optionInitialInteractions = initialInteractions,
                initialConfigurations = initialConfigurations,

                sampleIDs = initialConfIDs,

                sigma_initial_J = 1.0e-3,             // Jijab is initialized by the Gaussian of this standard deviation and mean 0.0 for optionInitialInteractions == None/null

                nInitialIterationsPerUnit = nInitialIterationsPerUnit,
                everyNIterationsPerUnit = everyNIterationsPerUnit,

                miniBatchSize = miniBSize,             // if miniBatchSize > 0, then nSamples = 1

                nSamples = nSamples,                   // if nSamples == 0, initialNSamples, incrementSamples, and maxNSamples are used.

                initialT = 1.0,         // 1.2
                annealingRate = 0.99,
              //maxExtendedIterations = 1,
                mcmcKernel = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

                lambdaPhi = lambdaPhi, lambdaPhij = lambdaPhij,

              //gauge = "phi_zeroSum",  // "phi_zeroSum" for L2 and GL1L2, "ungauged" for all
                gauge = "ungauged",  // "phi_zeroSum" for L2 and L2GL1, "ungauged" for all
                gradientDescentMethod = gradientDescentMethod,

                betaV = 0.999, betaM = 0.9, eps = 1.0e-6,    // for ModAdam (= ModAdamMax)

                learningRate = learningRate,

                pseudoNCounts = 10.0,

                nBestKLs = 1,

                minLearningStepsForBestKL = minLearningStepsForBestKL,
                minLearningSteps = minLearningSteps,
                maxNoLearnings = maxNoLearnings,

                nMiniBatchesForEnsembleAve = nMiniBatchesForEnsembleAve,
                logInterval = logInterval
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

