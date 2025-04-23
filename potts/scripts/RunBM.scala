#! /bin/sh
#! LANG=${LANG:-en_US.UTF-8} /usr/local/bin/scala-cli shebang
# "$0" is passed as the first argument.
SCALA=scala3
if [ "$SCALA" = "" ] ; then
        #exec scala -savecompiled "$0" "$0" "$@"        # for scala <= 3.4
        #exec scala-cli "$0" -- "$0" "$@"               #for .sc with CLASSPATH
        #exec /usr/local/bin/scala-cli "$0" -- "$0" "$@" #for .sc
        #exec scala-cli "$0" -M RunBM -- "$0" "$@"      #for .scala
         LANG=${LANG:-en_US.UTF-8} exec scala "$0" -M RunBM -- "$0" "$@"        #for .scala
else
        #exec $SCALA -savecompiled "$0" "$0" "$@"       # for scala <= 3.4
         LANG=${LANG:-en_US.UTF-8} exec $SCALA "$0" -M RunBM -- "$0" "$@" #for .scala
fi

!#

object RunBM {
  import java.io.File
  import java.io.PrintStream

  import scala.collection.immutable.Seq
  import scala.collection.immutable.ArraySeq
  import scala.collection.immutable.Vector
  import scala.collection.immutable.HashMap

//import org.biojava.nbio.core.util.InputStreamProvider

  import org.sanzo.potts.MCMC
  import org.sanzo.potts.MCMC.readFasta
  import org.sanzo.potts.MCMC.printPiPij
  import org.sanzo.potts.BM
  import org.sanzo.potts.LearningRate
  import org.sanzo.potts.zip.gunzipInputStream

  import breeze.linalg.DenseVector
  import breeze.linalg.DenseMatrix

  def processFiles(args: Array[String], 
        ioFiles: HashMap[String, java.io.File]) = {

// Use scala.util.Random.setSeed(seed: Long), if you want consistent behavior with the same random numbers from execution to execution.
    val seedForRandom = 19L
    scala.util.Random.setSeed(seedForRandom)
    scala.sys.process.stdout.print(s"# seed for scala.util.Random = ${seedForRandom}\n")
// Use org.sanzo.potts.SetSeedForRand.randBasis instead of breeze.stats.distributions.Rand.Fixed.randBasis
    org.sanzo.potts.SetSeedForRand(0)
    scala.sys.process.stdout.print(s"# seed for org.sanzo.potts.SetSeedForRand.randBasis = ${org.sanzo.potts.SetSeedForRand.seed}\n")

    val (effectiveNumberOfSamples, fia, fijab, stateOrderString) =
      if ( args.size < 1 ) {
        require( args.size > 0)
        (0.0, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]](), "")
      } else {
        val file = new java.io.File( args(0) )
        val fin = gunzipInputStream(file)

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
        val fastaIn = gunzipInputStream(fastaFile)
        val fastaLines = scala.io.Source.fromInputStream(fastaIn).getLines()

        val (idArraySeq, seqArraySeq, numArraySeq) =
            readFasta(fastaLines, stateOrderString = stateOrderString, caseSensitive = false)

        fastaIn.close()

        assert(idArraySeq.size == seqArraySeq.size && idArraySeq.size == numArraySeq.size )

        val initialConfigurations = numArraySeq

        val initialConfIDs = idArraySeq.map( id => id.split("[ \t]")(0) )

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
        val fin = gunzipInputStream(file)
        val bufSource = scala.io.Source.fromInputStream(fin)
        val lines = bufSource.getLines()

      //val (hia, jijab, orderString) = MCMC.readhJ(lines)
        val (interactions, orderString) = MCMC.readInteractions(lines)
        bufSource.close()
        assert(stateOrderString == orderString)
      //Option(MCMC.Interactions(hia, jijab))
        Option(interactions)
      }

    val protein = "PF00018uniq"

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
              //maxLR = 1.0E-2 //0.5E-2, //1.0E-2,        // for minibatch, ModAdam (= ModAdamMax)
              //maxLR = 1.0E-2                  // for ModRAdamMax

    val ( (firstSGD, maxLRFor1stSGD),  (secondSGD, maxLRFor2ndSGD), (lambdaPhij, lambdaPhi) ) = protein match {

      case "PF00018uniq" => ( ("ModAdam", 0.006), ("Adam", 0.003), (2.0e-03, 5.0e-04) ) 
      case "PF00127uniq" => ( ("ModAdam", 0.003), ("Adam", 0.002), (4.0e-03, 4.0e-04) )
      case "PF00153uniq" => ( ("ModAdam", 0.006), ("Adam", 0.003), (1.4e-03, 2.5e-04) )
      case "PF00290uniq" => ( ("ModAdam", 0.003), ("Adam", 0.0006),(8.0e-04, 4.5e-04) )
      case "PF00565uniq" => ( ("ModAdam", 0.01),  ("Adam", 0.005), (1.7e-03, 3.5e-04) ) 
      case "PF00595uniq" => ( ("ModAdam", 0.006), ("Adam", 0.003), (3.0e-03, 5.0e-04) )
      case "PF00887uniq" => ( ("ModAdam", 0.006), ("Adam", 0.003), (2.0e-03, 2.5e-04) )
      case "PF00959uniq" => ( ("ModAdam", 0.005), ("Adam", 0.002), (6.0e-04, 1.0e-04) )
      case _ => {
                  ( ("ModAdam", 0.006), ("Adam", 0.003), (2.0e-03, 2.0e-04) ) 
                }
    }

  //val (nEff, nRep) = (effectiveNumberOfSamples, initialConfigurations.size)

    val fullBatchSize = initialConfigurations.size
    val miniBatchSize = { import scala.util.Sorting

                          val minSize = 100
                          if ( fullBatchSize <= minSize ) {
                              fullBatchSize
                          } else {
                              val r = Vector.range(minSize, (minSize * 1.1).toInt )
                              val t = r.map{i =>  (i, (fullBatchSize % i) / i.toDouble) }
                              t.fold(t(0)){ (x, y) => if (x._2 <= y._2) x else y }._1
                          }
                        }

    val miniBSize = BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
    val nMBsInFB = BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )
    val stepsPerEpoch = nMBsInFB

    val logInterval = {
                          val interval = ( (math.max(100, stepsPerEpoch) / stepsPerEpoch.toDouble + 0.5).floor * stepsPerEpoch).toInt
                          if (stepsPerEpoch <= 1) {
                                 if ( interval > 0 ) interval else 1
                          } else if ( stepsPerEpoch <= interval ) {
                                 (stepsPerEpoch * (interval / stepsPerEpoch.toDouble).floor ).toInt
                          } else {
                                 stepsPerEpoch
                          }
                      }

    val mBoltzman = new BM(ioFiles, stateOrderString,
                        effectiveNConfigs = effectiveNumberOfSamples, fia, fijab, 
                        regTerm = "GL1L2", propL1h = 0.0, propL1J = 1.0 )

                        //regTerm = "L1L2", propL1h = 0.0, propL1J = 0.9,
                        //regTerm = "L1L2", propL1h = 0.9, propL1J = 0.9,
                        //regTerm = "L2", propL1h = 0.0, propL1J = 0.0,

    val (betaVFor1stSGD, betaMFor1stSGD, epsFor1stSGD) = firstSGD match {
      case "Adam"     => (0.999, 0.9, 1.0e-6)
      case "ModAdam"  => (0.999, 0.9, 1.0e-6)
      case "RAdam"    => (0.999, 0.9, 1.0e-6)
      case "ModRAdam" => (0.999, 0.9, 1.0e-6)
      case "NAG"      => (0.95,  0.0, 0.0 ) // (0.9,  0.0, 0.0) // (0.99, 0.0, 0.0)  // for NAG: learningRate / (1 - betaV) should be constant.
      case "Adadelta" => (0.9, 0.9, 1.0e-5)
      case "RPROP-LR" => (0.5, 1.2, 1.0e-6)        // for RPROP-LR; betaV=rateDecrease , betaM=rateIncrease
      case _ => {
                  require( "" != "", s"${firstSGD} is not supported for SGD." ) 
                  (0.0, 0.0, 0.0)
                }
    }
    val (betaV, betaM, eps) = (betaVFor1stSGD, betaMFor1stSGD, epsFor1stSGD)

    val (betaVFor2ndSGD, betaMFor2ndSGD, epsFor2ndSGD) = secondSGD match {
      case "Adam"     => (0.999, 0.9, 1.0e-6)
      case "ModAdam"  => (0.999, 0.9, 1.0e-6)
      case "RAdam"    => (0.999, 0.9, 1.0e-6)
      case "ModRAdam" => (0.999, 0.9, 1.0e-6)
      case "NAG"      => (0.95,  0.0, 0.0 ) // (0.9,  0.0, 0.0) // (0.99, 0.0, 0.0)  // for NAG: learningRate / (1 - betaV) should be constant.
      case "Adadelta" => (0.9, 0.9, 1.0e-5)
      case "RPROP-LR" => (0.5, 1.2, 1.0e-6)        // for RPROP-LR; betaV=rateDecrease , betaM=rateIncrease
      case _ => {
                  require( "" != "", s"${secondSGD} is not supported for SGD." ) 
                  (0.0, 0.0, 0.0)
                }
    }

    val gradientDescentMethod = firstSGD  // + secondSGD by jobControlParam.dat    // NAG, Adadelta, RPROP-LR, RAdam, Adam, ModRAdam, ModAdam

    val maxLR = maxLRFor1stSGD    //0.5E-02 0.6E-02 // 1.0E-2  // 5.0E-2,        // for ModAdam (= ModAdamMax)

    val warmupEpochs = 100.0
    val minLearningSteps = {
                           //val m = (2200.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
                           //val m = (1800.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
                             (1800 / logInterval.toDouble).ceil.toInt * logInterval
                           }

    val maxLREpochs = minLearningSteps - warmupEpochs
    assert( minLearningSteps > warmupEpochs )
  //val minLearningStepsForBestKL = minLearningSteps + (500.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
  //val minLearningStepsForBestKL = minLearningSteps + (1000.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
    val minLearningStepsForBestKL = minLearningSteps + (2000.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
    val maxNoLearnings = (200.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch 

    val learningRate = new LearningRate(warmupEpochs = warmupEpochs,                //minLearningSteps * 0.05,
                     maxLR = maxLR, maxLREpochs = maxLREpochs,  //minLearningSteps * 0.95, //scala.Int.MaxValue.toDouble , // minLearningSteps.toDouble,
                     coolingMethod = "(1+at)^b", constPerEpoch = 0.01, power = -0.5,
                   //coolingMethod = "1/(1+at)", epochsUntilHalf = 10.0,
                   //coolingMethod = "rate", coolRate = 0.99, minLR = maxLR * 0.01,
                     stepsPerEpoch = 1.0 )

  // The second SGD method following the first SGD method defined in the RunBM(...).
    val minLearningStepsFor2ndSGD = {
                           //(1800.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
                           //(2800.0/stepsPerEpoch).ceil.toInt * stepsPerEpoch
                             (1800 / logInterval.toDouble).ceil.toInt * logInterval
                           //(2800 / logInterval.toDouble).ceil.toInt * logInterval
                           }

    val outJobControlParamFile = new PrintStream(ioFiles("jobControlParamFile") )
    outJobControlParamFile.print(s"minLearningSteps\t${minLearningStepsFor2ndSGD}\n")
    outJobControlParamFile.print(s"gradientDescentMethod\t${secondSGD}\n")
    outJobControlParamFile.print(s"betaV\t${betaVFor2ndSGD}\n")
    outJobControlParamFile.print(s"betaM\t${betaMFor2ndSGD}\n")
    outJobControlParamFile.print(s"maxLR\t${maxLRFor2ndSGD}\n")  //("maxLR\t0.1E-02\n")  //("maxLR\t0.15E-02\n")  //("maxLR\t0.2E-2\n")  //("maxLR\t0.25E-02\n")  //("maxLR\t0.5E-2\n")
    outJobControlParamFile.close()

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
                maxExtendedIterations = 1,
                mcmcKernel = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

                lambdaPhi = lambdaPhi, lambdaPhij = lambdaPhij,

                gauge = "ungauged",  // "ungauged" or "phi_zeroSum" for L2 and L2GL1, "ungauged" for all
                gradientDescentMethod = gradientDescentMethod,

                betaV = betaV, betaM = betaM, eps = eps, //betaV = 0.999, betaM = 0.9, eps = 1.0e-6,    // for ModAdam (= ModAdamMax)

                learningRate = learningRate,

                pseudoNCountsForBM = 0.0,
                pseudoNCountsForhia = 10.0,
                pseudoNCountsForProposedPia = 10.0,
                pseudoNCountsForKL = 1.0,

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

