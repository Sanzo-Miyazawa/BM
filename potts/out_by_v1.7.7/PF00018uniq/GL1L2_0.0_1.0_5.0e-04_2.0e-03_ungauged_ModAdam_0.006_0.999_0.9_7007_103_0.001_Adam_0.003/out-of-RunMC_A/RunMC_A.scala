#! /bin/sh
# "$0" is passed as the first argument.
if [ "$SCALA" = "" ] ; then
        exec scala -savecompiled "$0" "$0" "$@"
else
        exec $SCALA -savecompiled "$0" "$0" "$@"
fi

!#


object RunMC {
  import java.io.File
  import java.io.PrintStream
  import scala.collection.immutable.Seq
  import scala.collection.immutable.HashMap
  import scala.collection.immutable.ArraySeq
  import scala.collection.immutable.Vector
  import scala.collection.parallel.immutable.ParVector

  //import scala.collection.parallel.CollectionConverters._     // for par in scala 2.13
  import scala.collection.parallel.CollectionConverters.*       // for par in scala 3
  import scala.math.Ordering.Double.TotalOrdering               // for scala 2.13

//import org.biojava.nbio.core.util.InputStreamProvider

  import org.sanzo.potts.MCMC
  import org.sanzo.potts.MCMC.readFasta
  import org.sanzo.potts.MCMC.readInteractions
  import org.sanzo.potts.MCMC.Interactions
  import org.sanzo.potts.MCMC.State
  import org.sanzo.potts.MCMC.EnsembleAverages
  import org.sanzo.potts.MCMC.printPiPij
  import org.sanzo.potts.BM
  import org.sanzo.potts.BM.BMState
  import org.sanzo.potts.BM.bayesianCorrection
//import org.sanzo.potts.BM.calcKL
  import org.sanzo.potts.BM.printKL
  import org.sanzo.potts.BM.printLog
  import org.sanzo.potts.BM.miniBatchSizeCorrected
  import org.sanzo.potts.BM.nMiniBatchesInFullBatch
  import org.sanzo.potts.BM.nMiniBatchesForEnsembleAverage
  import org.sanzo.potts.BM.startPosInFB
  import org.sanzo.potts.BM.endPosOfMiniBatch
  import org.sanzo.potts.Util.toProb
  import org.sanzo.potts.zip.gunzipInputStream

  import breeze.linalg.DenseVector
  import breeze.linalg.DenseMatrix
  import breeze.stats.distributions.Rand
  import breeze.stats.distributions.RandBasis

  def processFiles(args: Array[String], 
        ioFiles: HashMap[String, java.io.File]) = {

// Use scala.util.Random.setSeed(seed: Long), if you want consistent behavior with the same random numbers from execution to execution.
    val seedForRandom = 19L
    scala.util.Random.setSeed(seedForRandom)
    scala.sys.process.stdout.print(s"# seed for scala.util.Random = ${seedForRandom}\n")
// Use org.sanzo.potts.SetSeedForRand.randBasis instead of breeze.stats.distributions.Rand.Fixed.randBasis
    org.sanzo.potts.SetSeedForRand(0)
    scala.sys.process.stdout.print(s"# seed for org.sanzo.potts.SetSeedForRand.randBasis = ${org.sanzo.potts.SetSeedForRand.seed}\n")

    val (effectiveNConfigs, fia, fijab, stateOrderString) =
      if ( args.size < 1 ) {
        require( args.size > 0)
        (0.0, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]](), "")
      } else {
      //val file = new java.io.File( args(0) )
        val file = ioFiles("fiafijabFile")
        val fin = gunzipInputStream(file)

        val bufSource = scala.io.Source.fromInputStream(fin)
        val lines = bufSource.getLines()
        val (effectiveNumberOfSamples, fia, fijab, stateOrderString) = MCMC.readPiPij(lines)
        bufSource.close()
        (effectiveNumberOfSamples, fia, fijab, stateOrderString)
      }

    val fa = toProb(fia)    

    val nSamples = 1

  //val initialMCStates: MCMC.InitialMCStates  = {
    val nativeConfigs = {
      //val fastaFile = new java.io.File( args(1) )
        val fastaFile = ioFiles("configFile")
        val fastaIn = gunzipInputStream(fastaFile)
        val fastaLines = scala.io.Source.fromInputStream(fastaIn).getLines()

        val (idArraySeq, seqArraySeq, numArraySeq) =
            MCMC.readFasta(fastaLines, stateOrderString = stateOrderString, caseSensitive = false)

        fastaIn.close()

        assert(idArraySeq.size == seqArraySeq.size && idArraySeq.size == numArraySeq.size )

        val initialConfigurations = numArraySeq

        val initialConfIDs = idArraySeq.map( id => id.split("[ \t]")(0) )

        val miniBatchSize = initialConfigurations.size
        val fullBatchSize = miniBatchSize
        val miniBSize = miniBatchSize //BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
        val nMBsInFB = 1              //BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )

        val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) * nMBsInFB  + 0.5 ).floor.toInt
        val nMBsForEnsembleAve = BM.nMiniBatchesForEnsembleAverage(miniBSize, math.max(fullBatchSize, 50000), nMiniBatchesForEnsembleAve )
        
        val multipleInitialConfigurations = ArraySeq.range(0, nMBsForEnsembleAve).map{ a => initialConfigurations }.flatten
        val multipleInitialConfIDs = ArraySeq.range(0, nMBsForEnsembleAve).map{ i => initialConfIDs.map{id => id + s"_${i}"} }.flatten

        ArraySeq.range(0, multipleInitialConfigurations.size).map{ i => MCMC.NativeConfig( i , multipleInitialConfIDs(i), multipleInitialConfigurations(i) ) }
        
     // (multipleInitialConfigurations , multipleInitialConfIDs )
     // val initialMCStates = initialConfigurations.map{ c => MCMC.State( configuration = c )
    }

    val mcmcInteractions: MCMC.Interactions = {

      //val file2 = new java.io.File( args(2) )
        val file2 = ioFiles("interactionFile")
        val fin2 = gunzipInputStream(file2)

        val bufSource2 = scala.io.Source.fromInputStream(fin2)
        val lines2 = bufSource2.getLines()

      //val (hia, jijab, stateOrderString2) = MCMC.readhJ(lines2)
        val (interactions, stateOrderString2) = MCMC.readInteractions(lines2)

        bufSource2.close()

        require( stateOrderString == stateOrderString2)
      //MCMC.Interactions(hia, jijab)
        interactions
    }

//
//

    val initialMCStates: ParVector[MCMC.State] = {
        val configs = nativeConfigs.map{ n => MCMC.State( nativeConfig = n,
                                                        configuration = n.sample, interactions = mcmcInteractions ) }
        Vector(configs*).par
    }

    val nUnits = fia.size
    val nStatesOfUnit = fia(0).size

    val nIndependentMC = initialMCStates.size

    val miniBatchSize = initialMCStates.size
    val fullBatchSize = miniBatchSize
    val miniBSize = miniBatchSize //BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
    val nMBsInFB = 1              //BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )


    val logInterval = {
                          val stepsPerEpoch = nMBsInFB
                          val interval = ( (math.max(100, stepsPerEpoch) / stepsPerEpoch.toDouble + 0.5).floor * stepsPerEpoch).toInt
                          if (stepsPerEpoch <= 1) {
                                 if ( interval > 0 ) interval else 1
                          } else if ( stepsPerEpoch <= interval ) {
                                 (stepsPerEpoch * ((interval / stepsPerEpoch.toDouble).floor) ).toInt
                          } else {
                                 stepsPerEpoch
                          }
                      }

//

  //val pseudoNCounts = 10.0
    val pseudoNCountsForBM = 0.0
  //val pseudoNCountsForhia = 10.0
    val pseudoNCountsForProposedPia = 10.0
    val pseudoNCountsForKL = 1.0

    val nInitialIterationsPerUnit = 0
  //val everyNIterationsPerUnit = 10
    val everyNIterationsPerUnit = 1000
    val maxExtendedIterations = 1
//

    val stepsPerEpoch = nMBsInFB.toDouble

    val logInterval_rev = if (stepsPerEpoch <= 1.0) {
                                 if ( logInterval > 0 ) logInterval else 1
                              } else if ( stepsPerEpoch <= logInterval ) {
                                 (stepsPerEpoch * (logInterval / stepsPerEpoch).floor ).toInt
                              } else {
                                 stepsPerEpoch.toInt
                              }

  /*
    val (interactions, bmState) = BM.initialize(fia, mcmcInteractions, 
                  // The following arguments are dummy.
                  optMethod = "Adam",
                  learningRate = new LearningRate(), learningRateForRPROPLR = new LearningRateForRPROPLR(),
                  betaV = 0.0, betaM = 0.0 )
                  //
    assert(interactions == mcmcInteractions)
  */
    val bmState = BM.BMState()

    val (fiaForBM, fijabForBM, faForBM,  fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1) =
            BM.freqsForBMandKL( fia, fijab, fa,
                effectiveNConfigs,
                pseudoNCountsForBM,        // If the regularization is not used, this may be used to avoid 0 counts in this.fia and this.fijab
                pseudoNCountsForKL         // used to calculate KL valunes from ensemble averages of Pij(a,b)
            )

    val calcKL = BM.calcKL(fiaForKL, fijabForKL, observedNP = effectiveNConfigs, pseudoNP = 0.0,
                               _: ArraySeq[DenseVector[Double]], _: ArraySeq[DenseMatrix[Double]], _: Double, _: Double,
                               fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1 )

    val proposedPia = BM.bayesianCorrection( fia, effectiveNConfigs, pseudoNCountsForProposedPia )

    val mcmc = new MCMC( ioFiles, stateOrderString,
                mcmcInteractions, proposedPia = proposedPia ) 

    val mcmcKernel = "GibbsWithMHStep"

    val kernelForMCMC: (MCMC.State, RandBasis) => Rand[MCMC.State] =
          if ( mcmcKernel == "MH" ) {
             // mcmc.kernelMH(_)(_)
                mcmc.kernelMH(_)(using _: RandBasis)
          } else if (  mcmcKernel == "GibbsWithMHStep" || mcmcKernel == "MultiBlockMH" ) {
             // mcmc.kernelGibbsWithMHStep(_)(_)
                mcmc.kernelGibbsWithMHStep(_)(using _: RandBasis)
          } else {      // if ( mcmcKernel == "Gibbs" ) {
             // mcmc.kernelGibbs(_)(_)
                mcmc.kernelGibbs(_)(using _: RandBasis)
          }

    val (revInitialStates, nNonEquil, nExtendedIterations, independentSamplings) =
          mcmc.runMC(
                initialStates = initialMCStates,              
                nInitialIterationsPerUnit = nInitialIterationsPerUnit,
                everyNIterationsPerUnit = everyNIterationsPerUnit,
                nSamples = nSamples,
                initialT = 1.0,
                annealingRate = 0.99,
                finalT = 1.0,
                maxExtendedIterations = maxExtendedIterations,
              //kernel = mcmc.kernelGibbsWithMHStep
                kernel = kernelForMCMC
                )

  /*
    val fileIndpendentMCsamplings = new java.io.File(ioFiles("outputDir"), s"out_samples.fasta.gz")
    val outSamples =  new java.io.PrintStream(fileIndpendentMCsamplings)

    MCMC.printIndependentMCsamplings(outSamples, stateOrderString,
                independentSamplings )
  */

    val pia = mcmc.frequenciesAtUnitInSamples( independentSamplings )
    val pijab = mcmc.pairwiseFrequenciesInSamples( independentSamplings)
    val ensembleAverages = MCMC.EnsembleAverages( pia = pia , pijab = pijab )

    val bmStateWithMCSamples = bmState.copy( step = 0, initialMCStates = revInitialStates,
              independentSamplings = independentSamplings, ensembleAverages = ensembleAverages
              )

  //val lastMiniBatchSamplings = independentSamplings
    type IndependentSamplings = MCMC.IndependentSamplings            //ParVector[ArraySeq[MCMC.State]]

    val nMiniBatchesForEnsembleAve = 1
    val nMBsForEnsembleAve = nMiniBatchesForEnsembleAve
    
    val miniBacthesForEnsembleAve: Vector[BMState] = Vector.fill(1)(bmStateWithMCSamples)
    val revisedLastAndCurrentFullBatch: Vector[BMState] = Vector.fill(1)(bmStateWithMCSamples)

    val fileTE = new File(ioFiles("outputDir"),s"out_Energy_distribution_of_MC_samples.txt")
       //(f"Energy_distribution_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outTE  = new PrintStream(fileTE)
    val fileKL = new File(ioFiles("outputDir"),s"out_KL_of_MC_samples.txt")
       //      (f"KL_of_each_step_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outKL  = new PrintStream(fileKL)

    val totalNSamples = (independentSamplings.size * independentSamplings(0).size).toDouble

    val ensembleKL = calcKL( ensembleAverages.pia, ensembleAverages.pijab, totalNSamples, pseudoNCountsForKL )

    outKL.print(s"# nIndependentMC= ${nIndependentMC}  nSamples for 1 MC = ${nSamples}\n")
    outKL.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outKL.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${ioFiles("interactionFile").getName()}\n")

    outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and those over %d minibatch ensembles with pseudoNCountsForKL= %g\n".format(
                    nMiniBatchesForEnsembleAve, pseudoNCountsForKL ) )

    BM.printKL(outKL, bmStep = bmStateWithMCSamples.step, nNonEquil = nNonEquil, nExtendedIterations = nExtendedIterations, ensembleKL = ensembleKL, aveKL = ensembleKL )

    outTE.print(s"# nIndependentMC= ${nIndependentMC}  nSamples/MC= ${nSamples}\n")
    outTE.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outTE.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${ioFiles("interactionFile").getName()}\n#\n")
    outTE.print("# TE is calculated in the Ising gauge for comparison.\n#\n")
    outTE.print("# For <TE>_m and <(TE -<TE>)^2>_m, %d minibatch ensembles are employed for averaging.\n#\n".format(nMBsForEnsembleAve))

    val logFiles: Option[Tuple2[Option[java.io.File], Option[java.io.File]]] =
                  if ( (bmStateWithMCSamples.step + 1) % logInterval_rev == 0 || bmStateWithMCSamples.step == 0 ) {
                    val logFileIndpendentMCsamplings = new File(ioFiles("outputDir"), s"out_MC_samples_by_MCMC_0.fasta.gz")
                    Option( (None, Option(logFileIndpendentMCsamplings)) )
                  //val logFileInteractions = new File(ioFiles("outputDir"), s"out_Interactions_hJ_for_MCMC_0.txt.gz")
                  //Option( (Option(logFileInteractions), Option(logFileIndpendentMCsamplings)) )
                  } else {
                    None
                  }

    val (interactions_Ising, vecIndependentSamplings_Ising_ForEnsembleAve,
                     ensembleTE_MeanAndVariance, sampleTE_MeanAndVariance, allRandomizedSampleTE_MeanAndVariance,
                     allEnsembleTE_MeanAndVariance, allSampleTE_MeanAndVariance, randomSampleTE_MeanAndVariance) = {
                        val interactions = mcmc.interactions
                        BM.printLog( bmStateWithMCSamples.step, logFiles, mcmc.stateOrderString,
                          interactions, miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch,
                          outTE, faForBM, startPos = 0, miniBatchSize = miniBSize, nMiniBatchesForEnsembleAve = nMBsForEnsembleAve )
                     }

    outKL.close()
    outTE.close()
  }

  def main(argsP: Array[String]): Unit = {
    //val outfile = new File("/tmp/output.txt")
    //val fout = new FileOutputStream(outfile)

      scala.sys.process.stdout.print("#")
      argsP.foreach( i => scala.sys.process.stdout.print(" " + i) )
      scala.sys.process.stdout.print("\n")

      val args = argsP.drop(1)
      val ioFiles = HashMap( 
                ("outputDir", new java.io.File(args(0)) ),
                ("fiafijabFile", new java.io.File(args(1)) ),
                ("configFile", new java.io.File(args(2)) ),
                ("interactionFile", new java.io.File(args(3)) )
                )
      args(1) match {
        case _ => {
                //val fout = scala.sys.process.stdout   //System.out 
                processFiles(args.drop(1), ioFiles)
                }
      }
   }

}

