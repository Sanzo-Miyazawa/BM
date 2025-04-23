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

  import scala.annotation.tailrec
  import scala.util.Random

  import org.biojava.nbio.core.util.InputStreamProvider

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
  import org.sanzo.potts.BM.freqsForBMandKL
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

    val (effectiveNumberOfSamples, fia, fijab, stateOrderString) =
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
    val nUnits = fia.size
    val nStatesOfUnit = fia(0).size

    //val (initialConfIDs: ArraySeq[String], initialConfigurations: ArraySeq[ArraySeq[Byte]]) = {
    val nativeConfigs = {
     // val fastaFile = new java.io.File( args(1) )
        val fastaFile = ioFiles("configFile")
     // val fastaIn = isp.getInputStream(fastaFile)
        val fastaIn = gunzipInputStream(fastaFile)
        val fastaLines = scala.io.Source.fromInputStream(fastaIn).getLines()

        val (idArraySeq, seqArraySeq, numArraySeq) =
            MCMC.readFasta(fastaLines, stateOrderString = stateOrderString, caseSensitive = false)

        fastaIn.close()

        assert(idArraySeq.size == seqArraySeq.size && idArraySeq.size == numArraySeq.size )

        val initialConfigurations = numArraySeq

        val initialConfIDs = idArraySeq.map( id => id.split("[ \t]")(0) )

        ArraySeq.range(0, initialConfigurations.size).map{ i => MCMC.NativeConfig( index = i , id = initialConfIDs(i), sample = initialConfigurations(i) ) }
    }

    val mcmcInteractions: MCMC.Interactions = {

      //val file2 = new java.io.File( args(2) )
        val file2 = ioFiles("interactionFile")
      //val fin2 = isp.getInputStream(file2)
        val fin2 = gunzipInputStream(file2)

        val bufSource2 = scala.io.Source.fromInputStream(fin2)
        val lines2 = bufSource2.getLines()

      //val (hia, jijab, stateOrderString2) = MCMC.readhJ(lines2)
        val (interactions, stateOrderString2) = MCMC.readInteractions(lines2)
        require( stateOrderString == stateOrderString2)

        bufSource2.close()
      //MCMC.Interactions(hia, jijab)
        interactions
    }

    val fullBatchSize = nativeConfigs.size
  /*
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
  */
    val miniBatchSize = fullBatchSize

    val miniBSize = BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
    val nMBsInFB = BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )

    val nSamples = 1                // if nSamples == 0, initialNSamples, incrementSamples, and maxNSamples are used.

  //val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) * nMBsInFB  + 0.5 ).floor.toInt
  //val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) + 0.5 ).floor.toInt * nMBsInFB
    val nMiniBatchesForEnsembleAve = 1

    val stepsPerEpoch = nMBsInFB.toDouble
  //val logInterval = ( (math.max(100.0, stepsPerEpoch) / stepsPerEpoch + 0.5).floor * stepsPerEpoch).toInt
    val logInterval = 1

  //val results = runMCMCs(
    runMCMCs(
                ioFiles = ioFiles,
                stateOrderString = stateOrderString,

                effectiveNConfigs = effectiveNumberOfSamples,
                fia = fia,
                fijab = fijab,

                mcmcInteractions = mcmcInteractions,
                nativeConfigs = nativeConfigs,

                nInitialIterationsPerUnit = 0,
              //everyNIterationsPerUnit = 10,
                everyNIterationsPerUnit = 100,

                miniBatchSize = miniBSize,             // if miniBatchSize > 0, then nSamples = 1

                nSamples = nSamples,                   // if nSamples == 0, initialNSamples, incrementSamples, and maxNSamples are used.

                initialT = 1.0,         // 1.2
                finalT = 1.0,
                annealingRate = 0.99,
                maxExtendedIterations = 1,
                mcmcKernel = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

                pseudoNCountsForBM = 0.0,
                pseudoNCountsForhia = 10.0,
                pseudoNCountsForProposedPia = 10.0,
                pseudoNCountsForKL = 1.0,

                nMiniBatchesForEnsembleAve = 1,
                
                nSteps = 10,

                logInterval = logInterval
                )

  }

  def runMCMCs(
    ioFiles: HashMap[String, File],
    stateOrderString: String,

    effectiveNConfigs: Double,
    fia: ArraySeq[DenseVector[Double]], 
    fijab: ArraySeq[DenseMatrix[Double]],

    mcmcInteractions: MCMC.Interactions,
    nativeConfigs: ArraySeq[MCMC.NativeConfig],

    nInitialIterationsPerUnit: Int = 0,
    everyNIterationsPerUnit: Int = 10,

    miniBatchSize: Int = 100,      // if miniBatchSize != 0, then initialConfigurations.size >= miniBatchSize

    nSamples: Int = 1,

    initialT: Double = 1.0,         // 1.2
    finalT: Double = 1.0,
    annealingRate: Double = 0.99,
    maxExtendedIterations: Int = 1,
    mcmcKernel: String = "GibbsWithMHStep",   //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

    pseudoNCountsForBM: Double = 0.0,
    pseudoNCountsForhia: Double = 10.0,
    pseudoNCountsForProposedPia: Double = 10.0,
    pseudoNCountsForKL: Double = 1.0,

    nMiniBatchesForEnsembleAve: Int = 1,    // The ensemble average of TE will be done over max nMiniBatchesForEnsembleAve minibatch ensembles; 0 means over-epoch.

    nSteps: Int = 10,

    logInterval: Int = 100
  ) = {

  //val initialConfigurations = nativeConfigs.map{ native => native.sample }
  //val nIndependentMC = initialMCStates.size

    val fullBatchSize = nativeConfigs.size
    val miniBSize = BM.miniBatchSizeCorrected( miniBatchSize, fullBatchSize )
    val nMBsInFB = BM.nMiniBatchesInFullBatch( miniBSize, fullBatchSize )

    require( fullBatchSize >= miniBSize )

    val nMBsForEnsembleAve = BM.nMiniBatchesForEnsembleAverage(miniBSize, math.max(fullBatchSize, 50000), nMiniBatchesForEnsembleAve )
    require( nMBsForEnsembleAve == 1, "Required: nMiniBatchesForEnsembleAve == 1 for this program\n" )           

    val nFBs = (nMBsForEnsembleAve.toDouble / nMBsInFB).ceil.toInt
    assert( nFBs == 1 )

    val stepsPerEpoch = nMBsInFB.toDouble
    val logInterval_rev = if (stepsPerEpoch <= 1.0) {
                                 if ( logInterval > 0 ) logInterval else 1
                              } else if ( stepsPerEpoch <= logInterval ) {
                                 (stepsPerEpoch * (logInterval / stepsPerEpoch).floor ).toInt
                              } else {
                                 stepsPerEpoch.toInt
                              }

    val nativeMCStates: ArraySeq[MCMC.State] = {
        val configs = nativeConfigs.map{ n => MCMC.State( nativeConfig = n,
                                                        configuration = n.sample, interactions = mcmcInteractions ) }
        configs
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

    val fa = toProb(fia)

  //val (fiaWithPseudo, fijabWithPseudo) = BM.bayesianCorrection(fia, fijab, effectiveNConfigs, pseudoNCounts)
    val (fiaForBM, fijabForBM, faForBM,  fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1) =
            BM.freqsForBMandKL( fia, fijab, fa,
                effectiveNConfigs,
                pseudoNCountsForBM,        // If the regularization is not used, this may be used to avoid 0 counts in fia and fijab
                pseudoNCountsForKL         // used to calculate KL valunes from ensemble averages of Pij(a,b)
            )

    val calcKL = BM.calcKL(fiaForKL, fijabForKL, observedNP = effectiveNConfigs, pseudoNP = 0.0,
                               _: ArraySeq[DenseVector[Double]], _: ArraySeq[DenseMatrix[Double]], _: Double, _: Double,
                               fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1 )

    val proposedPia = BM.bayesianCorrection(fia, effectiveNConfigs, pseudoNCountsForProposedPia)

    val fileTE = new File(ioFiles("outputDir"),s"out_Energy_distribution_of_MC_samples.txt")
       //(f"Energy_distribution_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outTE  = new PrintStream(fileTE)

    val fileKL = new File(ioFiles("outputDir"),s"out_KL_of_MC_samples.txt")
       //      (f"KL_of_each_step_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outKL  = new PrintStream(fileKL)

    outKL.print(s"# miniBatchSize= ${miniBSize}  fullBatchSize= ${fullBatchSize}  nSamples/MC= ${nSamples}\n")
    outKL.print(s"# pseudoNCountsForProposedPia/KL/BM/hia= ${pseudoNCountsForProposedPia} / ${pseudoNCountsForKL} / ${pseudoNCountsForBM} / ${pseudoNCountsForhia}}\n")
    outKL.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outKL.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${ioFiles("interactionFile").getName()}\n")
    outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and %d minibatch ensembles with pseudocountsForKL= %g\n".format(
                    nMBsForEnsembleAve, pseudoNCountsForKL ) )

    outTE.print(s"# miniBatchSize= ${miniBSize}  fullBatchSize= ${fullBatchSize}  nSamples/MC= ${nSamples}\n")
    outTE.print(s"# pseudoNCountsForProposedPia/KL/BM/hia= ${pseudoNCountsForProposedPia} / ${pseudoNCountsForKL} / ${pseudoNCountsForBM} / ${pseudoNCountsForhia}}\n")
    outTE.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outTE.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${ioFiles("interactionFile").getName()}\n#\n")
    outTE.print("# TE is calculated in the Ising gauge for comparison.\n#\n")
    outTE.print("# For <TE>_m and <(TE -<TE>)^2>_m, %d minibatch ensembles are employed for averaging.\n#\n".format(nMBsForEnsembleAve))

  //outTE.print("# step\tL=%d\t(all randomized sample mean - all randomized sample variance/L)  TE_/L  <TE>/L\tall randomized sample mean\tall randomized sample variance/L  TE sample variance/L  <(TE-<TE>)^2>/L\n".format(nUnits) )
  //outTE.print("# step\tL=%d\t(random sample mean - random sample variance)/L  TE_all/L  <TE>/L  <TE>_m/L\trandom sample mean\trandom sample variance/L  TE all sample variance/ L  <(TE-<TE>)^2>/L  <(TE-<TE>)^2>_m/L ; <...>_m over max %d minibatch ensembles\n".format(nUnits, nMiniBatchesForEnsembleAve) )

    @annotation.tailrec
    def iterate(i_th: Int, nIterations: Int, nFullBatches: Vector[ArraySeq[MCMC.State]] ): Unit = {

    val nFBsOfBMStatesIntInt: Vector[ArraySeq[Tuple3[BM.BMState, Int, Int ]]] = {

        assert( nFullBatches.size == nFBs )
        val delta = nFBs * nMBsInFB  - nMBsForEnsembleAve

        Vector.range(0 , nFBs).map{ iFB =>

            val mcStates = Random.shuffle( nFullBatches(iFB) ).toVector

            assert(mcStates.size == fullBatchSize)

            ArraySeq.range(0 , nMBsInFB ).map{iMB =>

                val startPos = startPosInFB(iMB, miniBSize, fullBatchSize)
                val ( _, endPos, untilEndPos)  = endPosOfMiniBatch(startPos, miniBSize, fullBatchSize)

                val step = if ( iFB * nMBsInFB + iMB < delta ) {
                             if(  nFBs <= 1 )
                               iFB * nMBsInFB + iMB
                             else
                             //bmState.step      // indicates this is not a target for MCMC.
                               iFB * nMBsInFB + iMB
                           } else {
                             iFB * nMBsInFB + iMB
                           }

                if ( step < 0 ) {
                  (bmState, 0, 0)
                } else {
                  val initialMCStates = mcStates.slice(startPos, untilEndPos).par

                  val mcmc = new MCMC( ioFiles, stateOrderString,
                     mcmcInteractions, proposedPia = proposedPia ) 

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

                  val (revInitialStates: ParVector[MCMC.State], nNonEquil, nExtendedIterations, independentSamplings: ParVector[ArraySeq[MCMC.State]]) =
                    mcmc.runMC(
                      initialStates = initialMCStates,
                      nInitialIterationsPerUnit = nInitialIterationsPerUnit,
                      everyNIterationsPerUnit = everyNIterationsPerUnit,
                      nSamples = nSamples,
                      initialT = initialT,
                      annealingRate = annealingRate,
                      finalT = finalT,
                      maxExtendedIterations = maxExtendedIterations,
                    //kernel = mcmc.kernelGibbsWithMHStep     //mcmc.kernelGibbsWithMHStep==mcmc.kernelMultiBlockMH
                      kernel = kernelForMCMC
                      )

                  val pia = mcmc.frequenciesAtUnitInSamples( independentSamplings )
                  val pijab = mcmc.pairwiseFrequenciesInSamples( independentSamplings )
                  val ensembleAverages = MCMC.EnsembleAverages( pia = pia , pijab = pijab )

                  val revStep = step + i_th * nFBs * nMBsInFB   

                  val bmStateWithMCSamples = bmState.copy( step = revStep,
                          initialMCStates = revInitialStates, independentSamplings = independentSamplings, ensembleAverages = ensembleAverages )

                  (bmStateWithMCSamples, nNonEquil, nExtendedIterations )
                }

            }
        }
    }

    val nFBsOfBMStates = nFBsOfBMStatesIntInt.map{ FB => FB.map{ t3 => t3._1 } }
    val nFBsOfNnonEquiv = nFBsOfBMStatesIntInt.map{ FB => FB.map{ t3 => t3._2 } }
    val nFBsOfNExtend = nFBsOfBMStatesIntInt.map{ FB => FB.map{ t3 => t3._3 } }

    val bmStep = nFBsOfBMStates(nFBs - 1)(nMBsInFB - 1).step 
    val startPos = startPosInFB( nMBsInFB - 1, miniBatchSize, fullBatchSize )

    val nMBsOfBMState = nFBsOfBMStates.flatten
    val nMBsOfNnonEquiv = nFBsOfNnonEquiv.flatten
    val nMBsOfNExtend = nFBsOfNExtend.flatten

    assert( nFBsOfBMStates.size == 1 )

    assert( nFBs == nFBsOfBMStates.size )
    assert( nMBsInFB == nFBsOfBMStates(0).size )

    val (miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch) = if(  nFBs <= 1 ) {
            (nFBsOfBMStates(0).slice(nMBsInFB - nMBsForEnsembleAve, nMBsInFB).toVector, nMBsOfBMState) 
          } else {
            (nMBsOfBMState.slice( nFBs * nMBsInFB - nMBsForEnsembleAve, nFBs * nMBsInFB ), nFBsOfBMStates(nFBs - 1).toVector )
          }


  //val fa = toProb(fia)

    for ( iMB <- Range(0, nMBsOfBMState.size) ) {
        if ( nMBsOfBMState(iMB).step >= 0 ) {

          val nNonEquil = nMBsOfNnonEquiv(iMB)
          val nExtendedIterations = nMBsOfNExtend (iMB) 

          val ensembleKL = {
              val ensembleAverages = nMBsOfBMState(iMB).ensembleAverages
              val independentS = nMBsOfBMState(iMB).independentSamplings
              val totalNSamples = (independentS.size * independentS(0).size).toDouble
              calcKL( ensembleAverages.pia, ensembleAverages.pijab, totalNSamples, pseudoNCountsForKL )
          }

          val bmStep = nMBsOfBMState(iMB).step

          val miniBatches = nMBsOfBMState.slice( math.max(0, iMB - nMBsForEnsembleAve), iMB + 1) 

          val aveKL = if ( iMB <= 0 ) {
                        ensembleKL
                      } else {
                        val independentS: ParVector[ArraySeq[MCMC.State]] =
                                miniBatches.map{bms => bms.independentSamplings}.flatten.par
                        val nStatesOfUnit = fia(0).size
                        val pia = MCMC.frequenciesAtUnitInSamples( independentS, nStatesOfUnit )
                        val pijab = MCMC.pairwiseFrequenciesInSamples( independentS, nStatesOfUnit )
                        val totalNSamples = (independentS.size * independentS(0).size).toDouble

                      //BM.calcKL(fiaWithPseudo, fijabWithPseudo, observedNP = Double.PositiveInfinity, pseudoNP = 0.0,
                      //  pia, pijab, observedNQ = totalNSamples, pseudoNQ = pseudoNCounts)
                        calcKL( pia, pijab, totalNSamples, pseudoNCountsForKL)

                      }

          BM.printKL(outKL, bmStep = bmStep, nNonEquil = nNonEquil, nExtendedIterations = nExtendedIterations, ensembleKL = ensembleKL, aveKL = aveKL )

          {
            val revisedLastAndCurrentFullBatch = {
                val fb = nMBsOfBMState.slice( (iMB / nMBsInFB) * nMBsInFB, (iMB / nMBsInFB + 1) * nMBsInFB )
                fb.zipWithIndex.map{ (bms, i) => if ( i <= iMB % nMBsInFB ) bms else bms.copy(step = scala.Int.MinValue) }
            }
            val startPos = startPosInFB( iMB % nMBsInFB, miniBatchSize, fullBatchSize )

            val miniBacthesForEnsembleAve = miniBatches
            val bmStateWithMCSamples = miniBacthesForEnsembleAve(miniBacthesForEnsembleAve.size - 1)
            val independentSamplings = bmStateWithMCSamples.independentSamplings

            val logFiles: Option[Tuple2[Option[java.io.File], Option[java.io.File]]] =
                  if ( (bmStateWithMCSamples.step + 1) % logInterval_rev == 0 || bmStateWithMCSamples.step == 0 ) {
                    val logFileIndpendentMCsamplings = new File(ioFiles("outputDir"), s"out_MC_samples_by_MCMC_${bmStep}.fasta.gz")
                    Option( (None, Option(logFileIndpendentMCsamplings)) )
                  //val logFileInteractions = new File(ioFiles("outputDir"), s"out_Interactions_hJ_for_MCMC_${bmStep}.txt.gz")
                  //Option( (Option(logFileInteractions), Option(logFileIndpendentMCsamplings)) )
                  } else {
                    None
                  }

            val (interactions_Ising, vecIndependentSamplings_Ising_ForEnsembleAve,
                     ensembleTE_MeanAndVariance, sampleTE_MeanAndVariance, allRandomizedSampleTE_MeanAndVariance,
                     allEnsembleTE_MeanAndVariance, allSampleTE_MeanAndVariance, randomSampleTE_MeanAndVariance) =
                        BM.printLog( bmStateWithMCSamples.step, logFiles, stateOrderString,
                          independentSamplings(0)(0).interactions, miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch,
                          outTE, fa, startPos = startPos, miniBatchSize = miniBSize, nMiniBatchesForEnsembleAve = nMBsForEnsembleAve )
          }
        }
    }

    if ( i_th < nIterations ) {
      iterate( i_th + 1, nSteps, nFBsOfBMStates.map( b => b.map(bms => bms.initialMCStates.toVector ).flatten ) )
    }

    } // iterate

    val nFBsOfMCStates = {
       Vector.range(0, nFBs).map( _ => nativeMCStates )
    }
    iterate(1 , nSteps, nFBsOfMCStates )

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

