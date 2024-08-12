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
  import org.sanzo.potts.BM.calcKL
  import org.sanzo.potts.BM.printKL
  import org.sanzo.potts.BM.printLog
  import org.sanzo.potts.BM.miniBatchSizeCorrected
  import org.sanzo.potts.BM.nMiniBatchesInFullBatch
  import org.sanzo.potts.BM.nMiniBatchesForEnsembleAverage
  import org.sanzo.potts.BM.startPosInFB
  import org.sanzo.potts.BM.endPosOfMiniBatch
  import org.sanzo.potts.Util.toProb

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
        
  /*
    val (fia, fijab, stateOrderString) =
    {
                import breeze.linalg.DenseVector
                import breeze.linalg.DenseMatrix
                val nUnits = 10
                val nStatsOfUnit = 21
                val fia = (new Array[DenseVector[Double]](nUnits)).map(
                        pi => DenseVector.fill(nStatsOfUnit){1.0/nStatsOfUnit} )

                val fijab = ArraySeq.range(0, (nUnits * (nUnits - 1)) / 2).map{ ij =>
                                val (i, j) = MCMC.inversePairIndex(ij)
                                val (ij_ , a, b) = MCMC.pairIndex(i, j, 0, 1)
                                assert( ij == ij_ )
                                if ( a == 0 )
                                    fia(i) * fia(j).t
                                else
                                    fia(j) * fia(i).t
                            }
                  } )
                } )

                val stateOrderString = "ARNDCQEGHILKMFPSTWYV-"
                (fia, fijab, stateOrderString)
    }
  */
    val nUnits = fia.size
    val nStatesOfUnit = fia(0).size

    //val (initialConfIDs: ArraySeq[String], initialConfigurations: ArraySeq[ArraySeq[Byte]]) = {
    val nativeConfigs = {
        val fastaFile = new java.io.File( args(1) )
        val fastaIn = isp.getInputStream(fastaFile)
        val fastaLines = scala.io.Source.fromInputStream(fastaIn).getLines()

        val (idArraySeq, seqArraySeq, numArraySeq) =
            MCMC.readFasta(fastaLines, stateOrderString = stateOrderString, caseSensitive = false)

        fastaIn.close()

        assert(idArraySeq.size == seqArraySeq.size && idArraySeq.size == numArraySeq.size )

        val initialConfigurations = numArraySeq

        val initialConfIDs = idArraySeq.map( id => id.split("[ \t]")(0) )

        ArraySeq.range(0, initialConfigurations.size).map{ i => MCMC.NativeConfig( index = i , id = initialConfIDs(i), sample = initialConfigurations(i) ) }
    }

    val fileOfhJ = args(2)
    val initialInteractions: MCMC.Interactions = {

        val file2 = new java.io.File( fileOfhJ )
        val fin2 = isp.getInputStream(file2)

        val bufSource2 = scala.io.Source.fromInputStream(fin2)
        val lines2 = bufSource2.getLines()

      //val (hia, jijab, stateOrderString2) = MCMC.readhJ(lines2)
        val (interactions, stateOrderString2) = MCMC.readInteractions(lines2)
        require( stateOrderString == stateOrderString2)

        bufSource2.close()
      //MCMC.Interactions(hia, jijab)
        interactions
    }

    val pseudoNCounts = 10.0

    val nSamples = 1

  //val initialConfigurations = nativeConfigs.map{ native => native.sample }
  //val nIndependentMC = initialMCStates.size

    val nInitialIterationsPerUnit = 0
    val everyNIterationsPerUnit = 100
    val maxExtendedIterations = 1

    val fullBatchSize = nativeConfigs.size
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

    val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) * nMBsInFB  + 0.5 ).floor.toInt
  //val nMiniBatchesForEnsembleAve = ( math.max( 1.0, ( 10000.0 / (nSamples * miniBSize.toDouble) ) / nMBsInFB ) + 0.5 ).floor.toInt * nMBsInFB
    val nMBsForEnsembleAve = BM.nMiniBatchesForEnsembleAverage(miniBSize, math.max(fullBatchSize, 50000), nMiniBatchesForEnsembleAve )

    val nFBs = (nMBsForEnsembleAve.toDouble / nMBsInFB).ceil.toInt

    val stepsPerEpoch = nMBsInFB.toDouble

    val logInterval = ( (math.max(100.0, stepsPerEpoch) / stepsPerEpoch + 0.5).floor * stepsPerEpoch).toInt
    val logInterval_rev = if (stepsPerEpoch <= 1.0) {
                                 if ( logInterval > 0 ) logInterval else 1
                              } else if ( stepsPerEpoch <= logInterval ) {
                                 (stepsPerEpoch * (logInterval / stepsPerEpoch).floor ).toInt
                              } else {
                                 stepsPerEpoch.toInt
                              }

    val nativeMCStates: Vector[MCMC.State] = {
        val configs = nativeConfigs.map{ n => MCMC.State( nativeConfig = n,
                                                        configuration = n.sample, interactions = initialInteractions ) }
        Vector(configs*)
    }

  /*
    val (interactions, bmState) = BM.initialize(fia, initialInteractions, 
                  // The following arguments are dummy.
                  optMethod = "Adam",
                  learningRate = new LearningRate(), learningRateForRPROPLR = new LearningRateForRPROPLR(),
                  betaV = 0.0, betaM = 0.0 )
                  //
    assert(interactions == initialInteractions)
  */
    val bmState = BMState()

    val (fiaWithPseudo, fijabWithPseudo) = BM.bayesianCorrection(fia, fijab, effectiveNumberOfSamples, pseudoNCounts)

    val nFBsOfBMStatesIntInt: Vector[ArraySeq[Tuple3[BM.BMState, Int, Int]]] = {
        val delta = nFBs * nMBsInFB  - nMBsForEnsembleAve
        Vector.range(0 , nFBs).map{ iFB =>
            val mcStates = Random.shuffle(nativeMCStates)
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
                     initialInteractions, proposedPia = fiaWithPseudo ) 

                  val (revInitialStates: ParVector[MCMC.State], nNonEquil, nExtendedIterations, independentSamplings: ParVector[ArraySeq[MCMC.State]]) =
                    mcmc.runMC(
                      initialStates = initialMCStates,
                      nInitialIterationsPerUnit = nInitialIterationsPerUnit,
                      everyNIterationsPerUnit = everyNIterationsPerUnit,
                      nSamples = nSamples,
                      initialT = 1.0,
                      annealingRate = 0.99,
                      finalT = 1.0,
                      maxExtendedIterations = maxExtendedIterations,
                      kernel = mcmc.kernelGibbsWithMHStep     //mcmc.kernelGibbsWithMHStep==mcmc.kernelMultiBlockMH
                      )

                  val pia = mcmc.frequenciesAtUnitInSamples( independentSamplings )
                  val pijab = mcmc.pairwiseFrequenciesInSamples( independentSamplings )
                  val ensembleAverages = MCMC.EnsembleAverages( pia = pia , pijab = pijab )

                  val bmStateWithMCSamples = bmState.copy( step = step,
                          initialMCStates = revInitialStates, independentSamplings = independentSamplings, ensembleAverages = ensembleAverages )

                  (bmStateWithMCSamples, nNonEquil, nExtendedIterations)
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

    assert( nFBs == nFBsOfBMStates.size )
    assert( nMBsInFB == nFBsOfBMStates(0).size )

    val (miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch) = if(  nFBs <= 1 ) {
            (nFBsOfBMStates(0).slice(nMBsInFB - nMBsForEnsembleAve, nMBsInFB).toVector, nMBsOfBMState) 
          } else {
            (nMBsOfBMState.slice( nFBs * nMBsInFB - nMBsForEnsembleAve, nFBs * nMBsInFB ), nFBsOfBMStates(nFBs - 1).toVector )
          }

    val fileTE = new File(ioFiles("outputDir"),s"out_Energy_distribution_of_MC_samples.txt")
       //(f"Energy_distribution_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outTE  = new PrintStream(fileTE)

    val fileKL = new File(ioFiles("outputDir"),s"out_KL_of_MC_samples.txt")
       //      (f"KL_of_each_step_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
    val outKL  = new PrintStream(fileKL)

    outKL.print(s"# miniBatchSize= ${miniBSize}  fullBatchSize= ${fullBatchSize}  nSamples/MC= ${nSamples}  pseudoNCounts= ${pseudoNCounts}\n")
    outKL.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outKL.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${fileOfhJ}\n")
    outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and %d minibatch ensembles with pseudocounts= %g\n".format(
                    nMBsForEnsembleAve, pseudoNCounts ) )

    outTE.print(s"# miniBatchSize= ${miniBSize}  fullBatchSize= ${fullBatchSize}  nSamples/MC= ${nSamples}  pseudoNCounts= ${pseudoNCounts}\n")
    outTE.print(s"# nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}  maxExtendedIterations= ${maxExtendedIterations}\n")
    outTE.print(s"# Interactions (hia, Jijab) are the same for all steps; read from ${fileOfhJ}\n#\n")
    outTE.print("# TE is calculated in the Ising gauge for comparison.\n#\n")
    outTE.print("# For <TE>_m and <(TE -<TE>)^2>_m, %d minibatch ensembles are employed for averaging.\n#\n".format(nMBsForEnsembleAve))

  //outTE.print("# step\tL=%d\t(all randomized sample mean - all randomized sample variance/L)  TE_/L  <TE>/L\tall randomized sample mean\tall randomized sample variance/L  TE sample variance/L  <(TE-<TE>)^2>/L\n".format(nUnits) )
  //outTE.print("# step\tL=%d\t(random sample mean - random sample variance)/L  TE_all/L  <TE>/L  <TE>_m/L\trandom sample mean\trandom sample variance/L  TE all sample variance/ L  <(TE-<TE>)^2>/L  <(TE-<TE>)^2>_m/L ; <...>_m over max %d minibatch ensembles\n".format(nUnits, nMiniBatchesForEnsembleAve) )


    val fa = toProb(fia)

    for ( iMB <- Range(0, nMBsOfBMState.size) ) {
        if ( nMBsOfBMState(iMB).step >= 0 ) {

          val nNonEquil = nMBsOfNnonEquiv(iMB)
          val nExtendedIterations = nMBsOfNExtend (iMB) 

          val ensembleKL = {
              val ensembleAverages = nMBsOfBMState(iMB).ensembleAverages
              val independentS = nMBsOfBMState(iMB).independentSamplings
              val totalNSamples = (independentS.size * independentS(0).size).toDouble
                  BM.calcKL(fiaWithPseudo, fijabWithPseudo, observedNP = Double.PositiveInfinity, pseudoNP = 0.0,
                      ensembleAverages.pia, ensembleAverages.pijab, observedNQ = totalNSamples, pseudoNQ = pseudoNCounts)
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

                        BM.calcKL(fiaWithPseudo, fijabWithPseudo, observedNP = Double.PositiveInfinity, pseudoNP = 0.0,
                          pia, pijab, observedNQ = totalNSamples, pseudoNQ = pseudoNCounts)
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

            val logFiles: Option[Tuple2[java.io.File, java.io.File]] =
                  if ( (bmStateWithMCSamples.step + 1) % logInterval_rev == 0 || bmStateWithMCSamples.step == 0 ) {
                    val logFileInteractions = new File(ioFiles("outputDir"), s"out_Interactions_hJ_for_MCMC_${bmStep}.txt.gz")
                    val logFileIndpendentMCsamplings = new File(ioFiles("outputDir"), s"out_MC_samples_by_MCMC_${bmStep}.fasta.gz")
                    Option( (logFileInteractions, logFileIndpendentMCsamplings) )
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
                ("outputDir", new java.io.File(args(0) ) ) )
      args(1) match {
        case _ => {
                //val fout = scala.sys.process.stdout   //System.out 
                processFiles(args.drop(1), ioFiles)
                }
      }
   }

}

