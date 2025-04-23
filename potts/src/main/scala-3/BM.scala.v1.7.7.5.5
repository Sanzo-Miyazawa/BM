
//import org.{sanzo => miyazawa}                                // for scala 2.13
import org.{sanzo as miyazawa}                                  // for scala 3


package org.sanzo.potts {

// Use scala.util.Random.setSeed(seed: Long), if you want consistent behavior with the same random numbers from execution to execution.

  import scala.util.Random

  //import scala.collection.parallel.CollectionConverters._     // for par in scala 2.13
  import scala.collection.parallel.CollectionConverters.*       // for par in scala 3
  import scala.math.Ordering.Double.TotalOrdering               // for scala 2.13

  import scala.reflect.ClassTag

  import breeze.linalg
  import breeze.linalg.DenseMatrix
  import breeze.linalg.DenseVector
//import breeze.linalg.`*`
//import breeze.linalg.{max, sum}
//import breeze.numerics.{exp, log}

  import breeze.stats.meanAndVariance
  import breeze.stats.meanAndVariance.MeanAndVariance
//import breeze.stats.distributions.Uniform
  import breeze.stats.distributions.Process
  import breeze.stats.distributions.Rand
  import breeze.stats.distributions.RandBasis

/* This part is moved at functions where the default randBasis is needed; MCMC.randomConfiguration , bm.initializehJ
 
// As of Breeze 2, import this if you want "matlab"-like behavior with different random numbers from execution to execution
//import breeze.stats.distributions.Rand.VariableSeed.*

// As of Breeze 2, import this if you want consistent behavior with the same random numbers from execution to execution (modulo threading or other sources of nondeterminacy)
  //import breeze.stats.distributions.Rand.FixedSeed._          //for scala 2.13
  import breeze.stats.distributions.Rand.FixedSeed.*            //for scala 3
*/

  //import spire.math.Integral
  import spire.syntax.cfor.cfor                                 //for scala 2.13
  //import spire.syntax.fastaFor.{fastFor as cfor}              //for scala 3

  import java.io.File 
  import java.io.{FileInputStream, FileOutputStream}
  import java.io.PrintStream
  import java.io.PrintWriter
  import java.util.zip.{GZIPInputStream, GZIPOutputStream}

  import scala.sys.process.stdin
  import scala.sys.process.stdout
  import scala.sys.process.stderr
  import scala.io.Source
  import scala.annotation.tailrec
  import scala.collection.parallel.mutable.ParArray
  import scala.collection.parallel.immutable.ParVector
  //import scala.collection.mutable.HashMap
  import scala.collection.mutable
  import scala.collection.mutable.ArrayBuffer
  import scala.collection.immutable
  import scala.collection.immutable.ArraySeq
  import scala.collection.immutable.HashMap
  import scala.collection.immutable.Queue
  import scala.collection.immutable.Vector
  import scala.math
  import scala.util.matching.Regex
  import scala.util.Sorting.stableSort

  //import miyazawa.sequence.SequenceArraySeq.pairIndex
  import miyazawa.potts.MCMC.pairIndex
  import miyazawa.potts.MCMC.{inversePairIndex, fromPairIndex}
  import miyazawa.potts.Util.{energyOfRandomSeq, toProb}
  import miyazawa.potts.{LearningRate, LearningRateForRPROPLR} 

  import miyazawa.potts.zip.gzipOutputStream

  object SetSeedForRand {
    import breeze.stats.distributions.Rand
    import breeze.stats.distributions.RandBasis
  //private val defaultSeed = Int.MinValue
  //private var fixedRandBasis = Rand.VariableSeed.randBasis
    private val defaultSeed = -1
    private var fixedRandBasis = Rand.FixedSeed.randBasis
    private var pSeed = defaultSeed
    private var randBasisUsed = false
  //private var nCalled = 0
    def apply(i: Int = 0) = {
      assert( randBasisUsed == false, "SetSeedForRand(seed) cannot be called after randBasis is used.\n")
  //  assert( nCalled <= 0, s"Already a seed for Rand was set; ${pSeed}\n")
  //  nCalled = nCalled + 1 
      pSeed = i
      fixedRandBasis =
        if ( i == defaultSeed ) {
          Rand.FixedSeed.randBasis
        } else if ( i == Int.MinValue ) {
          Rand.VariableSeed.randBasis
        } else {
          assert ( i >= 0, "An argument for SetSeedForRand must be a non-negative integer to be passed to RandBasis.withSeed or equal to ${defaultSeed} to use Rand.FixedSeed.randBasis or Int.MinValue to use Rand.VariableSeed.randBasis." )
          RandBasis.withSeed(i)
        }
    //assert( seed == pSeed, "SetSeedForRand.seed cannot be accessed before SetSeedForRand(i) is called." )
      assert( seed == pSeed, "SetSeedForRand.seed was accessed before SetSeedForRand(i) is called or two different seeds, ${seed} and ${pSeed}, have been set." )
    }
    lazy val seed = pSeed
    lazy given randBasis: RandBasis = {
               randBasisUsed = true
               fixedRandBasis
             }
  }

// Call directly scala.util.Random.setSeed(seed) in the main.
//def setSeedForRandom(seed: Long): Unit = {
//  Random.setSeed(seed)
//}

// A new seed for each of markovChain that is executed in parallel to reproduce the same result.
  def setRandBasisIArrayWithRandomSeeds( n: Int ): IArray[RandBasis] = {
    require( n > 0 )
    val randBasisIArray = IArray.range(0,n)
    randBasisIArray.map{ _ =>
      val seed = scala.util.Random.nextInt()
      RandBasis.withSeed(seed)
    }
  }


  object MCMC {

  /*
    val randBasisIArray: IArray[RandBasis] = 
  */

    type Configuration        = IArray[Byte]
    type IndependentMC        = ParVector[Iterator[MCMC.State]]       //type   IndependentMC[T] = ParVector[T]
    type IndependentSamplings = ParVector[ArraySeq[MCMC.State]]       //type   IndependentSamplings[T] = ParVector[T]
    type InitialMCStates      = ParVector[MCMC.State]

    case class NativeConfig(val index: Int = 0, val id: String = "", val sample: Configuration = IArray[Byte]() )

    //{
    //   require ( kT >= finalkT )
    //}
    case class Interactions(val hia: ArraySeq[DenseVector[Double]] = ArraySeq[DenseVector[Double]](),
                            val jijab: ArraySeq[DenseMatrix[Double]] = ArraySeq[DenseMatrix[Double]](),
                            val gauge: String = "unused")

    case class State(val nativeConfig: NativeConfig = NativeConfig(),
                 val configuration: Configuration, val configurationArray: Array[Byte] = Array[Byte](),
                 val interactions: Interactions = Interactions(),
                 val energy: Double = 0.0, val kT: Double = 0.0, val annealingRate: Double = 0.0,
                 val step: Int = 0, val finalkT: Double = 0.0)

    case class ProposedDistributions(val proposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]],
        val logProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ,
        val log1_ProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] )

    case class EnsembleAverages(val pia: ArraySeq[DenseVector[Double]] = ArraySeq[DenseVector[Double]](),
                                val pijab: ArraySeq[DenseMatrix[Double]] = ArraySeq[DenseMatrix[Double]]() )


    def pairIndexL(i: Long, j: Long): Long = {
                require( i >= 0 && j >= 0 )
                if(i > j) {
                        val ij = (i * (i - 1) / 2 + j)
                        assert( ij >= 0, "pairIndex is too large for Long.\n" )
                        ij
                } else if ( j > i ) {
                        val ji = (j * (j - 1) / 2 + i)
                        assert( ji >= 0, "pairIndex is too large for Long.\n" )
                        ji
                } else {
                        sys.error("*** Error: pairIndex must be for i != j\n")
                        -1L
                }
    }

    def pairIndexL(i: Long, j: Long, ai: Int, aj: Int): Tuple3[Long, Int, Int] = {
                require( i >= 0 && j >= 0 )
                val pIndex = pairIndexL(i, j)
                if(i > j) {
                        val ij = pIndex
                        (ij, ai, aj)
                } else { // j > i
                        assert( i != j , "*** Error: pairIndex must be for i != j\n")
                        val ji = pIndex
                        (ji, aj, ai)
                }
    }

  //def pairIndex(i: Long, j: Long): Int = {
    def pairIndex(i: Int, j: Int): Int = {
                require( i >= 0 && j >= 0 )
                val pIndex = pairIndexL(i, j).toInt
                assert( pIndex >= 0, "pairIndex is too large for Int.\n" )
                pIndex
    }

  //def pairIndex(i: Long, j: Long, ai: Int, aj: Int): Tuple3[Int, Int, Int] = {
    def pairIndex(i: Int, j: Int, ai: Int, aj: Int): Tuple3[Int, Int, Int] = {
                require( i >= 0 && j >= 0 )
                val (pIndexL, a1, a2)  = pairIndexL(i, j, ai, aj)
                val pIndex = pIndexL.toInt
                assert( pIndex >= 0, "pairIndex is too large for Int.\n" )
                (pIndex, a1, a2)
    }

  //def inversePairIndex(ij: Int) = {
    def inversePairIndex(ij: Long): Tuple2[Int, Int] = {
                require( ij >= 0 )
                val i = {
                // (i-1) <= sqrt(i*(i-1)) <= sqrt(2*ij)=sqrt(i*(i-1)+ 2*j) < sqrt(i*(i+1)) < (i+1)
                // (i-1) <= floor ( sqrt(2*ij ) < (i+1)
                // floor ( sqrt(2*ij ) = (i-1) or i
                //
                          val i = math.sqrt(ij * 2.0).floor.toLong
                          if ( ij - (i * ( i - 1)) / 2 < i )
                                i
                          else
                                i + 1
                        }
                val j = ij - (i * ( i - 1)) / 2  
                require( pairIndexL(i, j) == ij )
                val iInt = i.toInt
                val jInt = j.toInt
                assert( iInt >= 0 && jInt >= 0, "pairIndex is too large for Int pairs." )
                
                (iInt, jInt)
    }

  //def fromPairIndex(ij: Int) = inversePairIndex(ij)
    def fromPairIndex(ij: Long) = inversePairIndex(ij)


    def readNSamples(lines: collection.Iterator[String] ) = {
        val (lines0, lines1) = lines.duplicate

        val nHeaderLine = ("""^#[ \t]*Effective_number_of_samples:[ \t]*([0-9.eE+]+).*$""").r
        val lines11 = lines1.dropWhile(line => nHeaderLine.findFirstMatchIn(line).isEmpty )
        val effectiveNumberOfSamples =
          if ( lines11.hasNext ) {
            val line = lines11.next()
            //val nHeaderLine(effectiveNumberOfSamplesString) = line : @unchecked
            val effectiveNumberOfSamples = nHeaderLine.findFirstMatchIn(line).get.group(1).toDouble
            effectiveNumberOfSamples
          } else {
            0.0
          }

        (effectiveNumberOfSamples, lines0)
    }

    def readGauge(lines: collection.Iterator[String]) = {
        val (lines0, lines1) = lines.duplicate
        val nHeaderLine = ("""^#[ \t]*Gauge:[ \t]*([^ \t,:;]*).*$""").r
        val lines11 = lines1.dropWhile(line => nHeaderLine.findFirstMatchIn(line).isEmpty )
        val gauge =
          if ( lines11.hasNext ) {
            val line = lines11.next()
            //val nHeaderLine(gaugeString) = line : @unchecked
            val gauge = nHeaderLine.findFirstMatchIn(line).get.group(1)
            gauge
          } else {
            "unused"
          }

        (gauge, lines0)
    }

    def readN(lines: collection.Iterator[String]) = {
        val (lines0, lines1) = lines.duplicate
        val nHeaderLine = ("""^#[ \t]*n_units: *([0-9]+)[ \t]*n_states_of_unit:[ \t]*([0-9]+).*$""").r
        val lines11 = lines1.dropWhile(line => nHeaderLine.findFirstMatchIn(line).isEmpty )
        val (nUnits, nStatesOfUnit, stateOrderString) =
          if ( lines11.hasNext ) {
            val line = lines11.next()
            //val nHeaderLine(nUnitsString, nStatesOfUnitString) = line : @unchecked
            //val (nUnits, nStatesOfUnit) = (nUnitsString.toInt, nStatesOfUnitString.toInt ) 
            //val matchN = nHeaderLine.findFirstMatchIn(line) 
            //val (nUnits, nStatesOfUnit) = ( matchN.get.group(1).toInt, matchN.get.group(2).toInt )
            val (nUnits, nStatesOfUnit) = line match {
                                                case nHeaderLine(nUnitsStr, nStatesOfUnitStr) => ( nUnitsStr.toInt, nStatesOfUnitStr.toInt )
                                          }

            val sHeaderLine = ("""^.*[ \t]+State_Order_String:[ \t]*([^ \t]+).*$""").r
            /*
            val matchS = sHeaderLine.findFirstMatchIn(line)
            val stateOrderString: String = matchS match {
                          case Some(x) => x.group(1) 
                          case None => ""       // null
                        }
            */
            val stateOrderString: String = line match {
                          case sHeaderLine(stateOrderStr) => stateOrderStr
                          case _ => ""          // null
                        }
            (nUnits, nStatesOfUnit, stateOrderString)
          } else {
            (0, 0, "")  // (0, 0, null)
          }

        (nUnits, nStatesOfUnit, stateOrderString, lines0)
    }

    def readhOrPi(lines: collection.Iterator[String], nUnits: Int, nStatesOfUnit: Int, whichhOrPi: String ) = {
        //Error: val lines1 = lines     //lines can not be used later because this is not a duplicate and lines1.dropWhile is done.
        val (lines0, lines1) = lines.duplicate
        
        val commentline = """^#.*""".r
        val hPiHeaderLine = ("""^""" + whichhOrPi + """[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9eE.+-]+).*$""").r
        val hPiLine = hPiHeaderLine
        val lines11 = lines1.dropWhile(line => hPiLine.findFirstMatchIn(line).isEmpty )

      //val h = ArraySeq.fill(nUnits)(0).map( _ => new DenseVector[Double](nStatesOfUnit) )
        val h = ArraySeq.fill(nUnits)(0).map( _ => DenseVector.fill(nStatesOfUnit)(0.0) )

        @annotation.tailrec
        def processhPiLine(lines: collection.Iterator[String], h: ArraySeq[DenseVector[Double]], n: Int): Tuple2[Int, collection.Iterator[String]] = {
                  if(lines.hasNext) {
                    val line = lines.next()
                    if ( (hPiLine.findFirstMatchIn(line) ).nonEmpty ) {
                        val field = line.split("""[ \t]+""")
                        val ir = field(1).toInt + 1
                        val ia_orig = field(2).toInt

                        if ( ir <= nUnits ) {
                        //cfor(0)(i => i < nAADelTypes, i => i + 1) ( i => {
                        //val i = ia_orig
                        //      val ia = aaDelOrderMap.compoundOrderNumber( aaOrderString(i), -1 )      /* error generator */
                        //      h(ir - 1)(ia) = field(i + 1).toDouble
                          val ia = ia_orig
                          h(ir - 1)(ia) = field(3).toDouble
                        //} )
                          if( n + 1 >= nUnits * nStatesOfUnit) { 
                                (n + 1, lines)
                          } else {
                                processhPiLine(lines, h, n + 1)
                          }
                        } else {
                          sys.error("*** Error: the number of residues is larger than %d\n".format(nUnits))
                          (- n, lines)
                        }
                    } else if ( (commentline.findFirstMatchIn(line) ).nonEmpty ) {
                        processhPiLine(lines, h, n)
                    } else {
                        if ( n > 0 )
                          sys.error("*** Error in %s lines: %s".format(whichhOrPi, line))
                        //(- n, lines)
                        processhPiLine(lines, h, n)
                    }
                  } else {
                        sys.error("*** Error: no enough number of %s lines.\n".format(whichhOrPi))
                        (- n, lines)
                  }
        }
        val (nh, lines12) = processhPiLine(lines11, h, 0)
        if( nh != nUnits * nStatesOfUnit )
                sys.error("*** Error: no enough number of %s lines: %d\n".format(whichhOrPi, nh))

        (h, lines12)
        //Or (h, lines0)
    }

    def readJOrPij(lines: collection.Iterator[String], nUnits: Int, nStatesOfUnit: Int, whichJOrPij: String ) = {
        //Error: val lines1 = lines     //lines can not be used later because this is not a duplicate and lines1.dropWhile is done.
        val (lines0, lines1) = lines.duplicate
        
        val commentline = """^#.*""".r
        val JPijHeaderLine = ("""^""" + whichJOrPij + 
                """[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9eE.+-]+).*$""").r
        val JPijLine = JPijHeaderLine
        val lines11 = lines1.dropWhile(line => JPijLine.findFirstMatchIn(line).isEmpty )

        val nPairs = (nUnits * (nUnits - 1)) / 2
      //val J = ArraySeq.fill(nPairs)(0).map( _ => new DenseMatrix[Double](nStatesOfUnit, nStatesOfUnit) )
        val J = ArraySeq.fill(nPairs)(0).map( _ => DenseMatrix.fill(nStatesOfUnit, nStatesOfUnit)(0.0) )

        @annotation.tailrec
        def processJPijLine(lines: collection.Iterator[String], J: ArraySeq[DenseMatrix[Double]], n: Int): Tuple2[Int, collection.Iterator[String]] = {
                  if(lines.hasNext) {
                    val line = lines.next()
                    if ( (JPijLine.findFirstMatchIn(line) ).nonEmpty ) {
                        val field = line.split("""[ \t]+""")
                        val ir = field(1).toInt + 1
                        val jr = field(2).toInt + 1
                        val ia_orig = field(3).toInt
                      //val aa = aaOrderString(ia_orig)
                      //val ia = aaDelOrderMap.compoundOrderNumber( aa, -1 )
                        val ia = ia_orig

                        val ja_orig = field(4).toInt
                        val ja = ja_orig

                        if ( ir <= nUnits && jr <= nUnits) {
                        //cfor(0)(j => j < nAADelTypes, j => j + 1) ( j => {
                        //val j = ja_orig
                        //      val ia = aaDelOrderMap.compoundOrderNumber( aaOrderString(i), -1 )      /* error generator */
                        val (k, a, b) = pairIndex(jr -1, ir -1, ja, ia)
                        //      J(k)(a) = field(i + 3).toDouble
                        J(k)(a, b) = field(5).toDouble
                        //} )
                          val nn = (nUnits * ( nUnits - 1)) / 2 * nStatesOfUnit * nStatesOfUnit
                          if( n + 1 >= nn ) { 
                                //(ir == nUnits && jr == nUnits -1) ||
                                //(ir == nUnits - 1 && jr == nUnits) ) {
                                (n + 1, lines)
                          } else {
                                processJPijLine(lines, J, n + 1)
                          }
                        } else {
                          sys.error("*** Error: the number of residues is larger than %d\n".format(nUnits))
                          (- n, lines)
                        }
                    } else if ( (commentline.findFirstMatchIn(line) ).nonEmpty ) {
                        processJPijLine(lines, J, n)
                    } else {
                        if ( n > 0 )
                          sys.error("*** Error in %s lines: %s".format(whichJOrPij, line))
                        //(- n, lines)
                        processJPijLine(lines, J, n)
                    }
                  } else {
                        sys.error("*** Error: no enough number of %s lines.\n".format(whichJOrPij))
                        (- n, lines)
                  }
        }
        val (nJ, lines12 )= processJPijLine(lines11, J, 0)
        if( nJ != nUnits * (nUnits - 1) / 2 * nStatesOfUnit * nStatesOfUnit )
                sys.error("*** Error: no enough number of %s lines: %d\n".format(whichJOrPij, nJ))
        
        (J, lines12)
        //Or (J, lines0)
    }

    def readhJ(lines: collection.Iterator[String] ) = {

        val (nUnits, nStatesOfUnit, stateOrderString, lines1) = readN(lines)
        val (gauge, lines2) = readGauge(lines1)
        val ( h, lines3 ) = readhOrPi(lines2, nUnits, nStatesOfUnit, "h" )
        val ( j, lines4 ) = readJOrPij(lines3, nUnits, nStatesOfUnit, "J" )

        (h, j, stateOrderString)
    }

    def readInteractions(lines: collection.Iterator[String] ) = {

        val (nUnits, nStatesOfUnit, stateOrderString, lines1) = readN(lines)
        val (gauge, lines2) = readGauge(lines1)
        val ( h, lines3 ) = readhOrPi(lines2, nUnits, nStatesOfUnit, "h" )
        val ( j, lines4 ) = readJOrPij(lines3, nUnits, nStatesOfUnit, "J" )

        (Interactions(h, j, gauge), stateOrderString)
    }

    def readPiPij(lines: collection.Iterator[String] ) = {

        val (effectiveNumberOfSamples, lines0) = readNSamples(lines)
        val (nUnits, nStatesOfUnit, stateOrderString, lines1) = readN(lines0)
        val ( pia, lines2 ) = readhOrPi(lines1, nUnits, nStatesOfUnit, "Pi" )
        val ( pijab, lines3 ) = readJOrPij(lines2, nUnits, nStatesOfUnit, "Pij" )

        (effectiveNumberOfSamples, pia, pijab, stateOrderString)
    }

    def printNSamples(out: PrintStream, effectiveNumberOfSamples: Double) = {
        out.print("# Effective_number_of_samples: %g\n".format( effectiveNumberOfSamples ) )
    }

    def printhOrPi(outhOrPi: PrintStream, stateOrderString: String, 
                hOrPi: ArraySeq[DenseVector[Double]], whichhOrPi: String ): Unit = {

        val nUnits = hOrPi.size
        val nStatesOfUnit = hOrPi(0).size

        outhOrPi.print("# n_units: %d  n_states_of_unit: %d".format( nUnits, nStatesOfUnit ) )

        if( stateOrderString != null && stateOrderString.size != 0 )
          outhOrPi.print("  State_Order_String: %s\n".format(stateOrderString) )
        else
          outhOrPi.print("\n")

        cfor(0)(i => i < nUnits, i => i + 1)( i => {
                cfor(0)(ia => ia < nStatesOfUnit, ia => ia + 1)( ia => {
                  //if( stateOrderString != null && stateOrderString.size != 0 ) {
                  //  outhOrPi.print("h %d %c %g\n".format(i,
                  //    stateOrderString(ia), hOrPi(i)(ia) ) )
                  //} else
                  {
                    outhOrPi.print("%s %d %d %g\n".format(whichhOrPi, i,
                        ia, hOrPi(i)(ia) ) )
                  }
                } )
        } ) 
    }

  /*
    def printhOrPi(outhOrPi: PrintStream, stateOrderString: String, 
                hOrPi: ArraySeq[ArraySeq[Double]], whichhOrPi: String ): Unit = {

        printhOrPi(outhOrPi, stateOrderString,
                hOrPi.map( x => DenseVector(x.toArray) ), whichhOrPi)
     }
   */

    def printJOrPij(outJOrPij: PrintStream, stateOrderString: String, 
                JOrPij: ArraySeq[DenseMatrix[Double]], whichJOrPij: String ) = {

        val nPairs = JOrPij.size
        val nUnits = {
                val n = math.sqrt(nPairs * 2.0).toInt
                val np = (n * (n - 1)) / 2
                val nu = if ( np == nPairs ) n else n + 1
                if ( (nu * ( nu - 1)) / 2 != nPairs ) {
                        sys.error("algorithm error")
                        nu
                } else
                        nu
            }
        
        val nStatesOfUnit = JOrPij(0).rows

        outJOrPij.print("# n_units: %d  n_states_of_unit: %d".format( nUnits, nStatesOfUnit ) )

        if( stateOrderString != null && stateOrderString.size != 0 )
          outJOrPij.print("  State_Order_String: %s\n".format(stateOrderString) )
        else
          outJOrPij.print("\n")

        cfor(0)(i => i < nUnits, i => i + 1)( i => {
          cfor(i + 1)(j => j < nUnits, j => j + 1)( j => {
                cfor(0)(ia => ia < nStatesOfUnit, ia => ia + 1)( ia => {
                cfor(0)(jb => jb < nStatesOfUnit, jb => jb + 1)( jb => {
                  val (ij, a, b) = pairIndex(i, j, ia, jb) 
                  //if( stateOrderString != null && stateOrderString.size != 0 ) {
                  //  outJOrPij.print("J %d %d %c %c %g\n".format(i, j, 
                  //    stateOrderString(ia), stateOrderString(jb), interactions.jijab(ij)(a, b) ) )
                  //} else {
                    outJOrPij.print("%s %d %d %d %d %g\n".format(whichJOrPij, i, j, 
                        ia, jb, JOrPij(ij)(a, b)   ))
                  //}
                } )
                } )
          } ) 
        } ) 
    }

  /*
    def printJOrPij(outJOrPij: PrintStream, stateOrderString: String, 
                JOrPij: ArraySeq[ArraySeq[ArraySeq[Double]]], whichJOrPij: String ): Unit = {
        import breeze.linalg.DenseMatrix

        printJOrPij(outJOrPij, stateOrderString, 
                JOrPij.map(pij => {
                        //DenseMatrix( pij: _*)         // for 2.11
                          DenseMatrix( pij* )           // for 2.13
                        //new DenseMatrix( pij(0).size, pij.size, pij.flatten.toArray).t        // for 2.13
                } ), 
                whichJOrPij )
    }
  */

    def printInteractions(outhJ: PrintStream, stateOrderString: String, 
                interactions: Interactions ) = {

        val h = interactions.hia
        val J = interactions.jijab

        outhJ.println("# Gauge: " + interactions.gauge )
        printhOrPi(outhJ, stateOrderString, h, "h" )
        printJOrPij(outhJ, stateOrderString, J, "J")
    }

    def printhJ(outhJ: PrintStream, stateOrderString: String, 
                hia: ArraySeq[DenseVector[Double]],
                jijab: ArraySeq[DenseMatrix[Double]] ) = {

        printhOrPi(outhJ, stateOrderString, hia, "h" )
        printJOrPij(outhJ, stateOrderString, jijab, "J")
    }

  /*
    def printhJ(outhJ: PrintStream, stateOrderString: String, 
                hia: ArraySeq[ArraySeq[Double]],
                jijab: ArraySeq[ArraySeq[ArraySeq[Double]]] ) = {

        printhOrPi(outhJ, stateOrderString, hia, "h" )
        printJOrPij(outhJ, stateOrderString, jijab, "J")
    }
  */

    def printPiPij(outPiPij: PrintStream,
                stateOrderString: String, 
                pia: ArraySeq[DenseVector[Double]],
                pijab: ArraySeq[DenseMatrix[Double]] ) = {

        printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
        printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }

    def printPiPij(outPiPij: PrintStream,
                effectiveNumberOfSamples: Double, stateOrderString: String, 
                pia: ArraySeq[DenseVector[Double]],
                pijab: ArraySeq[DenseMatrix[Double]] ) = {

        printNSamples(outPiPij, effectiveNumberOfSamples)
        printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
        printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }

  /*
    def printPiPij(outPiPij: PrintStream,
                stateOrderString: String, 
                pia: ArraySeq[ArraySeq[Double]],
                pijab: ArraySeq[ArraySeq[ArraySeq[Double]]] ) = {

        printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
        printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }

    def printPiPij(outPiPij: PrintStream,
                effectiveNumberOfSamples: Double, stateOrderString: String, 
                pia: ArraySeq[ArraySeq[Double]],
                pijab: ArraySeq[ArraySeq[ArraySeq[Double]]] ) = {

        printNSamples(outPiPij, effectiveNumberOfSamples)
        printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
        printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }
  */

    def printIndependentMCsamplings(outSamples: PrintStream, stateOrderString: String, 
          independentSamplings: immutable.IndexedSeq[ArraySeq[State]] ) : Unit = {
        //sampleIndex: immutable.IndexedSeq[Int], sampleIDs: immutable.IndexedSeq[String] ) : Unit = {

        val independentSamplingsArraySeq = independentSamplings

        val nSamplings = independentSamplingsArraySeq.size
      //val nSamples = independentSamplingsArraySeq(0).size
        val nUnits = independentSamplingsArraySeq(0)(0).configuration.size

        if( stateOrderString == null || stateOrderString.size == 0) {
          //outSamples.print("# State_Order_String: %s\n".format(stateOrderString) )
                cfor(0)(i => i < nSamplings, i => i + 1)( i => {
                        //outSamples.print("# Independent_Sampling: %d / %d\n".format(i, nSamplings))
                        val nSamples = independentSamplingsArraySeq(i).size
                        cfor(0)(j => j < nSamples, j => j + 1)( j => {
                                val state = independentSamplingsArraySeq(i)(j)
                                outSamples.print("> %s %d/%d_%d/%d  sampling=%d/%d sample=%d/%d step=%d gauge=%s E=%g kT=%g annealing_rate=%g\n".format(
                                        state.nativeConfig.id,
                                        state.nativeConfig.index, nSamplings, j, nSamples,
                                        state.nativeConfig.index, nSamplings, j, nSamples,
                                        state.step, state.interactions.gauge, state.energy, state.kT, state.annealingRate) )
                                cfor(0)(k => k < nUnits, k => k + 1)( k => {
                                        outSamples.print("%d ".format( state.configuration(k) ) ) } )
                                outSamples.print("\n")
                        } )
                } )
          
        } else {
                cfor(0)(i => i < nSamplings, i => i + 1)( i => {
                        val nSamples = independentSamplingsArraySeq(i).size
                        cfor(0)(j => j < nSamples, j => j + 1)( j => {
                                val state = independentSamplingsArraySeq(i)(j)
                                outSamples.print("> %s %d/%d_%d/%d  sampling=%d/%d sample=%d/%d step=%d gauge=%s E=%g kT=%g annealing_rate=%g\n".format(
                                        state.nativeConfig.id,
                                        state.nativeConfig.index, nSamplings, j, nSamples,
                                        state.nativeConfig.index, nSamplings, j, nSamples,
                                        state.step, state.interactions.gauge, state.energy, state.kT, state.annealingRate) )
                                outSamples.print("%s\n".format( state.configuration.map(s 
                                        => stateOrderString(s).toString).reduceLeft(_ + _) ) )
                        } )
                } )
        }
        
    }

    def printIndependentMCsamplings(outSamples: PrintStream, stateOrderString: String, 
                independentSamplings: IndependentSamplings ): Unit = {
              //sampleIndex: Vector[Int], sampleIDs: ArraySeq[String] ): Unit = {

        printIndependentMCsamplings(outSamples, stateOrderString,
                independentSamplings.toVector )
              //independentSamplings.toVector, sampleIndex, sampleIDs )
    }

/*
    def printIndependentMCsamplings(outSamples: PrintStream, stateOrderString: String, 
                independentSamplingsVector: IndependentSamplings,
                sampleIndex: Vector[Int], sampleIDs: ArraySeq[String] ): Unit = {

        printIndependentMCsamplings(outSamples, stateOrderString,
                independentSamplings.toVector, sampleIndex, sampleIDs )
    }
*/

    def proposeState(prob: Double, currentState: Int, proposedDistributionAtUnit: DenseVector[Double]): Byte = {
        val nStatesofUnit =  proposedDistributionAtUnit.size

        //val prob = Uniform(0.0, 1.0 - proposedDistributionsAtAllUnit(currentS))
        @annotation.tailrec
        def searchState(current: Int, i: Int, prob: Double): Byte  = {
          if( i == current )
                if ( i == 0 )
                        1.toByte
                else
                        searchState(current, i - 1, prob) 
          else if( i <= 0 )
                0.toByte
          else if( prob > proposedDistributionAtUnit(i) )
                searchState(current, i - 1, prob - proposedDistributionAtUnit(i))
          else 
                i.toByte
        }
        searchState(currentState, nStatesofUnit - 1, prob ) 
    }

    def proposeState(currentState: Int, proposedDistributionAtUnit: DenseVector[Double])(using rand: RandBasis): Byte = {

        import breeze.stats.distributions.Uniform

        // x <= Uniform(x, y) < y ; see breeze.stats.distributions.{Uniform, Rand.uniform}
        // Just in case, it is better to code   u = Uniform(x, y) if u != y
        def uniformGTxLEy(x: Double, y: Double): Rand[Double] = 
          for {
                r <- Uniform(x, y)      if r < y
          } yield { y + x - r }
        val nStatesofUnit =  proposedDistributionAtUnit.size
        def uniformP(): Double = {
          val (prob, maxP) = if ( currentState < 0 || currentState >= nStatesofUnit ) {
                          (uniformGTxLEy(0.0, 1.0).sample(), 1.0)
                        } else {
                          val maxP = math.max(1.0 - proposedDistributionAtUnit(currentState), 0.0)
                          (uniformGTxLEy(0.0, maxP).sample(), maxP)
                        }
          prob
        }
        proposeState(uniformP(), currentState, proposedDistributionAtUnit)
    }

  //def totalE(configuration: mutable.IndexedSeq[Byte],
    def totalE(configuration: IArray[Byte],
        interactions: Interactions ): Double = {

        val nUnits = configuration.size
        val hia = interactions.hia
        val Jijab = interactions.jijab

        /*
        val ne = new Array[Double](nUnits)
        val nei = new Array[Double](nUnits)

        cfor(0)(i => i < nUnits, i => i + 1)( i => {
          val ia = configuration(i)
          cfor(0)(j => j < nUnits, j => j + 1)( j => {
            if ( i != j ) { 
              val jb = configuration(j)
              val (ij, a, b) = pairIndex(i, j, ia, jb)
              nei(j) = Jijab(ij)(a, b) * 0.5
            } else {
              nei(i) = hia(i)(ia)
            }
          } )
          ne(i) = nei.sum 
        } )
        (- ne.sum)
        */

        val ne = IArray.range(0, nUnits).map ( i => {
                    val ia = configuration(i)
                    val nei = IArray.range(0, nUnits).map ( j => {
                                if ( i != j ) { 
                                        val jb = configuration(j)
                                        val (ij, a, b) = pairIndex(i, j, ia, jb)
                                        Jijab(ij)(a, b) * 0.5
                                } else {
                                        hia(i)(ia)
                                }
                              } )
                    nei.sum
                 } )
        (- ne.sum)
    }

    def totalE(configuration: IndexedSeq[Byte],
        interactions: Interactions ): Double = {

        val nUnits = configuration.size
        val hia = interactions.hia
        val Jijab = interactions.jijab

        /*
        val ne = new Array[Double](nUnits)
        val nei = new Array[Double](nUnits)

        cfor(0)(i => i < nUnits, i => i + 1)( i => {
          val ia = configuration(i)
          cfor(0)(j => j < nUnits, j => j + 1)( j => {
            if ( i != j ) { 
              val jb = configuration(j)
              val (ij, a, b) = pairIndex(i, j, ia, jb)
              nei(j) = Jijab(ij)(a, b) * 0.5
            } else {
              nei(i) = hia(i)(ia)
            }
          } )
          ne(i) = nei.sum 
        } )
        (- ne.sum)
        */

        val ne = IArray.range(0, nUnits).map ( i => {
                    val ia = configuration(i)
                    val nei = IArray.range(0, nUnits).map ( j => {
                                if ( i != j ) { 
                                        val jb = configuration(j)
                                        val (ij, a, b) = pairIndex(i, j, ia, jb)
                                        Jijab(ij)(a, b) * 0.5
                                } else {
                                        hia(i)(ia)
                                }
                              } )
                    nei.sum
                 } )
        (- ne.sum)
    }

    def totalE(configurations: ParVector[IArray[Byte]],
        interactions: Interactions ): ArraySeq[Double] = {
        val vec = configurations.map(conf => totalE(conf, interactions)).toVector
        ArraySeq(vec*)
    }

    def totalE(configurations: Vector[IArray[Byte]],
        interactions: Interactions ): ArraySeq[Double] = {
        totalE(configurations.par, interactions)
    }
/*
    def totalE(configurations: ArraySeq[IArray[Byte]],
        interactions: Interactions ): ArraySeq[Double] = {
        val vec = totalE(configurations.toVector.par, interactions)
        ArraySeq(vec*)
    }
*/
    def totalE[C <: immutable.Seq[IArray[Byte]]] (configurations: C, 
        interactions: Interactions ) (using classTagC: ClassTag[C]): ArraySeq[Double] = {
        totalE(configurations.toVector.par, interactions)
    }


  //def deltaE(configuration: mutable.IndexedSeq[Byte], 
    def deltaE(configuration: Array[Byte], 
        interactions: Interactions,
        position: Int, proposedState: Byte ) = {

        val i = position
        val currentState = configuration(i)
        val nUnits = configuration.size

        val hia = interactions.hia
        val Jijab = interactions.jijab
        val dNE =
            ArraySeq.range(0, nUnits).map ( j => {
                if ( i != j ) { 
                        val jb = configuration(j)
                        val (ij, a, b) = pairIndex(i, j, currentState, jb)
                        val (ij_, c, d) = pairIndex(i, j, proposedState, jb)
                        assert(ij == ij_)
                        Jijab(ij)(c, d) - Jijab(ij)(a, b)
                } else {
                        hia(i)(proposedState) - hia(i)(currentState)
                }
                
            } )  
        (- dNE.sum) 
    }

    def deltaE(configuration: IndexedSeq[Byte], 
        interactions: Interactions,
        position: Int, proposedState: Byte ) = {

        val i = position
        val currentState = configuration(i)
        val nUnits = configuration.size

        val hia = interactions.hia
        val Jijab = interactions.jijab
        val dNE =
            ArraySeq.range(0, nUnits).map ( j => {
                if ( i != j ) { 
                        val jb = configuration(j)
                        val (ij, a, b) = pairIndex(i, j, currentState, jb)
                        val (ij_, c, d) = pairIndex(i, j, proposedState, jb)
                        assert(ij == ij_)
                        Jijab(ij)(c, d) - Jijab(ij)(a, b)
                } else {
                        hia(i)(proposedState) - hia(i)(currentState)
                }
                
            } )  
        (- dNE.sum) 
    }

  //def dLogProposedP(configuration: mutable.IndexedSeq[Byte], 
    def dLogProposedP(configuration: Array[Byte], 
        position: Int, proposedState: Int,
        logProposedDistributionAtUnit: DenseVector[Double],
        log1_ProposedDistributionAtUnit: DenseVector[Double]
        ) = {
        val currentState = configuration(position)

        (logProposedDistributionAtUnit(proposedState) - 
                log1_ProposedDistributionAtUnit(currentState) ) - 
                (logProposedDistributionAtUnit(currentState) - 
                        log1_ProposedDistributionAtUnit(proposedState) )
    }

    def dLogProposedP(configuration: IndexedSeq[Byte], 
        position: Int, proposedState: Int,
        logProposedDistributionAtUnit: DenseVector[Double],
        log1_ProposedDistributionAtUnit: DenseVector[Double]
        ) = {
        val currentState = configuration(position)

        (logProposedDistributionAtUnit(proposedState) - 
                log1_ProposedDistributionAtUnit(currentState) ) - 
                (logProposedDistributionAtUnit(currentState) - 
                        log1_ProposedDistributionAtUnit(proposedState) )
    }

    def deltaEs(configuration: Array[Byte], 
        interactions: Interactions,
        position: Int) = {

        val i = position
        val currentState = configuration(i)
        val nUnits = configuration.size
        val nStates = interactions.hia(0).size

        val hia = interactions.hia
        val Jijab = interactions.jijab
        val dNEs = ArraySeq.range(0, nUnits).map ( j => {
                        if ( i != j ) { 
                                val jb = configuration(j)
                                val (ij, ia, ja) = pairIndex(i, j, 0, 1)
                                ArraySeq.range(0, nStates).map ( s => {
                                            if ( ia == 0 ) {
                                                Jijab(ij)(s, jb) - Jijab(ij)(currentState, jb)
                                            } else {
                                                Jijab(ij)(jb, s) - Jijab(ij)(jb, currentState)
                                            }
                                          } )
                        } else {
                                ArraySeq.range(0, nStates).map ( s => {
                                        hia(i)(s) - hia(i)(currentState)
                                  } )
                        }
                
                } )  

        ArraySeq.range(0, nStates).map ( s => {
                val dNE = ArraySeq.range(0, nUnits).map ( j => { 
                                dNEs(j)(s)
                        } )
                (- dNE.sum)
                } )
    }

    def deltaEs(configuration: IndexedSeq[Byte], 
        interactions: Interactions,
        position: Int) = {

        val i = position
        val currentState = configuration(i)
        val nUnits = configuration.size
        val nStates = interactions.hia(0).size

        val hia = interactions.hia
        val Jijab = interactions.jijab
        val dNEs = ArraySeq.range(0, nUnits).map ( j => {
                        if ( i != j ) { 
                                val jb = configuration(j)
                                val (ij, ia, ja) = pairIndex(i, j, 0, 1)
                                ArraySeq.range(0, nStates).map ( s => {
                                            if ( ia == 0 ) {
                                                Jijab(ij)(s, jb) - Jijab(ij)(currentState, jb)
                                            } else {
                                                Jijab(ij)(jb, s) - Jijab(ij)(jb, currentState)
                                            }
                                          } )
                        } else {
                                ArraySeq.range(0, nStates).map ( s => {
                                        hia(i)(s) - hia(i)(currentState)
                                  } )
                        }
                
                } )  

        ArraySeq.range(0, nStates).map ( s => {
                val dNE = ArraySeq.range(0, nUnits).map ( j => { 
                                dNEs(j)(s)
                        } )
                (- dNE.sum)
                } )
    }

    // x <= Uniform(x, y) < y
    def kernelMH1Step(
        currentState: State,
        proposedDistributions: ProposedDistributions,
        //proposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]],
        //logProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ,
        //log1_ProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ,
        interactions: Interactions,
        siteLocation: Int = -1  // In the case of -1, the siteLocation will be randomly determined.
        )( using randBasis: RandBasis ): Rand[State] = {

      import breeze.stats.distributions.Uniform

      assert( interactions == currentState.interactions )
      val energy = currentState.energy
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val proposedDistributionsAtAllUnits = proposedDistributions.proposedDistributionsAtAllUnits
      val nUnits = proposedDistributionsAtAllUnits.size
      val newkT = if ( kT == finalkT ) finalkT else {
                     val t = currentState.annealingRate * (kT/finalkT - 1.0)
                     if ( t.abs < 0.000001 ) finalkT else (t + 1.0) * finalkT
                }
      for {
        u <- Uniform(0.0, 1.0) if (u < 1.0)

        configuration = currentState.configurationArray      //currentState.configuration.clone
        uniformP <- Uniform(0.0, 1.0) if (uniformP < 1.0) 
        position = if ( siteLocation < 0 ) {
                        ( (uniformP * nUnits).floor.toInt ) % nUnits
                } else {
                        siteLocation
                }

        currentS = configuration(position)
      //maxP = math.max(1.0 - proposedDistributionsAtAllUnits(position)(currentS), 0.0)
      //uniformS <- Uniform(0.0, maxP) if ( maxP - uniformS > 0.0 || maxP <= 0.0 )
      //proposedS: Byte =  
      //                proposeState(maxP - uniformS, currentS, proposedDistributionsAtAllUnits(position))
        proposedS: Byte = proposeState(currentS, proposedDistributionsAtAllUnits(position))
        dE = deltaE(configuration, interactions, position, proposedS)  
        dLoglik = (- dE / newkT ) - 
                dLogProposedP(configuration, position, proposedS,
                   proposedDistributions.logProposedDistributionsAtAllUnits(position),
                   proposedDistributions.log1_ProposedDistributionsAtAllUnits(position) )

      } yield if ( math.log(1.0 - u) <= dLoglik) 
            { 
        /* Be careful; when a sample is taken, the configuration in a sample must be cloned. */
                configuration(position) = proposedS; 
                currentState.copy( configurationArray = configuration, 
                                   energy = energy + dE, kT = newkT, step = currentState.step + 1 )
              //State(configuration, energy + dE, newkT, 
              //        currentState.annealingRate, currentState.step + 1, currentState.finalkT)
                //println("position: %d , state:  %d => %d".format(position, currentS, configuration(position) ));
                //proposedState
            } else {
                currentState.copy( configurationArray = configuration,
                                   energy = energy, kT = newkT, step = currentState.step + 1 )
              //State(configuration, energy, newkT, 
              //        currentState.annealingRate, currentState.step + 1, currentState.finalkT)
            }
    }

    def kernelMH(
        currentState: State,
        proposedDistributions: ProposedDistributions,
        //proposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]],
        //logProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ,
        //log1_ProposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ,
        interactions: Interactions )( using randBasis: RandBasis ): Rand[State] = {

      import breeze.stats.distributions.Uniform

      assert( interactions == currentState.interactions )

      val nUnits = proposedDistributions.proposedDistributionsAtAllUnits.size

      val kernel1Step: (State, RandBasis) => Rand[State] = 
          kernelMH1Step(_: State, proposedDistributions, interactions, siteLocation = -1 )(using _: RandBasis)    // In the case of -1, the siteLocation will be randomly determined.
        //kernelMH1Step(_: State, proposedDistributions, interactions, siteLocation = -1 )( _: RandBasis)    // In the case of -1, the siteLocation will be randomly determined.

    //val dupCurrentState = currentState.copy(configuration = currentState.configuration.clone)
      /* The following must be within the for block;  the thread stops before <-.
      val iter = markovChain(currentState, kernel1Step)
      val state = if( nUnits > 1 )
                        iter.drop( nUnits - 1 ).next()
                  else
                        iter.next()
      */
      for{ r <- Uniform(0.0, 1.0)       // dummy statement
         // The configurationArray.clone is not needed but just in case.
         //dupCurrentState = currentState.copy(configurationArray = currentState.configurationArray.clone)
         //mc = markovChain(dupCurrentState, kernel1Step)
           mc = markovChain(currentState, kernel1Step)
           state = mc.drop( nUnits - 1 ).next()
      } yield { state }
    }

    def kernelGibbsWithMHStep(
        currentState: State,
        proposedDistributions: ProposedDistributions,
        interactions: Interactions ) ( using randBasis: RandBasis ) : Rand[State] = {

      import breeze.stats.distributions.Uniform

      val kernel = kernelMH1Step(_: State, proposedDistributions, interactions, _: Int)

      @annotation.tailrec
      def recCallKernel(site: Int, state: State): State = {
        if (site >= 0) {
        //val kernel = kernelMH1Step(_: State, proposedDistributions, interactions, siteLocation = site)
        //val newState = kernel(state).draw()
          val newState = kernel(state, site).draw()
          recCallKernel(site - 1, newState)
        } else {
          state
        }
      }

      assert( interactions == currentState.interactions )

    //val configuration = currentState.configuration.clone
      val nUnits = currentState.configuration.size
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val newkT = if ( kT == finalkT ) finalkT else {
                     val t = math.pow(currentState.annealingRate, nUnits.toDouble) * (kT/finalkT - 1.0)
                     if ( t.abs < 0.000001 ) finalkT else (t + 1.0) * finalkT
                }
      for {
        xx <- Uniform(0.0, 1.0) // dummy statement

      // The configurationArray.clone is not needed but just in case.
      //dupCurrentState = currentState.copy(configuration = currentState.configuration.clone, kT = newkT, finalkT = newkT)
        state = recCallKernel(nUnits - 1, currentState.copy(kT = newkT, finalkT = newkT) )
      } yield { 
        /* Be careful; state.finalkT == newkT as specified in the call.*/
        state.copy(finalkT = currentState.finalkT)
        /**/
      }
    }

    def kernelGibbs(
        currentState: State,
        interactions: Interactions )( using randBasis: RandBasis ): Rand[State] = {

      import breeze.stats.distributions.Uniform

      assert( interactions == currentState.interactions )

    //val configuration = currentState.configuration.clone
      val energy = currentState.energy
      val nUnits = currentState.configuration.size
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val newkT = if ( kT == finalkT ) finalkT else {
                     val t = math.pow(currentState.annealingRate, nUnits.toDouble) * (kT/finalkT - 1.0)
                     if ( t.abs < 0.000001 ) finalkT else (t + 1.0) * finalkT
                }
      for {
        xx <- Uniform(0.0, 1.0) // dummy statement

      //configuration = currentState.configuration.clone
        configuration = currentState.configurationArray
        dE = ArraySeq.range(0, nUnits).map( i => {
                val dEs = deltaEs(configuration, interactions, i )
                val mindEs = dEs.min
                val unnormP = dEs.map( e => math.exp( - (e - mindEs ) / newkT ) )
                val z = unnormP.sum
                val p = unnormP.map( x => x / z )
                val pDVec = DenseVector(p.toArray)
                val stateChosen: Byte =  proposeState(-10, pDVec)
                configuration(i) = stateChosen          // Be caution
                dEs(stateChosen.toInt)
              } )
        newE = energy + dE.sum
      } yield { 
        /**/
        currentState.copy(configurationArray = configuration, 
                         energy = newE, kT = newkT, step = currentState.step + nUnits, finalkT = currentState.finalkT)
      //State(configuration, newE, newkT, 
      //        currentState.annealingRate, currentState.step + nUnits, currentState.finalkT)
        /**/
      }
    }

    def markovChain(initialState: State, kernel: (State, RandBasis) => Rand[State] )(using rand: RandBasis) : Iterator[State] = {
        //val init = State(initialState.configuration.clone, initialState.energy,
        //               initialState.kT, initialState.annealingRate, initialState.step, initialState.finalkT)
        val init =
            if ( initialState.configurationArray.size == 0 ) {
              assert( initialState.step >= 0 )
              initialState.copy( configurationArray = Array(initialState.configuration*) )
            } else { 
              assert( initialState.step > 0 )
              initialState
            //assert(initialState.configurationArray.size == 0, "Must be for markovChain: initialState.configurationArray.size == 0\n" )
            }
        //
        // import breeze.stats.distributions.MarkovChain
        // MarkovChain(initState)(kernel).steps                      // for breeze 0.13.1
        // or
        val kern = kernel(_: State, rand).draw()                      // return T rather than Rand[T]
        val iter = scala.collection.Iterator.iterate(init)(kern)
        val initValue = iter.next()             // remove the initial value.
        iter
        // or
        // val kern = kernel(_: State).draw()                   // return T rather than Rand[T]
        // scala.collection.Iterator.iterate(init)(kern).drop(1)        // remove the initial value.
    }

   /*
    def sampling(markovProcess: Process[State], 
                nInitialIterations: Int, 
                everyNIterations: Int,
                nSamples: Int 
                ): ArraySeq[State] = {

        val iterator = markovProcess.steps.drop(nInitialIterations)
        val samples = ArraySeq.fill(nSamples)(0).map( _ => {
                val s = iterator.drop(everyNIterations - 1).next()      // Error ?
                State(s.configuration.clone, s.energy, s.kT, s.annealingRate, s.step, s.finalkT )
        } )
        samples
    }
   */

    /*
    def sampling(markovChain: Iterator[State], 
                nInitialIterations: Int, 
                everyNIterations: Int,
                nSamples: Int 
                ): ArraySeq[State] = {

        //val iterator = markovChain.drop(nInitialIterations)
        val iterator = if ( nInitialIterations > 0 ) {
                val i = (nInitialIterations - 1 ) / 2
                val j = nInitialIterations - 1 - i
                val state = markovChain.drop(i).next()
                require( state.kT == state.finalkT )
                markovChain.drop(j)                             // Error; see the manual for Iterator
                //
        } else {
                markovChain
        }
        val samples = ArraySeq.fill(nSamples)(0).map( _ => {
                val s = iterator.drop(everyNIterations - 1).next()              // Error ?
                s.copy(configuration = s.configuration.clone)
                //State(s.configuration.clone, s.energy, s.kT, s.annealingRate, s.step, s.finalkT )
        } )
        samples
    }
   */

    def sampling(markovChain: Iterator[State], 
                nInitialIterationsPerUnit: Int, 
                everyNIterationsPerUnit: Int,
                nSamples: Int 
                ): Tuple2[ArraySeq[State], Iterator[State]] = {

        //val iterator = markovChain.drop(nInitialIterationsPerUnit)
        val iterator = if ( nInitialIterationsPerUnit > 0 ) {
              /*
                val i = (nInitialIterationsPerUnit - 1 ) / 2
                val j = nInitialIterationsPerUnit - 1 - i
                Range(0, i).foreach( _ => { val drop = markovChain.next() } )
                val state = markovChain.next()
                require( state.kT == state.finalkT )
                Range(0, j).foreach( _ => { val drop = markovChain.next() } )
              */
                Range(0, nInitialIterationsPerUnit - 1).foreach( _ => { val drop = markovChain.next() } )
                val state = markovChain.next()
                require( state.kT == state.finalkT, "MC: does not reach final kT.\n" )
                markovChain
              //Or
              //val markovChain2 = markovChain.drop(i)
              //val state = markovChain2.next()
              //require( state.kT == state.finalkT )
              //val markovChain3 = markovChain2.drop(j)
              //markovChain3
        } else {
                markovChain
        }
        val samples = ArraySeq.fill(nSamples)(0).map( _ => {
                Range(0, everyNIterationsPerUnit - 1).foreach( _ => { val drop = iterator.next() } )
                val s = iterator.next()
              //assert( s.kT == s.finalkT, "MC: does not reach final kT.\n")
                s.copy(configuration = IArray(s.configurationArray*), configurationArray = Array[Byte]() )
              //s.copy(configuration = ArraySeq(s.configurationArray*), configurationArray = s.configurationArray.clone )
              //State(s.configuration.clone, s.energy, s.kT, s.annealingRate, s.step, s.finalkT )
        } )
        // Or recursive call must be used.

      //assert( samples(0).configuration.zip(samples(0).configurationArray).map( (x, y) => if(x == y) 0 else 1 ).sum == 0 )
      //val samplesRevised = samples.map{ s => s.copy(configurationArray = Array[Byte]() ) }

        (samples, iterator)
    }

    def nNonEquilibrium(independentSamplings: IndependentSamplings ) = {

     import breeze.stats.meanAndVariance

     val nRuns = independentSamplings.size
     val nSamples = independentSamplings(0).size

     if ( nSamples <= 1 ) {
      0
     } else { 
      val firstSamples = independentSamplings.map{ s => s(0).energy }.toVector
      val lastSamples = independentSamplings.map{ s => s(nSamples - 1).energy }.toVector
      val statFirstSamples = meanAndVariance(firstSamples)
      val statLastSamples = meanAndVariance(lastSamples)
      val avFirstSamples = statFirstSamples.mean
      val avLastSamples = statLastSamples.mean

      //val stdDifAv = math.sqrt( (statFirstSamples.variance + statLastSamples.variance) / nRuns )
      //val nDifAv = if ( (avFirstSamples - avLastSamples).abs > stdDifAv * 2.0) 1 else 0

      val stdDif = math.sqrt( (statFirstSamples.variance + statLastSamples.variance) )
      val difFreq = firstSamples.zip(lastSamples).map( m => if ( (m._1 - m._2).abs > stdDif * 2.0) 1 else 0 )

    //val nDif = difFreq.take(5).sum
      val nDif = { //Normal distr. is assumed.
          if ( difFreq.sum / (difFreq.size).toDouble > 0.0455 ) 1 else 0  
      }

    //nDif

      // t-test for autocorrelation

    //val intervalAC = math.min(900, nSamples - 10)
    //val intervalAC = nSamples - 10
      val intervalAC = {
          val nSForAve = math.max( 2, math.min(nSamples - 1, 100) )  
          nSamples - nSForAve
      }

      if ( intervalAC <= 0 ) {
        nDif

      } else {

        def isCorrelated(r2: Double, n: Int, thresholdT2ForBothSides: Double = 14.691889 ) = {
                // t= 3.833 (t2 = 14.691889) : for 99.5 % for the t distribution of freedom = 8
                // t= 5.041 (t2 = 25.411681) : for 99.9 % for the t distribution of freedom = 8
                // t= 4.587 (t2 = 21.040569) : for 99.9 % for the t distribution of freedom = 10
                // t= 3.291 (t2 = 10.830681) : for 99.9 % for the t distribution of freedom = inf
                // t= 2.807 (t2 = 7.879249) : for 99.5 % for the t distribution of freedom = inf
                // t= 1.960 (t2 = 3.8416) : for 95 % for the t distribution of freedom = inf
                // r= 0.0520 (n=1000): for 95 % for the t distribution of freedom = inf
                
                        //val t2 = t * t 
                        //val r2 = t2 / (n - 2 + t2)
                        val t2 = r2 * ( n - 2) / ( 1.0 - r2 )

                        if ( t2 > thresholdT2ForBothSides )
                                true
                        else 
                                false
        }

        /* */
        val independentSampledEnergies1 = independentSamplings.map( s => s.take(nSamples - intervalAC).map( si => si.energy ) ) 
        val independentSampledEnergies2 = independentSamplings.map( s => s.drop(intervalAC).map( si => si.energy ) ) 
        val independentSampledAutoCorr = independentSampledEnergies1.zip(independentSampledEnergies2).map(x => 
                                {       val (e1, e2) = x
                                        val corrMat = breeze.stats.corrcoeff ( DenseMatrix(e1, e2).t )
                                        corrMat(0, 1)
                                } )
        val autoCorrelated  = independentSampledAutoCorr.map ( c => if (isCorrelated(c, nSamples - intervalAC)) 1 else 0 )
        val nCorr = autoCorrelated.take(5).sum
        /* */

        /*
        val nRunsArraySeq = ArraySeq.range(0, nRuns)
        val nACArraySeq = ArraySeq.range(0, nSamples - intervalAC)
        val independentSampledEnergies = independentSamplings.map( s => s.map( si => si.energy ) )
        val statIndependentSampledEnergies = independentSampledEnergies.map(s => meanAndVariance(s) )
        val independentSampledDEnergies = nRunsArraySeq.map( s => independentSampledEnergies(s).map
                                                (ei => ei - statIndependentSampledEnergies(s).mean) )
        val independentSampledAutoCorr = nRunsArraySeq.map( s => {
                        val corr = nACArraySeq.map( i => {
                          independentSampledDEnergies(s)(i) * independentSampledDEnergies(s)(i + intervalAC) 
                          } ).sum / (nSamples - intervalAC)
                        val corrCoeff = corr / statIndependentSampledEnergies(s).variance
                        corrCoeff * corrCoeff
                } )

        val autoCorrelated  = independentSampledAutoCorr.map ( c => if (isCorrelated(c, nSamples - intervalAC)) 1 else 0 )
        val nCorr = autoCorrelated.take(10).sum
        */

      //scala.sys.process.stderr.println( (nDifAv, nCorr) )

        nDif + nCorr
      }
     }
    }

    def randomConfiguration(nUnits: Int, nStatesofUnit: Int,
        proposedDistributionsAtAllUnits: ArraySeq[DenseVector[Double]] ) = {

    //import breeze.stats.distributions.Rand.FixedSeed.*            //for scala 3
      import SetSeedForRand.randBasis
      import breeze.stats.distributions.Uniform

      // Probably: x <= Uniform(x, y) < y
      def uniformGExLTy(x: Double, y: Double): Rand[Double] = 
        for {
                r <- Uniform(x, y)      if r < y
        } yield { r }

      val conf = IArray.range(0, nUnits)
      if ( proposedDistributionsAtAllUnits == null ) {
          conf.map( _ => { 
                ( ( ( uniformGExLTy(0.0, 1.0).sample() * nStatesofUnit ).floor.toInt ) % nStatesofUnit ).toByte
                } )
      } else {
        //conf.zipWithIndex.map( x => {
        //      val position = x._2
          if ( proposedDistributionsAtAllUnits.size <= 0 ) {
                    conf.map( position => {
                        ( ( ( uniformGExLTy(0.0, 1.0).sample() * nStatesofUnit ).floor.toInt ) % nStatesofUnit ).toByte
                    } )
          } else if ( proposedDistributionsAtAllUnits.size <= 1 ) {
                    conf.map( position => {
                        proposeState(currentState = -1 , proposedDistributionsAtAllUnits(0))
                    } )
          } else {
                    conf.map( position => {
                        proposeState(currentState = -1 , proposedDistributionsAtAllUnits(position))
                    } )
          }
      }
    }

    def frequenciesAtUnitInSamples( configurations: IndexedSeq[IArray[Byte]],                
                nStatesOfUnit: Int
                ): ArraySeq[DenseVector[Double]] = {

          val nUnits = configurations(0).size
        //val pia = ArraySeq.fill(nUnits)(0).map( _ => (new DenseVector[Double](nStatesOfUnit)) )
          val pia = ArraySeq.fill(nUnits)(0).map( _ => DenseVector.fill(nStatesOfUnit)(0.0) )

          val weight = 1.0 / configurations.size

        /* configurations cannot be ParVector. */
          configurations.foreach( s => {
                s.zipWithIndex.foreach( c => pia(c._2)(c._1) += weight )
          } )
          pia
    }

    def frequenciesAtUnitInSamples( independentSamplings: ArraySeq[ArraySeq[State]],                
                nStatesOfUnit: Int
                ): ArraySeq[DenseVector[Double]] = {

          val configurations = independentSamplings.flatten.map( s => s.configuration )

          frequenciesAtUnitInSamples( configurations, nStatesOfUnit )
    }

    def frequenciesAtUnitInSamples( independentSamplings: Vector[ArraySeq[State]],                
                nStatesOfUnit: Int
                ): ArraySeq[DenseVector[Double]] = {

          val configurations = independentSamplings.flatten.map( s => s.configuration )

          frequenciesAtUnitInSamples( configurations, nStatesOfUnit )
    }

    def frequenciesAtUnitInSamples( independentSamplings: IndependentSamplings,                
                nStatesOfUnit: Int
                ): ArraySeq[DenseVector[Double]] = {

        //val configurations = independentSamplings.toVector.flatten.map( s => s.configuration )
        //frequenciesAtUnitInSamples( configurations, nStatesOfUnit )

          frequenciesAtUnitInSamples( independentSamplings.toVector, nStatesOfUnit )
    }

    def pairwiseFrequenciesInSamples( configurations: IndexedSeq[IArray[Byte]],
                nStatesOfUnit: Int
                ): ArraySeq[DenseMatrix[Double]] = {

          val nUnits = configurations(0).size
          val nPairs = ( nUnits * ( nUnits - 1)) / 2

        //val pijab = ArraySeq.fill(nPairs)(0).map( _ => new DenseMatrix[Double](nStatesOfUnit, nStatesOfUnit) )
          val pijab = ArraySeq.fill(nPairs)(0).map( _ => DenseMatrix.fill(nStatesOfUnit, nStatesOfUnit)(0.0) )

          val weight = 1.0 / configurations.size

        /* configurations cannot be ParVector. */
          configurations.foreach( configuration => {
                cfor(0)(i => i < nUnits, i => i + 1)( i => {
                        cfor(0)(j => j < i, j => j + 1)( j => {
                          val c = configuration(i)
                          val d = configuration(j)
                          val (ij , a, b ) = pairIndex(i, j, c, d)
                          pijab(ij)(a, b) += weight
                        } )
                } )
          } )
          pijab
    }

    def pairwiseFrequenciesInSamples( independentSamplings: ArraySeq[ArraySeq[State]],               
                nStatesOfUnit: Int
                ): ArraySeq[DenseMatrix[Double]] = {

          val configurations = independentSamplings.flatten.map( s => s.configuration )

          pairwiseFrequenciesInSamples( configurations, nStatesOfUnit)
    }

    def pairwiseFrequenciesInSamples( independentSamplings: Vector[ArraySeq[State]],               
                nStatesOfUnit: Int
                ): ArraySeq[DenseMatrix[Double]] = {

          val configurations = independentSamplings.flatten.map( s => s.configuration )

          pairwiseFrequenciesInSamples( configurations, nStatesOfUnit)
    }

    def pairwiseFrequenciesInSamples( independentSamplings: IndependentSamplings,               
                nStatesOfUnit: Int
                ): ArraySeq[DenseMatrix[Double]] = {

        //val configurations = independentSamplings.toVector.flatten.map( s => s.configuration )
        //pairwiseFrequenciesInSamples( configurations, nStatesOfUnit)
 
          pairwiseFrequenciesInSamples( independentSamplings.toVector, nStatesOfUnit)
    }

    def paramFileToHashMap(file: File): HashMap[String, String] = {
        if ( file.exists() ) {
                val bufsource = Source.fromFile(file)
                val lines = bufsource.getLines().toList
                bufsource.close()
                val tuplelist = lines.map { line =>
                        val l = line.strip()
                        val items = l.split(Array(' ', '\t'))
                        val k = items(0)
                        val v = l.substring( k.length ).strip()
                        (k, v)
                }
                //HashMap[String, String]( tuplelist:_* )       // for scala 2.13
                HashMap[String, String]( tuplelist* )           // for scala 3
        } else {
                HashMap[String, String]()
        }
    }

    def compoundOrderStringToMap(compoundorderstring: String,
                caseSensitive: Boolean = false ) = {

                val comporderstr = compoundorderstring
                val case_sensitive = caseSensitive

                val compOrderMap =
                  if ( ! case_sensitive ) {
                    val compoundUpper = comporderstr.toUpperCase.toArray
                    val compoundLower = comporderstr.toLowerCase.toArray
                  //compoundUpper.zipWithIndex.toMap ++ compoundLower.zipWithIndex.toMap
                    HashMap(compoundUpper.zipWithIndex*) ++ HashMap(compoundLower.zipWithIndex*)
                  } else {
                  //comporderstr.toArray.zipWithIndex.toMap
                    HashMap(comporderstr.toArray.zipWithIndex*)
                  }

                val new_compoundOrderMap =
                  if (compOrderMap.contains('-') && ! compOrderMap.contains('.') )
                        compOrderMap + ('.' -> compOrderMap('-'))
                  else if (compOrderMap.contains('.') && ! compOrderMap.contains('-'))
                        compOrderMap + ('-' -> compOrderMap('.'))
                  else
                        compOrderMap
                new_compoundOrderMap
    }

    def compoundOrderNumber( compound: Char, other: Int, compoundOrderMap: Map[Char, Int] ) = {
          if(compoundOrderMap.contains(compound)) {
                compoundOrderMap(compound)
          } else {
                require( other >= 0, "Unknown compound character %c".format(compound))
                other
          }
    }

    def readFasta(lines: collection.Iterator[String],
                  stateOrderString: String, caseSensitive: Boolean = false ):
                  Tuple3[ArraySeq[String], ArraySeq[String], ArraySeq[IArray[Byte]] ] = {

        val (lines0, lines1) = lines.duplicate
        val headerLine = ("""^>[ \t]*(.*)[ \t]*$""").r

        @annotation.tailrec
        def readEntries (lines: collection.Iterator[String], line: String,
                         idList: List[String], seqList: List[String] ): 
                         Tuple3[collection.Iterator[String], List[String], List[String]] = {

          @annotation.tailrec
          def readSeq (lines: collection.Iterator[String],
                         seqList: List[String],
                         seq: String ): Tuple3[collection.Iterator[String], List[String], String] = {
              if ( lines.hasNext ) {
                // stip: removes kanji-blank as well as ascii blanks. trim: removes ascii blanks. 
                val line = lines.next().strip()
                if ( headerLine.findFirstMatchIn(line).nonEmpty ) {
                  val newSeqList = seqList.+:(seq)
                  (lines, newSeqList, line)
                } else {
                  val newSeq = seq + line
                  readSeq (lines, seqList, newSeq)
                }
              } else {
                val newSeqList = seqList.+:(seq)
                (lines, newSeqList, "")
              }
          }

          val id = headerLine.findFirstMatchIn(line).get.group(1)
          val newIdList = idList.+:(id)

          val (lines1, newSeqList, newLine) =
            readSeq (lines, seqList, seq = "")

          if ( newLine == "" ) {
            (lines1, newIdList, newSeqList)
          } else {
            readEntries (lines1, newLine, newIdList, newSeqList)
          }

        }

        val lines11 = lines1.dropWhile(line => headerLine.findFirstMatchIn(line).isEmpty )

        val (lines12, idList, seqList ) =
          if ( lines11.hasNext ) {
            val line = lines11.next()
            readEntries (lines11, line, idList = List[String](), seqList = List[String]())
          } else {
            (lines11, List[String](), List[String]() )
          }
        val idArraySeq = ArraySeq(idList*).reverse
        val seqArraySeq = ArraySeq(seqList*).reverse
        assert( idArraySeq.size == seqArraySeq.size )

        val stateOrderMap = compoundOrderStringToMap(stateOrderString, caseSensitive = caseSensitive )
        
        val numberArraySeq = seqArraySeq.map( s => 
                           IArray(s.toArray*).map( c =>
                            {
                              val b = compoundOrderNumber( c, other = stateOrderMap('-'), compoundOrderMap = stateOrderMap )
                              b.toByte
                            }
                           )
                          )
        (idArraySeq, seqArraySeq, numberArraySeq)
    }

    def readSeqWeights(lines: collection.Iterator[String], idArraySeq: ArraySeq[String] ): ArraySeq[Double] = {

        val (lines0, lines1) = lines.duplicate
      //val headerLine = ("""^#[ \t]*(id)[ \t]*$""").r
        val weightLine = ("""^[^#]""").r

        @annotation.tailrec
        def readEntries (lines: collection.Iterator[String], line: String,
                         idPos: Int, wList: List[String] ):
                         Tuple3[collection.Iterator[String], Int, List[String]] = {

          def readWeight(lines: collection.Iterator[String], idPos: Int,
                         wList: List[String],
                         seq: String ): Tuple3[collection.Iterator[String], List[String], String] = {
              if ( lines.hasNext ) {
                // stip: removes kanji-blank as well as ascii blanks. trim: removes ascii blanks. 

                val line = lines.next()
                val newW =
                  if ( weightLine.findFirstMatchIn(line).nonEmpty ) {
                    seq + line
                  } else {
                    require( line(0) == '#', "The next line of ID must be a weight line." )
                    seq
                  }
                val newWList = wList.+:(newW)
                if ( idPos + 1 >= idArraySeq.size ) {
                  val lines2 = lines.dropWhile(line => weightLine.findFirstMatchIn(line).isEmpty )
                  require(lines2.hasNext == true, "too many weight lines; idPos = %d".format(idPos) )
                  (lines2, newWList, "")
                } else {
                  val headerLine = ("""^#[ \t]*(""" + idArraySeq(idPos + 1) + """)[ \t]*$""").r
                  val lines2 = lines.dropWhile(line => headerLine.findFirstMatchIn(line).isEmpty )
                  if ( lines2.hasNext ) {
                    (lines2, newWList, lines2.next)
                  } else {
                    (lines2, newWList, "")
                  }
                }
              } else {
                assert(lines.hasNext == true, "Not ID line; idPos = %d".format(idPos) ) 
                val newWList = wList.+:(seq)
                (lines, newWList, "")
              }
          }

          val headerLine = ("""^#[ \t]*(""" + idArraySeq(idPos) + """)[ \t]*$""").r
          val id = headerLine.findFirstMatchIn(line).get.group(1)

          val (lines1, newWList, newLine) =
            readWeight (lines, idPos, wList, seq = "")

          if ( newLine == "" ) {
            (lines1, idPos + 1, newWList)
          } else {
            readEntries (lines1, newLine, idPos + 1, newWList)
          }

        }

        val idPos = 0
        val headerLine = ("""^#[ \t](""" + idArraySeq(idPos) + """)[ \t]$""").r
        val lines11 = lines1.dropWhile(line => headerLine.findFirstMatchIn(line).isEmpty )

        val (lines12, nIDs, wList ) =
          if ( lines11.hasNext ) {
            val line = lines11.next()
            readEntries (lines11, line, idPos, wList = List[String]())
          } else {
            (lines11, idPos, List[String]() )
          }
        val wArraySeq = ArraySeq(wList.toArray*).reverse.map( x => x.toDouble )
        
        assert( idArraySeq.size == nIDs && wArraySeq.size == nIDs, "Data is inconsistent; nIDs = %d".format(nIDs) )
        
        wArraySeq
    }

  } // End of object MCMC

  class MCMC( val ioFiles: HashMap[String, File], 
              val stateOrderString: String,
              val interactions: MCMC.Interactions,
              val proposedPia: ArraySeq[DenseVector[Double]] ) {

    import miyazawa.potts.MCMC
    import miyazawa.potts.MCMC.*

    //assert( nUnits == proposedDistributionAtUnit.size )
    //assert( nStatesOfUnit == proposedDistributionAtUnit(0).size )
    val nUnits = interactions.hia.size
    val nStatesOfUnit = interactions.hia(0).size
    require( interactions.jijab.size == (nUnits * (nUnits - 1)) / 2 &&
              interactions.jijab(0).size == nStatesOfUnit * nStatesOfUnit )
    require( proposedPia.size == nUnits && proposedPia(0).size == nStatesOfUnit )

    //val (proposedDistributionsAtAllUnits, 
    //  logProposedDistributionsAtAllUnits, log1_ProposedDistributionsAtAllUnits ) = 
    val proposedDistributions = this.proposeDistributions(proposedPia)

    def this(ioFiles: HashMap[String, File],
                stateOrderString: String,
                hia: ArraySeq[ArraySeq[Double]],
                jijab: ArraySeq[ArraySeq[ArraySeq[Double]]],
                proposedPia: ArraySeq[ArraySeq[Double]] ) = {
      this(ioFiles, stateOrderString,
                MCMC.Interactions( hia.map( hi => DenseVector(hi.toArray) ), 
                //jijab.map( Jij => DenseMatrix(Jij:_*) )  ),                                           // for scala 2.11
                  jijab.map( Jij => DenseMatrix(Jij*) )  ),                                             // for scala 2.13
                //jijab.map( Jij => new DenseMatrix(Jij(0).size, Jij.size, Jij.flatten.toArray ).t )  ),        // for scala 2.13
                proposedPia.map(pi => DenseVector(pi.toArray)) )
    }

    def this(ioFiles: HashMap[String, File], 
                stateOrderString: String,
                interactions: MCMC.Interactions  ) = {
        this( ioFiles, stateOrderString, interactions,
                {
                  import breeze.linalg.{max, sum}
                  import breeze.numerics.exp

                  val hia = interactions.hia
                  val p = 1.0           //0.6
                  val uniform = (1.0/ hia(0).size) * (1.0 - p)

                  hia.map( hi => {
                        val maxhi = breeze.linalg.max(hi)
                        val exphi = breeze.numerics.exp(hi - maxhi)     //hi.map( math.exp(_ - maxhi) )
                        val z = breeze.linalg.sum(exphi)
                        exphi.map{ x => (x / z) * p + uniform }
                  } )
                }
        )
    }

    def this(ioFiles: HashMap[String, File], 
                stateOrderString: String,
                hia: ArraySeq[ArraySeq[Double]],
                jijab: ArraySeq[ArraySeq[ArraySeq[Double]]] ) = {
        this(ioFiles, stateOrderString, hia, jijab,
                {
                  val p = 1.0  //0.6
                  val uniform = (1.0/ hia(0).size) * (1.0 - p)

                  hia.map( hi => {
                        val maxhi = hi.max
                        val exphi = hi.map(x => math.exp( x - maxhi) ) 
                        val z = exphi.sum
                        exphi.map{ x => (x / z) * p + uniform }
                  } )
                }
        )
    }

    def proposeDistributions(proposedPia: ArraySeq[DenseVector[Double]]) = {
      // must be corrected with pseudocounts.

      proposedPia.foreach(pi => pi.foreach(p => require( p > 0.0 , "proposedPia must be positive and less than 1." ) ) )
      /*
        val ratio = observedN / (observedN + pseudoN) 
        val pseudo = (1.0 - ratio) / proposedFia(0).size
      //val proposedP = proposedPia.map ( fi => fi.map( fia => fia * ratio + pseudo ) )
        val proposedP = proposedPia.map ( fi => fi * ratio + pseudo )
      */
      //val logProposedP = proposedPia.map(pi => pi.map(pia =>  math.log(pia)) )
        val logProposedP = proposedPia.map(pi => breeze.numerics.log(pi) )
      //val log1_proposedP = proposedPia.map(pi => pi.map( pia => math.log(1.0 - pia) ) )
        val log1_proposedP = proposedPia.map(pi => breeze.numerics.log(-pi + 1.0) )
        ProposedDistributions(proposedPia, logProposedP, log1_proposedP)
    }

    def kernelMH(currentState: State )( using randBasis: RandBasis ): Rand[State] =
          MCMC.kernelMH(currentState, 
                this.proposedDistributions,
                //proposedDistributionsAtAllUnits, 
                //logProposedDistributionsAtAllUnits,
                //log1_ProposedDistributionsAtAllUnits,
                this.interactions )

    def kernelGibbsWithMHStep(currentState: State )( using randBasis: RandBasis ): Rand[State] =
          MCMC.kernelGibbsWithMHStep(currentState,
                this.proposedDistributions, this.interactions )

    def kernelGibbs(currentState: State )( using randBasis: RandBasis ): Rand[State] =
          MCMC.kernelGibbs(currentState, 
                this.interactions )

  //def markovChain(initial: State, kernel: (State, RandBasis) => Rand[State] = kernelGibbsWithMHStep(_)(_) )( using randBasis: RandBasis ): Iterator[State] =
    def markovChain(initial: State, kernel: (State, RandBasis) => Rand[State] = kernelGibbsWithMHStep(_)(using _: RandBasis) )( using randBasis: RandBasis ): Iterator[State] =
          MCMC.markovChain(initial, kernel)
        
    def runMC(  initialStates: InitialMCStates,
                nInitialIterationsPerUnit: Int = 100, 
                everyNIterationsPerUnit: Int = 5,
                nSamples: Int = 1000,
              //nIndependentMC: Int = 10,
                initialT: Double = 1.2,
                annealingRate: Double = 0.99,
                finalT: Double = 1.0,
                maxExtendedIterations: Int = 10,
              //kernel: (State, RandBasis) => Rand[State] = kernelGibbsWithMHStep(_)(_)
                kernel: (State, RandBasis) => Rand[State] = kernelGibbsWithMHStep(_)(using _: RandBasis)
                ) = {

        assert( nInitialIterationsPerUnit >= 0 )
        assert( everyNIterationsPerUnit >= 1 )
        assert( nSamples > 0 )
        assert( maxExtendedIterations >= 1 )

        val nInitialIterations = nInitialIterationsPerUnit * this.nUnits
        val everyNIterations = everyNIterationsPerUnit * this.nUnits

        val nIndependentMC = initialStates.size

        // T(n) = r^n T(0) + ( 1 - r^n) T(final) = (1 + t) T(final)
        val (kT, rate) = {
                val n = (nInitialIterations * 0.8 ).floor
                if ( initialT == finalT || n <= 2 ) {
                  (finalT, 0.0)
                } else {
                  val maxRate = math.exp( (math.log(0.000001) - math.log( (initialT / finalT - 1.0).abs + 0.000001) ) / n )
                               // (nInitialIterations / 5.0) )  /* (nInitialIterations / 10.0) ) */
                  val r = if ( annealingRate > maxRate ) maxRate else  annealingRate 
                  (initialT, r)
                }
        }

        //val independentMarkovProcesses = ArraySeq.fill(nIndependentMC)(0).map( _ => {
        //      val configuration = MCMC.randomConfiguration(this.nUnits, this.nStatesOfUnit)

        val revInitialStates = {
                  val initialSVec = initialStates.toVector
                  // To reproduce the same result
                  val initialCVec = initialSVec.map( s => {
                        val x = s.configuration
                      //val configuration = if( x == null || x.size <= 0 ) {
                        val configuration = if( x.size <= 0 ) {
                                MCMC.randomConfiguration(this.nUnits, this.nStatesOfUnit, 
                                        this.proposedPia )
                                     // this.proposedDistributions.proposedDistributionsAtAllUnits )
                          } else {
                                require( x.size == this.nUnits )
                                x // x.clone for Array but not for immutable.ArraySeq
                          }
                        configuration
                  } )
                  val initialSPar = initialStates.zip(initialCVec).map( (s, configuration) => {
                        val energy = MCMC.totalE( configuration, this.interactions)
                        val initialS = s.copy(configuration = configuration, interactions = this.interactions,
                                            energy = energy, kT = kT, annealingRate = rate, step = 0, finalkT = finalT )
                        initialS
                  } )
                  initialSPar
        }
    
        //val independentMarkovProcesses = 
        //        initialStates.map( initial => markovChain(initial) )

        @annotation.tailrec
        def equilSampling(mthIteration: Int, 
                //independentMC: ArraySeq[Process[State]],
                initialStates: InitialMCStates,
                nInitialIterationsPerUnit: Int,
                everyNIterationsPerUnit: Int ): (Int, Int, IndependentSamplings ) = {

          val randBasisIA = setRandBasisIArrayWithRandomSeeds( initialStates.size )
          val independentMC = 
                    initialStates.zipWithIndex.map{ (initial, i) =>
                        given randBasis: RandBasis = randBasisIA(i)
                        MCMC.markovChain(initial, kernel)(using randBasis: RandBasis)
                      //MCMC.markovChain(initial, kernel)(randBasis)
                    }       // markovChain: Iterator[State]

                //(new IndependentMC(initialStates.size)).zipWithIndex.map( x => { MCMC.markovChain(initialStates(x._2), kernel) } )

          @annotation.tailrec
          def equilSamplingWithIterator(nthIteration: Int, 
                independentMC: IndependentMC,
                nInitialIterationsPerUnit: Int, 
                newNInitialIterationsPerUnit: Int,
                everyNIterationsPerUnit: Int ): (Int, Int, IndependentSamplings, IndependentMC) = {
        
            val independentSamplingsAndIterator = independentMC.map( mc => {
                MCMC.sampling(mc, nInitialIterationsPerUnit, everyNIterationsPerUnit * nthIteration, nSamples)
            } )
            val independentSamplings = independentSamplingsAndIterator.map( x => { x._1 } )
            val independentMCcont = independentSamplingsAndIterator.map( x => { x._2 } )

            val nNonEquil = nNonEquilibrium(independentSamplings)

            if ( newNInitialIterationsPerUnit > 0 || nNonEquil <= 0 || nthIteration >= maxExtendedIterations ) {

              (nNonEquil, nthIteration , independentSamplings, independentMCcont)

            } else {
              equilSamplingWithIterator(
                  nthIteration + 1, independentMCcont,
                  newNInitialIterationsPerUnit, newNInitialIterationsPerUnit, everyNIterationsPerUnit ) 
            }

          }

          val (nNonEquil, nthIteration, independentSamplings, independentMCcont) =
               equilSamplingWithIterator( mthIteration, independentMC,
                   nInitialIterationsPerUnit, newNInitialIterationsPerUnit = 0, everyNIterationsPerUnit = everyNIterationsPerUnit ) 

          val newNInitialIterationsPerUnit = 0  // 0 means that heating-again is not done.
          val newNInitialIterations = newNInitialIterationsPerUnit * this.nUnits

          //val nNonEquil = nNonEquilibrium(independentSamplings)

          if ( nNonEquil <= 0 || nthIteration >= maxExtendedIterations ) {

                (nNonEquil, nthIteration, independentSamplings)

          } else {
          // slightly heat 
            val (newkT, newRate) = {
                  if ( newNInitialIterations >= 1 ) {
                    val t = 1.2 * finalT
                    val maxRate = math.exp( (math.log(0.000001) - math.log( (t / finalT - 1.0).abs + 0.000001)) / 
                                        (newNInitialIterations / 5.0) ) /* (newNInitialIterations / 10.0) ) */
                    (t, maxRate)
                  } else {
                    (finalT, 0.0)
                  }
                }

            val lastStates = independentSamplings.map{ s => 
                                val last = s(s.size - 1)
                              //last.copy(configurationArray = last.configurationArray.clone, kT = newkT, annealingRate = newRate )
                                last.copy(kT = newkT, annealingRate = newRate )
                             }

            //val independentMarkovProcesses = lastStates.map ( s => markovChain(s))

            //equilSampling(nthIteration + 1, independentMarkovProcesses, nInitialIterations ) 
            //equilSampling(nthIteration + 1, independentMarkovProcesses, nInitialIterations = 0 ) 

            equilSampling(nthIteration + 1, lastStates,
                nInitialIterationsPerUnit = newNInitialIterationsPerUnit, everyNIterationsPerUnit = everyNIterationsPerUnit  ) 
          }

        }

        //equilSampling(1, independentMarkovProcesses, nInitialIterations ) 
        val (nNonEquil, nExtendedIterations, independentSamplings) =
            equilSampling(1, revInitialStates, nInitialIterationsPerUnit, everyNIterationsPerUnit )
        (revInitialStates, nNonEquil, nExtendedIterations, independentSamplings)
    }

    def frequenciesAtUnitInSamples( independentSamplings: IndependentSamplings ): ArraySeq[DenseVector[Double]] =
        MCMC.frequenciesAtUnitInSamples( independentSamplings, this.nStatesOfUnit)

    def pairwiseFrequenciesInSamples( independentSamplings: IndependentSamplings ): ArraySeq[DenseMatrix[Double]] =
        MCMC.pairwiseFrequenciesInSamples( independentSamplings, this.nStatesOfUnit)

  } // End of class MCMC


  object BM {

    case class BMInteractions(val phia: ArraySeq[DenseVector[Double]] = ArraySeq[DenseVector[Double]](),
               val phijab: ArraySeq[DenseMatrix[Double]] = ArraySeq[DenseMatrix[Double]](),
               val gauge: String = "unused")
    /*
    case class BMState (bmInteractions: BMInteractions, 
                        v: BMInteractions, m: BMInteractions, grad: BMInteractions, 
                        learningRate: LearningRate, betaV: Double, betaM: Double, eps: Double, step: Int )
        //betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6
    */
    case class BMState (val gradientDescentMethod: String = "",
                        val bmInteractions: BMInteractions = BMInteractions(), 
                        val grad: BMInteractions = BMInteractions(),
                        val v: BMInteractions = BMInteractions(), val m: BMInteractions = BMInteractions(),
                        val betaV: Double = 0.999, val betaM: Double = 0.9,
                        val eps: Double = 1.0e-6,
                        val learningRate: LearningRate = LearningRate(),
                        val learningRateForRPROPLR: LearningRateForRPROPLR = LearningRateForRPROPLR(), 
                      //val minLearningRate: Double, val maxLearningRate: Double,
                      //val rateDecrease: Double = 1.0, val rateIncrease: Double = 1.0,
                        val learningRates: BMInteractions = BMInteractions(),
                        val step: Int = scala.Int.MinValue, val stepOffset: Int = 0,   // step < 0 means this object is used as a container of initial MCMC.Statess by independentSamplings 
                        val initialMCStates: MCMC.InitialMCStates = ParVector[MCMC.State](),
                        val independentSamplings: MCMC.IndependentSamplings = ParVector[ArraySeq[MCMC.State]](),
                        val ensembleAverages: MCMC.EnsembleAverages = MCMC.EnsembleAverages(),
                        val leakyAveP: MCMC.EnsembleAverages = MCMC.EnsembleAverages(), val betaLAP: Double = 0.0,
                        val leakyAvePForKL: MCMC.EnsembleAverages = MCMC.EnsembleAverages(), val betaLAPForKL: Double = 0.0) 
    // 

    def cleanBMState ( bmState: BMState ) = {
        // Release unnecessary memories except memories for future uses.

        val cleanBMS = BMState() 

      //val cleanInitialMCStates = cleanBMS.initialMCStates
        val cleanInitialMCStates = {
            bmState.initialMCStates.map{ s => s.copy( interactions = MCMC.Interactions() ) }
        } 
        val cleanIndependentSamplings = {
            bmState.independentSamplings.map{ mc => 
                mc.map{ s => s.copy( interactions = MCMC.Interactions() ) }
            }
        } 

        bmState.copy( bmInteractions = cleanBMS.bmInteractions, grad = cleanBMS.grad,
                      v = cleanBMS.v, m = cleanBMS.m, 
                      learningRates = cleanBMS.learningRates,
                      initialMCStates = cleanInitialMCStates,
                      independentSamplings = cleanIndependentSamplings,
                      ensembleAverages = cleanBMS.ensembleAverages,
                      leakyAveP = cleanBMS.leakyAveP, leakyAvePForKL = cleanBMS.leakyAvePForKL )
    }

    def createBMState (gradientDescentMethod: String, bmInteractions: BMInteractions, 
                        grad: BMInteractions, 
                        v: BMInteractions, m: BMInteractions, 
                        betaV: Double, betaM: Double,
                        eps: Double, 
                        learningRate: LearningRate, 
                        step: Int ) = {

                require(betaV >= 0.0 && betaV < 1.0)
                require(betaM >= 0.0 && betaM < 1.0)
                require(step == 0.0)

                val phi: ArraySeq[DenseVector[Double]] = ArraySeq()
                val phij: ArraySeq[DenseMatrix[Double]] = ArraySeq()
                val learningRates = BMInteractions(phi, phij)

                BMState(gradientDescentMethod = gradientDescentMethod,
                        bmInteractions = bmInteractions, grad = grad, 
                        v = v, m = m, 
                        betaV = betaV, betaM = betaM, 
                        eps = eps,
                      //minLearningRate = learningRate.minLR, maxLearningRate = learningRate.maxLR,
                      //rateDecrease = 1.0, rateIncrease = 1.0,
                        learningRate = learningRate, 
                        learningRates = learningRates,
                        step = step )
    }

    // For Adam, ModAdamMax, ...
    def createBMState (gradientDescentMethod: String, bmInteractions: BMInteractions, 
                        grad: BMInteractions, 
                        v: BMInteractions, m: BMInteractions, 
                        betaV: Double, betaM: Double,
                        eps: Double,
                        learningRate: LearningRate, 
                        learningRates: BMInteractions,
                        step: Int ) = {

                require(betaV >= 0.0 && betaV < 1.0)
                require(betaM >= 0.0 && betaM < 1.0)
                require(step == 0.0)

                BMState(gradientDescentMethod = gradientDescentMethod,
                        bmInteractions = bmInteractions, grad = grad, 
                        v = v, m = m, 
                        betaV = betaV, betaM = betaM, 
                        eps = eps,
                        learningRate = learningRate, 
                      //minLearningRate = learningRate.minLR, maxLearningRate = learningRate.maxLR,
                      //rateDecrease = 1.0, rateIncrease = 1.0,
                        learningRates = learningRates,
                        step = step )
    }

    // For RPROP-LR
    def createBMState (gradientDescentMethod: String, bmInteractions: BMInteractions, 
                        grad: BMInteractions, 
                        eps: Double,
                        learningRateForRPROPLR: LearningRateForRPROPLR, 
                      //minLearningRate: Double, maxLearningRate: Double,
                      //rateDecrease: Double, rateIncrease: Double,
                        learningRates: BMInteractions,
                        step: Int ) = {

                require(step == 0.0)

                /*
                v = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
                    ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) )
                        ) 
                */

                val phi: ArraySeq[DenseVector[Double]] = ArraySeq()
                val phij: ArraySeq[DenseMatrix[Double]] = ArraySeq()
                val v = BMInteractions(phi, phij)

                BMState(gradientDescentMethod = gradientDescentMethod,
                        bmInteractions = bmInteractions, grad = grad, 
                        v = v, 
                        m = learningRates,       // obsolete
                        betaV = learningRateForRPROPLR.rateDecrease, betaM = learningRateForRPROPLR.rateIncrease,       // obsolete 
                        eps = eps,
                        learningRateForRPROPLR = learningRateForRPROPLR, 
                      //minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
                      //rateDecrease = rateDecrease, rateIncrease = rateIncrease,
                        learningRates = learningRates,
                        step = step )
    }

    def hJToPhi (hia: ArraySeq[DenseVector[Double]],  jijab: ArraySeq[DenseMatrix[Double]],
                fia: ArraySeq[DenseVector[Double]] ) = {
                
                val Jijab = jijab

                val nUnits = hia.size
                val nPairs = Jijab.size
                val phijab = Jijab
                val phia =  hia.zipWithIndex.map( x => {
                                val hi = x._1
                                val i = x._2
                                //val xj = ArraySeq.range(0, nUnits)
                                //xj.map( j => {
                                ArraySeq.range(0, nUnits).map( j => {
                                        if ( i == j ) {
                                                hi
                                        } else {
                                                val (ij, a, b) = pairIndex(i, j, 0, 1)
                                                  if ( a == 0 )
                                                        Jijab(ij) * fia(j)
                                                  else
                                                        Jijab(ij).t * fia(j)    
                                        }
                                    } ).reduce(_ + _)
                                
                        } )
                (phia, phijab)
    }

    def phiTohJ (phia: ArraySeq[DenseVector[Double]],  phijab: ArraySeq[DenseMatrix[Double]],
                fia: ArraySeq[DenseVector[Double]] ) = {
                
                val nUnits = phia.size
                val nPairs = phijab.size
                val Jijab = phijab
                val hia =  phia.zipWithIndex.map( x => {
                                val phi = x._1
                                val i = x._2
                                val xj = ArraySeq.range(0, nUnits)
                                xj.map( j => {
                                        if ( i == j ) {
                                                phi
                                        } else {
                                                val (ij, a, b) = pairIndex(i, j, 0, 1)
                                                  if ( a == 0 )
                                                        (- phijab(ij) * fia(j))
                                                  else
                                                        (- phijab(ij).t * fia(j))
                                        }
                                    } ).reduce(_ + _)
                                
                        } )
                (hia, Jijab)
    }

    def toIsingGauge( hia: ArraySeq[DenseVector[Double]], jijab: ArraySeq[DenseMatrix[Double]] ) = {
          //pairIndex: (Int, Int, Int, Int) => Tuple3[Int, Int, Int] 

          import breeze.linalg.{InjectNumericOps, sum, `*`}
        //import breeze.linalg.{Axis, BroadcastedColumns, BroadcastedRows}

          val Jijab = jijab

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
          val nPairs = Jijab.size
                
          val indexIJ = ArraySeq.range(0, nPairs).map( ij => inversePairIndex(ij) )   

          val fullJija_ = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                DenseVector.fill(nStatesOfUnit)(0.0)
                                        } else {
                                                val (ij, a, b) = pairIndex(i, j, 0, 1)
                                                val fullJijab = if( a == 0 ) Jijab(ij) else Jijab(ij).t
                                                // fullJijab * fia(j)
                                                breeze.linalg.sum(fullJijab(* , ::) ) / nStatesOfUnit.toDouble
                                                //breeze.linalg.sum(fullJijab, breeze.linalg.Axis._1) / nStatesOfUnit.toDouble
                                        }
                                } )     
                        } )

        // fullJij_b(i)(j)_(b) =  (fia(i).t * fullJijab(i)(j))(b) 
        //                      = ( fia(i).t * fullJijab(j)(i).t )(b) = (fullJijab(j)(i) * fia(i)).t(b) 
        //                      = fullJija_(j)(i).t(b)

          val fullJij_b = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                DenseVector.fill(nStatesOfUnit)(0.0).t
                                        } else {
                                                //fia(i).t * fullJijab(i)(j)
                                                fullJija_(j)(i).t
                                        }
                                } )     
                        } )
          val fullJij__ = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                0.0
                                        } else {
                                                //fia(i).t * fullJija_(i)(j)
                                                breeze.linalg.sum(fullJija_(i)(j)) / nStatesOfUnit
                                        }
                                } )     
                        } )

          val Jijab_Ising = ArraySeq.range(0, nPairs).map( ij =>
                        {
                                val (i, j) = indexIJ(ij)
                                val x = Jijab(ij)(::, *) - fullJija_(i)(j)
                                //val x = breeze.linalg.BroadcastedColumns[DenseMatrix[Double], DenseVector[Double]]( Jijab(ij) ) - fullJija_(i)(j)
                                val y = x.t(::, *) - fullJij_b(i)(j).t
                                //val y = breeze.linalg.BroadcastedColumns[DenseMatrix[Double], DenseVector[Double]] ( x.t ) - fullJij_b(i)(j).t
                                val yt = y.t
                                //val yt = x(*, ::) - fullJij_b(i)(j).t
                                //val yt = breeze.linalg.BroadcastedRows[DenseMatrix[Double], DenseVector[Double]] ( x ) - fullJij_b(i)(j).t
                                yt + fullJij__(i)(j)
                                //(Jijab(ij)(::,*) - fullJija_(i)(j))(*,::) - fullJij_b(i)(j) + fullJij__(i)(j)
                        } )

          val hia_Ising = ArraySeq.range(0, nUnits).map( i =>
                        {
                                //hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - fia(i).t * hia(i)
                                hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - breeze.linalg.sum(hia(i)) / nStatesOfUnit
                        } )

          (hia_Ising, Jijab_Ising)
    }

    def toCompositionGauge( hia: ArraySeq[DenseVector[Double]], jijab: ArraySeq[DenseMatrix[Double]],
                fia: ArraySeq[DenseVector[Double]] ) = {
          //pairIndex: (Int, Int, Int, Int) => Tuple3[Int, Int, Int] 

          import breeze.linalg.{InjectNumericOps, `*`}
        //import breeze.linalg.{BroadcastedColumns, BroadcastedRows}

          val Jijab = jijab

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
          val nPairs = Jijab.size
                
          val indexIJ = ArraySeq.range(0, nPairs).map( ij => inversePairIndex(ij) )   

          val fullJija_ = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                DenseVector.fill(nStatesOfUnit)(0.0)
                                        } else {
                                                val (ij, a, b) = pairIndex(i, j, 0, 1)
                                                val fullJijab = if( a == 0 ) Jijab(ij) else Jijab(ij).t
                                                fullJijab * fia(j)
                                        }
                                } )     
                        } )

        // fullJij_b(i)(j)_(b) =  (fia(i).t * fullJijab(i)(j))(b) 
        //                      = ( fia(i).t * fullJijab(j)(i).t )(b) = (fullJijab(j)(i) * fia(i)).t(b) 
        //                      = fullJija_(j)(i).t(b)

          val fullJij_b = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                DenseVector.fill(nStatesOfUnit)(0.0).t
                                        } else {
                                                //fia(i).t * fullJijab(i)(j)
                                                fullJija_(j)(i).t
                                        }
                                } )     
                        } )
          val fullJij__ = ArraySeq.range(0, nUnits).map( i =>
                        {
                                ArraySeq.range(0, nUnits).map( j =>
                                {
                                        if ( i == j ) {
                                                0.0
                                        } else {
                                                fia(i).t * fullJija_(i)(j)
                                        }
                                } )     
                        } )

          val Jijab_comp = ArraySeq.range(0, nPairs).map( ij =>
                        {
                                val (i, j) = indexIJ(ij)
                                val x = Jijab(ij)(::, *) - fullJija_(i)(j)
                                //val x = breeze.linalg.BroadcastedColumns[DenseMatrix[Double], DenseVector[Double]]( Jijab(ij) ) - fullJija_(i)(j)
                                val y = x.t(::, *) - fullJij_b(i)(j).t
                                //val y = breeze.linalg.BroadcastedColumns[DenseMatrix[Double], DenseVector[Double]] ( x.t ) - fullJij_b(i)(j).t
                                val yt = y.t
                                //val yt = x(*, ::) - fullJij_b(i)(j).t
                                //val yt = breeze.linalg.BroadcastedRows[DenseMatrix[Double], DenseVector[Double]] ( x ) - fullJij_b(i)(j).t
                                yt + fullJij__(i)(j)
                                //(Jijab(ij)(::,*) - fullJija_(i)(j))(*,::) - fullJij_b(i)(j) + fullJij__(i)(j)
                        } )

          val hia_comp = ArraySeq.range(0, nUnits).map( i =>
                        {
                                hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - fia(i).t * hia(i)
                        } )

          (hia_comp, Jijab_comp)
    }

    def toZeroSumGauge( hia: ArraySeq[DenseVector[Double]], jijab: ArraySeq[DenseMatrix[Double]] ) = {

          import breeze.linalg.sum

          val Jijab = jijab

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
          val nPairs = Jijab.size
        
          val nStatesOfUnitSquare = nStatesOfUnit * nStatesOfUnit
          val Jijab_ZeroSum = ArraySeq.range(0, nPairs).map( ij => {
                                val Jij__ = breeze.linalg.sum(Jijab(ij)) / nStatesOfUnitSquare.toDouble
                                Jijab(ij) - Jij__
                              } )
          val hia_ZeroSum = ArraySeq.range(0, nUnits).map( i => {
                                val hi_ = breeze.linalg.sum(hia(i)) / nStatesOfUnit.toDouble
                                hia(i) - hi_
                            } )

          (hia_ZeroSum, Jijab_ZeroSum)
    }

    def toGauge( hia: ArraySeq[DenseVector[Double]], jijab: ArraySeq[DenseMatrix[Double]], gauge: String, 
                 fia: ArraySeq[DenseVector[Double]] ):
                 Tuple2[ArraySeq[DenseVector[Double]], ArraySeq[DenseMatrix[Double]]] = {

            val Jijab = jijab 

            if ( gauge == "unused" || gauge == "ungauged" ) {
                (hia, Jijab)
            } else if ( gauge == "zeroSum" || gauge == "ZeroSum" || gauge == "composition" || gauge == "Ising" ) {
                val (newhia, newJijab) =
                    if ( gauge == "zeroSum" || gauge == "ZeroSum" ) 
                        BM.toZeroSumGauge(hia, Jijab)
                    else if ( gauge == "composition" ) 
                    // caution: h and J are converted to the composition Gauge.
                        BM.toCompositionGauge(hia, Jijab, fia )
                    else
                        BM.toIsingGauge(hia, Jijab)

                (newhia, newJijab)
            } else {
                sys.error("Not supported gauge: " + gauge)
                (hia, Jijab)
            }
    }

    def toGauge( hia: ArraySeq[DenseVector[Double]], jijab: ArraySeq[DenseMatrix[Double]], gauge: String ):
                 Tuple2[ArraySeq[DenseVector[Double]], ArraySeq[DenseMatrix[Double]]] = {
      toGauge( hia, jijab, gauge, fia = ArraySeq(DenseVector[Double]()) ) 
    }

    def toGauge( interactions: MCMC.Interactions, gauge: String, 
                 fia: ArraySeq[DenseVector[Double]] ): MCMC.Interactions = {

        if ( interactions.gauge == gauge ) {
          interactions
        } else {
          val (hia, jijab) = toGauge( interactions.hia, interactions.jijab, gauge, fia)
          MCMC.Interactions(hia, jijab, gauge)
        }
    }

    def toGauge( interactions: MCMC.Interactions, gauge: String ): MCMC.Interactions = {
      toGauge( interactions, gauge, fia = ArraySeq(DenseVector[Double]()) )
    }

    def newPhiThroughhJ(bmInteractions: BMInteractions, fia: ArraySeq[DenseVector[Double]], gauge: String):
                (MCMC.Interactions, BMInteractions) = {

            val phi_ = gauge.slice(0,4)

            if ( gauge == "unused" || gauge == "ungauged" ) {
                val hJ = phiTohJ (bmInteractions.phia, bmInteractions.phijab, fia )

                (MCMC.Interactions(hJ._1, hJ._2, gauge), bmInteractions)
            } else if ( phi_ == "phi_" ) {
                val phi_gauge = gauge.drop(4)
                val newPhiPhij = toGauge(bmInteractions.phia, bmInteractions.phijab, phi_gauge, fia)
                val newhJ = phiTohJ (newPhiPhij._1, newPhiPhij._2, fia )
                val newInteractions = MCMC.Interactions(newhJ._1, newhJ._2, gauge)
                val newBMInteractions = BMInteractions(newPhiPhij._1, newPhiPhij._2, gauge)

                (newInteractions, newBMInteractions)
            } else {
                val hJ_ = gauge.slice(0,3)
                val hJ_gauge = if ( hJ_ == "hJ_" ) gauge.drop(3) else gauge
                val hJ = phiTohJ (bmInteractions.phia, bmInteractions.phijab, fia )
                val newhJ = toGauge( hJ._1, hJ._2, hJ_gauge, fia )
                val newPhiPhij = hJToPhi ( newhJ._1, newhJ._2, fia )
                val newInteractions = MCMC.Interactions(newhJ._1, newhJ._2, gauge)
                val newBMInteractions = BMInteractions(newPhiPhij._1, newPhiPhij._2, gauge)

                (newInteractions, newBMInteractions)
            }
    }
    def newPhiThroughhJ(bmInteractions: BMInteractions, fia: ArraySeq[DenseVector[Double]]):
                (MCMC.Interactions, BMInteractions) = {
                newPhiThroughhJ(bmInteractions, fia, gauge = "unused")
    }

    def newPhiThroughhJ(bmInteractions: BMInteractions,  bmState: BMState, fia: ArraySeq[DenseVector[Double]], gauge: String):
                (MCMC.Interactions, BMState) = {

          //val bmInteractions = bmState.bmInteractions
            val (newInteractions, newBMInteractions) = newPhiThroughhJ(bmInteractions, fia, gauge)

          /*
            val newBMState = BMState( bmInteractions = newBMInteractions, 
                        grad = bmState.grad,
                        v = bmState.v, m = bmState.m,
                        betaV = bmState.betaV, betaM = bmState.betaM, 
                        learningRate = bmState.learningRate,
                        minLearningRate = bmState.minLearningRate, maxLearningRate = bmState.maxLearningRate,
                        rateDecrease = bmState.rateDecrease, rateIncrease = bmState.rateIncrease,
                        learningRates = bmState.learningRates,
                        eps = bmState.eps, step = bmState.step)
          */
            val newBMState = bmState.copy(bmInteractions = newBMInteractions)

            (newInteractions, newBMState)
    }

    def newPhiThroughhJ(bmState: BMState, fia: ArraySeq[DenseVector[Double]], gauge: String):
                (MCMC.Interactions, BMState) = {
            val bmInteractions = bmState.bmInteractions
            newPhiThroughhJ(bmInteractions, bmState, fia, gauge: String)
    }

    def newPhiThroughhJ(bmState: BMState, fia: ArraySeq[DenseVector[Double]]):
                (MCMC.Interactions, BMState) = {
                newPhiThroughhJ(bmState, fia, gauge = "unused")
    }

    def hJ_in_IsingGauge(interactions: MCMC.Interactions) = {
        if ( interactions.gauge == "Ising" || interactions.gauge == "hJ_Ising" ) { 
          interactions
        } else {
          val (hia_Ising, jijab_Ising) = toIsingGauge( interactions.hia, interactions.jijab )
          val interactions_Ising = MCMC.Interactions( hia_Ising, jijab_Ising, gauge = "Ising" )
          interactions_Ising
        }
    }


 // def hJ_and_TE_in_IsingGauge(interactions: MCMC.Interactions, configurations: ArraySeq[IArray[Byte]] ):
 //        Tuple2[MCMC.Interactions, ArraySeq[Double]] = {
 // def hJ_and_TE_in_IsingGauge[T, C <: scala.collection.Seq[T]](interactions: MCMC.Interactions,
 //        configurations: C & scala.collection.SeqOps[T, scala.collection.Seq, C] ):
 //        Tuple2[MCMC.Interactions, ArraySeq[Double]] = {
    def hJ_and_TE_in_IsingGauge[C <: immutable.Seq[IArray[Byte]]] (interactions: MCMC.Interactions,
           configurations: C ) (using classTagC: ClassTag[C]):
           Tuple2[MCMC.Interactions, ArraySeq[Double]] = {

      val interactions_Ising = hJ_in_IsingGauge(interactions)
      val energies_Ising = MCMC.totalE(configurations, interactions_Ising)
      //configurations.map ( configuration => MCMC.totalE(configuration, interactions_Ising) ) 

      (interactions_Ising, energies_Ising)
    }
/*
    def hJ_and_TE_in_IsingGauge(interactions: MCMC.Interactions, configurations: ArraySeq[IArray[Byte]] ):
           Tuple2[MCMC.Interactions, ArraySeq[Double]] = {

      val interactions_Ising = hJ_in_IsingGauge(interactions)
      val energies_Ising = MCMC.totalE(configurations, interactions_Ising)
      //configurations.map ( configuration => MCMC.totalE(configuration, interactions_Ising) ) 

      (interactions_Ising, energies_Ising)
    }
*/
/*
    def hJ_and_TE_in_IsingGauge(interactions: MCMC.Interactions, configurations: ParVector[IArray[Byte]] ):
           Tuple2[MCMC.Interactions, ArraySeq[Double]] = {

      val interactions_Ising = hJ_in_IsingGauge(interactions)
      val energies_Ising = MCMC.totalE(configurations, interactions_Ising)
      //configurations.map ( configuration => MCMC.totalE(configuration, interactions_Ising) ) 

      (interactions_Ising, energies_Ising)
    }
*/
    def hJ_and_TE_in_IsingGauge(interactions: MCMC.Interactions, independentSamplings: MCMC.IndependentSamplings):
           Tuple2[MCMC.Interactions, MCMC.IndependentSamplings] = {

      val interactions_Ising = hJ_in_IsingGauge(interactions)

      val independentSamplings_Ising =
          independentSamplings.map(samples =>
          {
            samples.map( s =>
              {
                val energy_Ising = MCMC.totalE(s.configuration, interactions_Ising) 
              //MCMC.State(s.configuration, energy_Ising, s.kT, s.annealingRate, s.step, s.finalkT)
                s.copy(interactions = interactions_Ising, energy = energy_Ising)
              } )
          } )

      (interactions_Ising, independentSamplings_Ising)
    }

    def leakyAveVec( beta: Double, lAia: ArraySeq[DenseVector[Double]],
                                               ia: ArraySeq[DenseVector[Double]] ) = {
    // To calculate leaky averages of hia and pia

          assert( beta < 1.0 && beta >= 0.0 )

          val newLAia = if ( beta == 0.0 || lAia.size == 0 || lAia(0).size == 0 ) {
              ia
          } else {
              lAia.zip(ia).map{ x => 
                val lA  = x._1
                val current = x._2
                lA * beta + current * (1.0 - beta)
              }
          }
          newLAia
    }

    def leakyAveMat( beta: Double, lAijab: ArraySeq[DenseMatrix[Double]],
                                               ijab: ArraySeq[DenseMatrix[Double]]  ) = {
    // To calculate leaky averages of jijab and pijab

          assert( beta < 1.0 && beta >= 0.0 )

          val newLAijab = if ( beta == 0.0 || lAijab.size == 0 || lAijab(0).size == 0 ) {
              ijab
          } else {
              lAijab.zip(ijab).map{ x => 
                val lA  = x._1
                val current = x._2
                lA * beta + current * (1.0 - beta)
              }
          }
          newLAijab
    }

    // softThresholdingForGL1: L1 for phi, GL1 for phij, because GL1 for phi is meaningless.
    def softThresholdingForGL1(oldBMInteractions: BMInteractions,
                optMethod: String, regTerm: String, 
                propL1h: Double, propL1J: Double, 
                lambdaPhi: Double, lambdaPhij: Double,
                newBMStateWithoutSoft: BMState, learningRates: BMInteractions, softThresholdingOnly: Boolean = true ) = {

          import breeze.linalg.{InjectNumericOps, sum}

        // The increments for interactions must be proportional to learningRates that depend on each parameter.

        // Not completed yet; Only for optMethod == "NAG", L1 corrections are included to new v.

          if(regTerm == "L2" || ( (regTerm == "L1L2" || regTerm == "GL1L2") && ( propL1h == 0.0 && propL1J == 0.0) ) ) {
                newBMStateWithoutSoft
          } else {
            val (pL1h, pL1J) = if ( regTerm == "GL1" ) {
                          assert( propL1h == 1.0 )
                          assert( propL1J == 1.0 )
                          (1.0, 1.0) 
                        } else if ( regTerm == "GL1L2" ) {
                          assert(propL1h >= 0.0 && propL1h <= 1.0)
                          assert(propL1J >= 0.0 && propL1J <= 1.0)
                          (propL1h, propL1J)
                        } else {
                          sys.error(s"Not supported: ${regTerm}")
                          (propL1h, propL1J)
                        }

            val newBMInteractionsWithoutSoft = newBMStateWithoutSoft.bmInteractions

            val nStatesOfUnit = learningRates.phia(0).size

            //val l2NormPhia = oldBMInteractions.phia.map( phi => math.sqrt( breeze.linalg.sum(phi.map(x => x * x)) ) ) 

            val phia = if ( pL1h <= 0.0 ) {
                    newBMInteractionsWithoutSoft.phia
                } else {
                    newBMInteractionsWithoutSoft.phia.zipWithIndex.map( xx => {
                        val (phi, i) = xx
                        val oldPhi = oldBMInteractions.phia(i)
                        val gammai = learningRates.phia(i) * ( lambdaPhi * pL1h )


                        // The contribution of the Regularization term is the same for minimization and maximization. 
                        /*
                        val signi = if (l2NormPhia(i) == 0.0 ) 
                                        DenseVector.zeros[Double](nStatesOfUnit)
                                     else oldPhi / l2NormPhia(i)
                        */
                      /* for L1 */
                        val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
                        val phiWithoutL1 = phi + gammai *:* signi
                      /**/ 
                      /* for GL1
                        val phiWithoutL1 = {
                            val l2NormPhi = math.sqrt( breeze.linalg.sum(oldPhi.map(x => x * x) ) )
                            val signi = if (l2NormPhi == 0.0 )
                                             DenseVector.zeros[Double](nStatesOfUnit)
                                         else oldPhi / l2NormPhi
                            phi + gammai *:* signi
                        }
                      */

                        /* V2: incorrect
                        val phiInSmallCase = phiWithoutL1.mapPairs( (s, p) => { 
                                                if(phiWithoutL1(s).abs <= gammai(s)) 0.0 else 1.0
                                        } )
                        val suffLarge = breeze.linalg.sum( phijInSmallCase )

                        if ( suffLarge == 0 ) {
                          phiInSmallCase
                        } else
                        */
                        {
                            phiWithoutL1.mapPairs( (s, p) => {
                                      /* for L1 */ 
                                        if ( p >= 0.0 ) {
                                          math.max( 0.0, p - gammai(s) )
                                        } else {
                                          //- math.max( 0.0, - p - gammai(s) )
                                          math.min( 0.0, p + gammai(s) )
                                        }
                                      /**/
                                      /* for GL1; refer to CPSC 540
                                        var newP = {
                                               val x = p * math.max( 0.0, l2NormPhi - gammai(s) )
                                               if ( x != 0.0 ) {
                                                   assert( l2NormPhi != 0.0 )
                                                   x / l2NormPhi
                                               } else {
                                                   x
                                               }
                                        }
                                        newP
                                      */
                            } )
                        }
                    } )
                }

          //val l2NormPhijab = oldBMInteractions.phijab.map( phij => math.sqrt( breeze.linalg.sum(phij.map(x => x * x) ) ) ) 
            val phijab = if (pL1J <= 0.0 ) {
                    newBMInteractionsWithoutSoft.phijab
                } else {
                    newBMInteractionsWithoutSoft.phijab.zipWithIndex.map( xx => {
                        val (phij, ij) = xx
                        val oldPhij = oldBMInteractions.phijab(ij)
                        val gammaij = learningRates.phijab(ij) * ( lambdaPhij * pL1J )

                        // The contribution of the Regularization term is the same for minimization and maximization.
                      /* for L1
                        val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
                        val phijWithoutL1 = phij + gammaij *:* signij

                        phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s

                                        // for L1
                                        if ( p >= 0.0 ) {
                                          math.max( 0.0, p - gammaij(s) )
                                        } else {
                                          //- math.max( 0.0, - p - gammaij(s) )
                                          math.min( 0.0, p + gammaij(s) )
                                        }
                                } )
                      */
                      /* for GL1 */
                        val phijWithoutL1 = {
                            val l2NormPhij = math.sqrt( breeze.linalg.sum(oldPhij.map(x => x * x) ) )
                            val signij = if (l2NormPhij == 0.0 ) 
                                             DenseMatrix.zeros[Double](nStatesOfUnit,nStatesOfUnit)
                                         else oldPhij / l2NormPhij
                            phij + gammaij *:* signij
                        }
                      /**/

                        /* V2: incorrect
                        val phijInSmallCase = phijWithoutL1.mapPairs( (s, p) => { 
                                                val (a, b) = s
                                                if(phijWithoutL1(a,b).abs <= gammaij(a,b)) 0.0 else 1.0
                                        } )
                        val suffLarge = breeze.linalg.sum( phijInSmallCase )

                        if ( suffLarge == 0 ) {
                          phijInSmallCase
                        } else
                        */

                        {
                                
                        /*
                          phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s
                                        /*
                                        if ( p > gammaij(a, b) ) p - gammaij(a,b)
                                        else if ( p < - gammaij(a,b) ) p + gammaij(a,b)
                                        else 0.0
                                        */
                                        if ( l2NormPhijab(ij) == 0.0 ) { 
                                                if( p > gammaij(a, b) )
                                                        p - gammaij(a, b)
                                                else if ( p < - gammaij(a, b) )
                                                        p + gammaij(a, b)
                                                else
                                                        0.0
                                        } else {
                                                p - gammaij(a, b) * (oldPhij(a, b) / l2NormPhijab(ij) )
                                        }
                                } )
                        */
                        /* for GL1 */
                          val l2NormPhij = math.sqrt( breeze.linalg.sum( phijWithoutL1.map(x => x * x) ) )
                          phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s
                                        /*
                                        if ( p > gammaij(a, b) ) p - gammaij(a,b)
                                        else if ( p < - gammaij(a,b) ) p + gammaij(a,b)
                                        else 0.0
                                        */
                                        /* Error
                                        val dp = if ( l2NormPhijab(ij) == 0.0 ) {
                                                   if ( p - gammaij(a, b) > 0.0 )
                                                        p - gammaij(a, b)
                                                   else if ( p + gammaij(a, b) < 0.0 )
                                                         p + gammaij(a, b)
                                                   else
                                                        0.0

                                                } else {
                                                        p - gammaij(a, b) * (oldPhij(a, b) / l2NormPhijab(ij) )
                                                }

                                        if ( oldPhij(a, b) > 0.0 ) {
                                                if ( dp > 0.0 )
                                                        dp
                                                else
                                                        0.0
                                        } else if ( oldPhij(a, b) < 0.0 ) {
                                                if ( dp < 0.0 )
                                                        dp
                                                else
                                                        0.0
                                        } else {
                                                dp
                                        }
                                        */

                                        // refer to CPSC 540
                                        var newP = {
                                               val x = p * math.max( 0.0, l2NormPhij - gammaij(s) ) 
                                               if ( x != 0.0 ) {
                                                   assert( l2NormPhij != 0.0 )
                                                   x / l2NormPhij 
                                               } else {
                                                   x
                                               }
                                        }
                                        newP
                          } )
                        /**/
                        }

                    } )
                }
                                        
            val newBMInteractionsWithL1 = BMInteractions(phia, phijab)

            val newV = if ( softThresholdingOnly ) {
                newBMStateWithoutSoft.v
            } else {
                sys.error("*** Corrections on BMState.v and BMState.m due to soft-thresholding are not completed; probably they are not needed.\n")
            } 
          /*
            } else if ( optMethod == "NAG" ) {
                /* not completed; probably the following correction is not needed. */

                val step = newBMStateWithoutSoft.step - newBMStateWithoutSoft.stepOffset  // newBMStateWithoutSoft.step was already incremented.
                assert( step >= 1 )
                val learnRate = newBMStateWithoutSoft.learningRate.lr(step)               // step was already incremented.  
                val vPhia = if (pL1h <= 0.0) {
                    newBMStateWithoutSoft.v.phia
                } else {
                    newBMStateWithoutSoft.v.phia.zipWithIndex.map( xx => {
                        val (v, i) = xx
                        val oldPhi = oldBMInteractions.phia(i)

                        val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

                        //val vWithoutL1 = v +  learnRate * lambdaPhi * pL1h * signi
                        val vWithoutL1 = v +  signi * (learnRate * lambdaPhi * pL1h)

                        val gammai = learningRates.phia(i) * ( lambdaPhi * pL1h )
                        val phi = newBMInteractionsWithoutSoft.phia(i)
                        val phiWithoutL1 = phi + gammai *:* signi

                        val vL1 = phiWithoutL1.mapPairs( (s, p) => {
                                        if ( p > gammai(s) ) {
                                                - gammai(s)
                                        } else if ( p < - gammai(s) ) {
                                                gammai(s)
                                        } else {
                                                - p
                                        }
                                } )
                        vWithoutL1 + vL1
                        } )
                }

                val vPhijab = if (pL1J <= 0.0) {
                    newBMStateWithoutSoft.v.phijab
                } else {
                    newBMStateWithoutSoft.v.phijab.zipWithIndex.map( xx => {
                        val (v, ij) = xx
                        val oldPhij = oldBMInteractions.phijab(ij)

                        val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

                        //val vWithoutL1 = v +  learnRate * lambdaPhij * pL1J * signij
                        val vWithoutL1 = v +  signij * ( learnRate * lambdaPhij * pL1J )

                        val gammaij = learningRates.phijab(ij) * ( lambdaPhij * pL1J )
                        val phij = newBMInteractionsWithoutSoft.phijab(ij)
                        val phijWithoutL1 = phij + gammaij *:* signij

                        val vL1 = phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s
                                        if ( p > gammaij(a, b) ) {
                                                - gammaij(a,b)
                                        } else if ( p < - gammaij(a,b) ) {
                                                gammaij(a,b)
                                        } else {
                                                - p
                                        }
                                } )
                        vWithoutL1 + vL1
                        } )
                }

                val vWithSoft = BMInteractions(vPhia, vPhijab) 
                vWithSoft
            } else { 
                sys.error("*** Corrections on BMState.v and BMState.m due to soft-thresholding are not completed; probably they are not needed.\n")

                /* The correction for v is ignored for other methods */
                newBMStateWithoutSoft.v
            }
          */
          //BMState(
            newBMStateWithoutSoft.copy(
                bmInteractions = newBMInteractionsWithL1, 
                grad = newBMStateWithoutSoft.grad,
                v = newV,
                m = newBMStateWithoutSoft.m, 
                betaV = newBMStateWithoutSoft.betaV, betaM = newBMStateWithoutSoft.betaM,
                learningRate = newBMStateWithoutSoft.learningRate,
              //minLearningRate = newBMStateWithoutSoft.minLearningRate, maxLearningRate = newBMStateWithoutSoft.maxLearningRate,
              //rateDecrease = newBMStateWithoutSoft.rateDecrease, rateIncrease = newBMStateWithoutSoft.rateIncrease,
                learningRates = newBMStateWithoutSoft.learningRates,
                eps = newBMStateWithoutSoft.eps, step = newBMStateWithoutSoft.step )

          }

    }

    def softThresholdingForL1(oldBMInteractions: BMInteractions,
                optMethod: String, regTerm: String,
                propL1h: Double, propL1J: Double,
                lambdaPhi: Double, lambdaPhij: Double,
                newBMStateWithoutSoft: BMState, learningRates: BMInteractions, softThresholdingOnly: Boolean = true ) = {

          import breeze.linalg.{InjectNumericOps, sum}

        // The increments for interactions must be proportional to learningRates that depend on each parameter.

        // Not completed yet; Only for optMethod == "NAG", L1 corrections are included to new v.

          if(regTerm == "L2" || ( (regTerm == "L1L2" || regTerm == "GL1L2") && (propL1h == 0.0 && propL1J == 0.0) ) ) {
                newBMStateWithoutSoft
          } else {
            val (pL1h, pL1J) = if ( regTerm == "L1" ) {
                          assert( propL1h == 1.0 )
                          assert( propL1J == 1.0 )
                          (1.0, 1.0) 
                        } else if ( regTerm == "L1L2" ) {
                          assert(propL1h >= 0.0 && propL1h <= 1.0)
                          assert(propL1J >= 0.0 && propL1J <= 1.0)
                          (propL1h, propL1J)
                        } else {
                          sys.error(s"Not supported: ${regTerm}")
                          (propL1h, propL1J)
                        }

            val newBMInteractionsWithoutSoft = newBMStateWithoutSoft.bmInteractions

            val nStatesOfUnit = learningRates.phia(0).size

            val phia = if (pL1h <= 0.0 ) {
                    newBMInteractionsWithoutSoft.phia
                } else {
                    newBMInteractionsWithoutSoft.phia.zipWithIndex.map( xx => {
                        val (phi, i) = xx
                        val oldPhi = oldBMInteractions.phia(i)
                        val gammai = learningRates.phia(i) * ( lambdaPhi * pL1h )

                        // The contribution of the Regularization term is the same for minimization and maximization.
                        val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
                        val phiWithoutL1 = phi + gammai *:* signi

                        phiWithoutL1.mapPairs( (s, p) => {

                                        // for L1
                                        if ( p >= 0.0 ) {
                                          math.max( 0.0, p - gammai(s) )
                                        } else {
                                          //- math.max( 0.0, - p - gammai(s) )
                                          math.min( 0.0, p + gammai(s) )
                                        }

                                      /*
                                        if ( p > gammai(s) ) {
                                                p - gammai(s)
                                        } else if ( p < - gammai(s) ) {
                                                p + gammai(s)
                                        } else {
                                                0.0
                                        }
                                      */
                                } )
                        } )
                }

            val phijab = if (pL1J <= 0.0 ) {
                    newBMInteractionsWithoutSoft.phijab
                } else {
                    newBMInteractionsWithoutSoft.phijab.zipWithIndex.map( xx => {
                        val (phij, ij) = xx
                        val oldPhij = oldBMInteractions.phijab(ij)
                        val gammaij = learningRates.phijab(ij) * ( lambdaPhij * pL1J )

                        // The contribution of the Regularization term is the same for minimization and maximization.
                        val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
                        val phijWithoutL1 = phij + gammaij *:* signij

                        phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s

                                        // for L1
                                        if ( p >= 0.0 ) {
                                          math.max( 0.0, p - gammaij(s) )
                                        } else {
                                          //- math.max( 0.0, - p - gammaij(s) )
                                          math.min( 0.0, p + gammaij(s) )
                                        }

                                      /*
                                        if ( p > gammaij(a, b) ) {
                                                p - gammaij(a,b)
                                        } else if ( p < - gammaij(a,b) ) {
                                                p + gammaij(a,b)
                                        } else {
                                                0.0
                                        }
                                      */
                                } )
                        } )
                }
                                        
            val newBMInteractionsWithL1 = BMInteractions(phia, phijab)

            val newV = if ( softThresholdingOnly ) {

                newBMStateWithoutSoft.v

              } else if ( optMethod == "NAG" ) {
                 /* Probably the following correction is not needed. */

                stderr.print("*** Corrections on BMState.v and BMState.m due to soft-thresholding are probably not needed.\n")

                val step = newBMStateWithoutSoft.step - newBMStateWithoutSoft.stepOffset  // newBMStateWithoutSoft.step was already incremented.
                assert( step >= 1 )
                val learnRate = newBMStateWithoutSoft.learningRate.lr( step )

                val vPhia = if (pL1h <= 0.0 ) {
                    newBMStateWithoutSoft.v.phia
                  } else {
                    newBMStateWithoutSoft.v.phia.zipWithIndex.map( xx => {
                        val (v, i) = xx
                        val oldPhi = oldBMInteractions.phia(i)

                        val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

                        //val vWithoutL1 = v +  ( learnRate * lambdaPhi * pL1h ) * signi
                        val vWithoutL1 = v +  signi * ( learnRate * lambdaPhi * pL1h )

                        val gammai = learningRates.phia(i) * ( lambdaPhi * pL1h )
                        val phi = newBMInteractionsWithoutSoft.phia(i)
                        val phiWithoutL1 = phi + gammai *:* signi

                        val vL1 = phiWithoutL1.mapPairs( (s, p) => {
                                        if ( p > gammai(s) ) {
                                                - gammai(s)
                                        } else if ( p < - gammai(s) ) {
                                                gammai(s)
                                        } else {
                                                - p
                                        }
                                } )
                        vWithoutL1 + vL1
                        } )
                  }

                val vPhijab = if (pL1J <= 0.0) {
                    newBMStateWithoutSoft.v.phijab
                  } else {
                    newBMStateWithoutSoft.v.phijab.zipWithIndex.map( xx => {
                        val (v, ij) = xx
                        val oldPhij = oldBMInteractions.phijab(ij)

                        val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

                        //val vWithoutL1 = v +  ( learnRate * lambdaPhij * pL1J ) * signij
                        val vWithoutL1 = v +  signij * ( learnRate * lambdaPhij * pL1J )

                        val gammaij = learningRates.phijab(ij) * ( lambdaPhij * pL1J )
                        val phij = newBMInteractionsWithoutSoft.phijab(ij)
                        val phijWithoutL1 = phij + gammaij *:* signij

                        val vL1 = phijWithoutL1.mapPairs( (s, p) => {
                                        val (a, b) = s
                                        if ( p > gammaij(a, b) ) {
                                                - gammaij(a,b)
                                        } else if ( p < - gammaij(a,b) ) {
                                                gammaij(a,b)
                                        } else {
                                                - p
                                        }
                                } )
                        vWithoutL1 + vL1
                        } )
                  }

                val vWithSoft = BMInteractions(vPhia, vPhijab) 
                vWithSoft
              } else { 
                sys.error("*** Corrections on BMState.v and BMState.m due to soft-thresholding are not completed; probably they are not needed.\n")

                /* The correction for v is ignored for other methods */
                newBMStateWithoutSoft.v
              }

          //BMState(

            newBMStateWithoutSoft.copy(
                bmInteractions = newBMInteractionsWithL1, 
                grad = newBMStateWithoutSoft.grad,
                v = newV,
                m = newBMStateWithoutSoft.m, 
                betaV = newBMStateWithoutSoft.betaV, betaM = newBMStateWithoutSoft.betaM,
                learningRate = newBMStateWithoutSoft.learningRate,
              //minLearningRate = newBMStateWithoutSoft.minLearningRate, maxLearningRate = newBMStateWithoutSoft.maxLearningRate,
              //rateDecrease = newBMStateWithoutSoft.rateDecrease, rateIncrease = newBMStateWithoutSoft.rateIncrease,
                learningRates = newBMStateWithoutSoft.learningRates,
                eps = newBMStateWithoutSoft.eps, step = newBMStateWithoutSoft.step)

          }

    }

    // GL1 means group L1 for J; for propL1= 1, R = sum_i sum_a |hia|  + sum_ij sqrt sum_ab Jijab^2 
    // GL1 for phi seems to be meaningless.
    def dLdphiWithGL1L2(propL1h: Double, propL1J: Double,
                fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]],
                pia: ArraySeq[DenseVector[Double]], pijab: ArraySeq[DenseMatrix[Double]],
                phia: ArraySeq[DenseVector[Double]],  phijab: ArraySeq[DenseMatrix[Double]], 
                //hia: ArraySeq[DenseVector[Double]],  jijab: ArraySeq[DenseMatrix[Double]], 
                lambdaPhi: Double, lambdaPhij: Double  ) = {

                import breeze.linalg.{InjectNumericOps, sum}

                //val (phia, phijab) = hJToPhi (hia,  jijab, fia)

                // Please notice that softThresholding will be taken care of later. 

                assert( propL1h >= 0.0 && propL1h <= 1.0 )
                assert( propL1J >= 0.0 && propL1J<= 1.0 )
                val nUnits = fia.size
                val nPairs = fijab.size
                val propL2h = 1.0 - propL1h
                val propL2J = 1.0 - propL1J

                //val l2NormPhia = phia.map( phi => math.sqrt( breeze.linalg.sum(phi.map(phia => phia * phia)) ) ) 
                val dldphia = phia.zipWithIndex.map( x => {
                                val (phi, i) = x
                                val pi = pia(i)
                                val fi = fia(i)
                                if ( propL1h == 0.0 ) {
                                  //fi - pi - lambdaPhi * ( propL2h * phi )
                                  fi - pi - phi * (propL2h * lambdaPhi)
                                } else {
                                /* for L1 */
                                  val signPhi = phi.map( x => {
                                        if (x > 0.0) 1.0 
                                        else if( x < 0.0) -1.0 
                                        else { 0.0 }
                                        } )

                                /* for GL1
                                  /*
                                  val signPhi = phi.map( x => {
                                        if ( l2NormPhia(i) == 0.0 ) {
                                                0.0
                                        } else {
                                                x / l2NormPhia(i)
                                        }
                                  } )
                                  */

                                  val signPhi = if ( l2NormPhia(i) == 0.0 ) {
                                                   phia(i).map{ x => 0.0 }
                                               } else {
                                                   phia(i) / l2NormPhia(i)
                                               }
                                */
                                                
                                  if (propL2h == 0.0 )
                                    //fi - pi - lambdaPhi * ( propL1h * signPhi )
                                    fi - pi - signPhi * ( propL1h * lambdaPhi )
                                  else
                                    //fi - pi - lambdaPhi * ( propL1h * signPhi + propL2h * phi )
                                    fi - pi - ( signPhi * ( propL1h * lambdaPhi ) + phi * ( propL2h * lambdaPhi ) )
                                }
                        } )

                val l2NormPhijab = phijab.map( phij => math.sqrt( breeze.linalg.sum(phij.map(x => x * x) ) ) ) 

                val dldphijab = ArraySeq.range(0, nPairs).map { ij =>
                          val (i, j) = inversePairIndex(ij)
                          val (ij_ , a, b) = pairIndex(i, j, 0, 1)
                          assert( ij == ij_ )
                        //dldphijab(ij) =
                          if ( propL1J == 0.0 ) {
                                if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - phijab(ij) * ( propL2J * lambdaPhij )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * ( propL2J * phijab(ij) )
                                } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - phijab(ij) * ( propL2J * lambdaPhij )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * ( propL2J * phijab(ij) )
                                }
                          } else {
                              /* for L1
                                val signPhij = phijab(ij).map( x => {
                                        if (x > 0.0) 1.0 
                                        else if( x < 0.0) -1.0 
                                        else { 0.0 }
                                        } )
                              */
                              /* for GL1 */
                               /*
                                val signPhij = phijab(ij).map( x => {
                                        if ( l2NormPhijab(ij) == 0.0 ) {
                                                0.0
                                        } else {
                                                x / l2NormPhijab(ij)
                                        }
                                } )
                               */
                                val signPhij = if ( l2NormPhijab(ij) == 0.0 ) {
                                                   phijab(ij).map{ x => 0.0 } 
                                               } else {
                                                   phijab(ij) / l2NormPhijab(ij)
                                               }
                               /**/

                                if ( propL2J == 0.0 ) 
                                  if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - signPhij * ( propL1J * lambdaPhij )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij )
                                  } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - signPhij * ( propL1J * lambdaPhij )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij )
                                  }
                                else
                                  if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - ( signPhij * (propL1J * lambdaPhij) + phijab(ij) * (propL2J * lambdaPhij) )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
                                  } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - ( signPhij * (propL1J * lambdaPhij) + phijab(ij) * (propL2J * lambdaPhij) )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
                                  }
                          }
                  }
                (dldphia, dldphijab)
    }

    def dLdphiWithL1L2(propL1h: Double, propL1J: Double,
                fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]],
                pia: ArraySeq[DenseVector[Double]], pijab: ArraySeq[DenseMatrix[Double]],
                phia: ArraySeq[DenseVector[Double]],  phijab: ArraySeq[DenseMatrix[Double]], 
                //hia: ArraySeq[DenseVector[Double]],  jijab: ArraySeq[DenseMatrix[Double]], 
                lambdaPhi: Double, lambdaPhij: Double  ) = {

                import breeze.linalg.InjectNumericOps

                //val (phia, phijab) = hJToPhi (hia,  jijab, fia)

                // Please notice that softThresholding will be taken care of later. 

                assert( propL1h >= 0.0 && propL1h <= 1.0 )
                assert( propL1J >= 0.0 && propL1J <= 1.0 )
                val nUnits = fia.size
                val nPairs = fijab.size
                val propL2h = 1.0 - propL1h
                val propL2J = 1.0 - propL1J

                val dldphia = phia.zip(pia).zip(fia).map( x => {
                                val phi = x._1._1
                                val pi = x._1._2
                                val fi = x._2
                                if ( propL1h == 0.0 ) {
                                  //fi - pi - lambdaPhi * ( propL2h * phi )
                                  fi - pi - phi * (propL2h * lambdaPhi)
                                } else {
                                  val signPhi = phi.map( x => {
                                        if (x > 0.0) 1.0 
                                        else if( x < 0.0) -1.0 
                                        else { 0.0 }
                                        } )
                                  if ( propL2h == 0.0 )
                                    //fi - pi - lambdaPhi * ( propL1h * signPhi )
                                    fi - pi - signPhi * ( propL1h * lambdaPhi )
                                  else
                                    //fi - pi - lambdaPhi * ( propL1h * signPhi + propL2h * phi )
                                    fi - pi - ( signPhi * (propL1h * lambdaPhi) + phi * (propL2h * lambdaPhi) )
                                }
                        } )

                val dldphijab = ArraySeq.range(0, nPairs).map { ij =>
                          val (i, j) = inversePairIndex(ij)
                          val (ij_ , a, b) = pairIndex(i, j, 0, 1)
                          assert( ij == ij_ )
                        //dldphijab(ij) = 
                          if ( propL1J == 0.0 ) {
                                if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - phijab(ij) * ( propL2J * lambdaPhij )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * ( propL2J * phijab(ij) )
                                } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - phijab(ij) * ( propL2J * lambdaPhij )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * ( propL2J * phijab(ij) )
                                }
                          } else {
                            val signPhij = phijab(ij).map( x => {
                                        if (x > 0.0) 1.0 
                                        else if( x < 0.0) -1.0 
                                        else { 0.0 }
                                        } )
                                if ( propL2J == 0.0 )
                                  if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - signPhij * ( propL1J * lambdaPhij )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij )
                                  } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - signPhij * ( propL1J * lambdaPhij )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij )
                                  }
                                else
                                  if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - ( signPhij * (propL1J * lambdaPhij) + phijab(ij) * (propL2J * lambdaPhij) )
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
                                  } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - ( signPhij * (propL1J * lambdaPhij) + phijab(ij) * (propL2J * lambdaPhij) )
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
                                  }
                          }
                }

                (dldphia, dldphijab)
    }

    def dLdphiWithL1(fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]],
                pia: ArraySeq[DenseVector[Double]], pijab: ArraySeq[DenseMatrix[Double]],
                phia: ArraySeq[DenseVector[Double]],  phijab: ArraySeq[DenseMatrix[Double]], 
                //hia: ArraySeq[DenseVector[Double]],  jijab: ArraySeq[DenseMatrix[Double]], 
                lambdaPhi: Double, lambdaPhij: Double ) = {

                import breeze.linalg.InjectNumericOps

                //val (phia, phijab) = hJToPhi (hia,  jijab, fia)
                
                // Please notice that softThresholding will be taken care of later. 

                val nUnits = fia.size
                val nPairs = fijab.size

                val dldphia = phia.zip(pia).zip(fia).map( x => {
                                        val phi = x._1._1
                                        val pi = x._1._2
                                        val fi = x._2
                                        val signPhi = phi.map( x => if (x > 0.0) 1.0 else if( x < 0.0) -1.0 else 0.0)
                                        //fi - pi - lambdaPhi * signPhi 
                                        fi - pi - signPhi * lambdaPhi
                                } )

                val dldphijab = ArraySeq.range(0, nPairs).map { ij =>
                        val (i, j) = inversePairIndex(ij)
                        val (ij_ , a, b) = pairIndex(i, j, 0, 1)
                        assert( ij == ij_ )
                        val signPhij = phijab(ij).map( x => if (x > 0.0) 1.0 else if( x < 0.0) -1.0 else 0.0)
                      //dldphijab(ij) =
                        if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - signPhij * lambdaPhij
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * signPhij
                        } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
                                                (pia(j) - fia(j)) * fia(i).t - signPhij * lambdaPhij
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * signPhij
                        }
                }
                (dldphia, dldphijab)
    }

    def dLdphiWithL2(fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]],
                pia: ArraySeq[DenseVector[Double]], pijab: ArraySeq[DenseMatrix[Double]],
                phia: ArraySeq[DenseVector[Double]],  phijab: ArraySeq[DenseMatrix[Double]], 
                //hia: ArraySeq[DenseVector[Double]],  jijab: ArraySeq[DenseMatrix[Double]], 
                lambdaPhi: Double, lambdaPhij: Double  ) = {

                import breeze.linalg.InjectNumericOps

                //val (phia, phijab) = hJToPhi (hia,  jijab, fia)
                
                val nUnits = fia.size
                val nPairs = fijab.size

                val dldphia = phia.zip(pia).zip(fia).map( x => {
                                        val phi = x._1._1
                                        val pi = x._1._2
                                        val fi = x._2
                                        //fi - pi - lambdaPhi * phi  
                                        fi - pi - phi * lambdaPhi
                                } )

                /*
                val cijab = new Array[DenseMatrix[Double]](nPairs)
                val gammaijab = new Array[DenseMatrix[Double]](nPairs)
                val deltaia = pia.zip(fia).map( x => x._1 - x._2 )
                cfor(0)( i => i < nUnits, i => i + 1)( i => {
                  cfor(0)( j => j < i, j => j + 1)( j => {
                                val (ij, a, b) = pairIndex(i, j, 0, 1)
                                cijab(ij) = if ( a == 0 ) {
                                                fijab(ij) - fia(i) * fia(j).t
                                        } else {
                                                fijab(ij) - fia(j) * fia(i).t
                                        }
                                gammaijab(ij) = if ( a == 0 ) {
                                                pijab(ij) - pia(i) * pia(j).t
                                        } else {
                                                pijab(ij) - pia(j) * pia(i).t
                                        }
                                
                  } )
                } )
                */

                val dldphijab = ArraySeq.range(0, nPairs).map { ij =>
                        val (i, j) = inversePairIndex(ij)
                        val (ij_ , a, b) = pairIndex(i, j, 0, 1)
                        assert( ij == ij_ )
                      //dldphijab(ij) =
                        if ( a == 0 ) {
                                        fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
                                                fia(i) * (pia(j) - fia(j)).t - phijab(ij) * lambdaPhij
                                        //      fia(i) * (pia(j) - fia(j)).t - lambdaPhij * phijab(ij)
                        } else {
                                        fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t + 
                                                (pia(j) - fia(j)) * fia(i).t - phijab(ij) * lambdaPhij
                                        //      (pia(j) - fia(j)) * fia(i).t - lambdaPhij * phijab(ij)
                        }
                }
                (dldphia, dldphijab)
    }

    //val (dldphia, dldphijab) = dLdphi(fia, fijab, pia, pijab, hia, Jijab, lambda )

    def learningByNAG(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ): 
                                                        Tuple2[BMState, Option[BMInteractions]] = {

      // v(t) = beta x(t-1) + ( 1 - beta ) rate grad 
      // x(t) = x(t-1) + minOrMax ( beta normalized-v(t) + (1 - beta) rate grad )
      // Rate = ( 1 - beta ) rate

          import breeze.linalg.{InjectNumericOps, max}

          assert(minOrMax.abs == 1)
          assert(0.0 <= bmState.betaV && bmState.betaV < 1.0)

          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )

        // Be careful of the definition of "leaningRate". 
          val learningRate = bmState.learningRate.lr( step ) * ( 1.0 - bmState.betaV )

        /* The following code was used until 191016
          val betaV = if ( step <= 10 ) {               // betaV == mu
                        val r = 1.0 - math.pow(0.5, step - 1)
                        if ( r > bmState.betaV ) bmState.betaV else r
                } else {
                        bmState.betaV
                }
        */
          val betaV = bmState.betaV
        // v for step = 0 is assumed to be 0.
        /*
          val normalizeBetaVPrevStep = if ( step > 1 ) {
                                        1.0 / (1.0 - math.pow(betaV, step - 1) )
                                } else {
                                        0.0
                                }
          val normalizeBetaVNextStep = 1.0 / (1.0 - math.pow(betaV, step + 1) )
        */
          val normalizeBetaVStep = 1.0 / (1.0 - math.pow(betaV, step) )

        //stderr.println( bmState.v.phia.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        betaV * v + learningRate * g
                } ) 
        //stderr.println( vphia.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
        //stderr.println( bmState.v.phijab.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        betaV * v + learningRate * g
                } ) 
        //stderr.println( vphijab.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )

          val phia = ArraySeq.range(0, vphia.size).map( i => {
                        val s = bmState.bmInteractions.phia(i)
                        s + ( vphia(i) * (betaV * normalizeBetaVStep * minOrMax.toDouble) + grad.phia(i) * (learningRate * minOrMax.toDouble ) )
                        //s + minOrMax.toDouble *
                        //      ( betaV * vphia(i) * normalizeBetaVStep + learningRate * grad.phia(i) )
                        /*
                        s + minOrMax.toDouble *
                                ( betaV * betaV * bmState.v.phia(i) * normalizeBetaVPrevStep + ( 1.0 + betaV ) * learningRate * grad.phia(i) )
                        */
                        /* BM.scala.NAG.revised.191024
                        s + minOrMax.toDouble *
                                ( normalizeBetaVStep + betaV * normalizeBetaVNextStep ) * vphia(i) - betaV * normalizeBetaVStep * bmState.v.phia(i)
                        */
                        /* BM.scala.NAG.revised.191017 
                        s + minOrMax.toDouble * 
                            ( (1.0 + betaV ) * vphia(i) * normalizeBetaVStep - betaV * bmState.v.phia(i) * normalizeBetaVPrevStep )   
                        */
                        /* BM.scala.NAG.revised.191017
                            (vphia(i) * normalizeBetaVStep + betaV * 
                                ( vphia(i) * normalizeBetaVStep - bmState.v.phia(i) * normalizeBetaVPrevStep) )   
                        */
                } ) 
          val phijab = ArraySeq.range(0, vphijab.size).map( ij => {
                        val s = bmState.bmInteractions.phijab(ij)
                        s + ( vphijab(ij) * (betaV * normalizeBetaVStep * minOrMax.toDouble) + grad.phijab(ij) * (learningRate * minOrMax.toDouble) )
                        //s + minOrMax.toDouble *
                        //      ( betaV * vphijab(ij) * normalizeBetaVStep + learningRate * grad.phijab(ij) )
                        /*
                        s + minOrMax.toDouble *
                                ( betaV * betaV * bmState.v.phijab(ij) * normalizeBetaVPrevStep + ( 1.0 + betaV ) * learningRate * grad.phijab(ij) )
                        */
                        /* BM.scala.NAG.revised.191024
                        s + minOrMax.toDouble *
                                ( normalizeBetaVStep + betaV * normalizeBetaVNextStep ) * vphijab(ij) - betaV * normalizeBetaVStep * bmState.v.phijab(ij)
                        */
                        /* BM.scala.NAG.revised.191017
                        s + minOrMax.toDouble * 
                            ( (1.0 + betaV ) * vphijab(ij) * normalizeBetaVStep - betaV * bmState.v.phijab(ij) * normalizeBetaVPrevStep )   
                        */
                        /* BM.scala.NAG.revised.191017
                            (vphijab(ij) * normalizeBetaVStep + betaV * 
                                ( vphijab(ij) * normalizeBetaVStep - bmState.v.phijab(ij) * normalizeBetaVPrevStep ) )   
                        */
                } ) 

          val nStatesOfUnit = phia(0).size
          /* BM.scala.NAG.revised.191017
          val lR = learningRate * (1.0 + betaV) * normalizeBetaVStep
          */
          /* BM.scala.NAG.revised.191024
          val lR = learningRate * (normalizeBetaVStep + betaV * normalizeBetaVNextStep)
          */
          /*
          val lR = learningRate * (1.0 + betaV)
          */
          val lR = learningRate * (betaV * normalizeBetaVStep + 1.0)

          val learningRatePhia = phia.map( phi => {
                        DenseVector.fill(nStatesOfUnit)( lR )
                        } )
          val learningRatePhijab = phijab.map( phij => {
                        DenseMatrix.fill(nStatesOfUnit, nStatesOfUnit)( lR )
                        } )

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)

          ( bmState.copy( 
                        bmInteractions = BMInteractions(phia, phijab),
                        grad = grad,
                        v = BMInteractions(vphia, vphijab), m = bmState.m,
                        betaV = bmState.betaV, betaM = bmState.betaM,
                        learningRate = bmState.learningRate,
                        learningRates = learningRates,
                        eps = bmState.eps, step = bmState.step + 1 ),
                        Option(learningRates) )

    }

    def learningByAdadelta(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {
        // reference: Dive into Deep Learning

          import breeze.linalg.InjectNumericOps
          import breeze.numerics.sqrt

          assert(minOrMax.abs == 1)
          val eps = bmState.eps
          val beta = bmState.betaV
          assert( bmState.betaM == bmState.betaV )
          val betaV = beta
          val betaM = beta

          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )

        /* 2024/04/15: Errors was corrected; learningRatePhi was incorrectly defined. */

          val nUnits = bmState.bmInteractions.phia.size
          val nPairs = (nUnits * ( nUnits - 1)) / 2 
          val nStatesOfUnit = bmState.bmInteractions.phia(0).size

       // bmState.v and bmState.m must be initialized in such a way that all elements are equal to zero.

          val g2phia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val g2phijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 

        /* error
        // Errors was corrected; learningRatePhi was incorrectly defined.
          val learningRatePhia = 
                if( step <= 1) {
                        ArraySeq.range(0, nUnits).map( _ => {
                                DenseVector.fill(nStatesOfUnit)(bmState.learningRate.lr(step))
                        } )
                } else {
                        bmState.m.phia.zip(bmState.v.phia).map( x => { 
                                val m = x._1
                                val v = x._2
                                breeze.numerics.sqrt(m + eps) /:/
                                  breeze.numerics.sqrt(v + eps)
                                } ) 
                }

          val learningRatePhijab =
                if( step <= 1) {
                        ArraySeq.range(0, nPairs).map( _ => {
                                DenseMatrix.fill(nStatesOfUnit, nStatesOfUnit)(bmState.learningRate.lr(step))
                        } )
                } else {
                        bmState.m.phijab.zip(bmState.v.phijab).map( x => { 
                                val m = x._1
                                val v = x._2
                                breeze.numerics.sqrt(m + eps) /:/
                                  breeze.numerics.sqrt(v + eps)
                                } ) 
                }
        */

          val learningRatePhia = {
                        bmState.m.phia.zip(g2phia).map( x => { 
                                val m = x._1
                                val v = x._2
                                breeze.numerics.sqrt(m + eps) /:/
                                  breeze.numerics.sqrt(v + eps)
                                } ) 
          }

          val learningRatePhijab = {
                        bmState.m.phijab.zip(g2phijab).map( x => { 
                                val m = x._1
                                val v = x._2
                                breeze.numerics.sqrt(m + eps) /:/
                                  breeze.numerics.sqrt(v + eps)
                                } ) 
          }

        /*
        */

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)

          val dphia = learningRatePhia.zip(grad.phia).map( x => {
                        val r = x._1
                        val g = x._2
                        minOrMax.toDouble * (r *:* g)
                        } )

          val dphijab = learningRatePhijab.zip(grad.phijab).map( x => {
                        val r = x._1
                        val g = x._2
                        minOrMax.toDouble * (r *:* g)
                        } )

          val phia = bmState.bmInteractions.phia.zip(dphia).map( x => {
                                x._1 + x._2
                        } )
          val phijab = bmState.bmInteractions.phijab.zip(dphijab).map( x => {
                                x._1 + x._2
                        } )

          val dphia2 = bmState.m.phia.zip(dphia).map( dx => {
                        val m = dx._1
                        val d = dx._2
                        //betaM * m + (1.0 - betaM) * (d *:* d)
                        m * betaM + (d *:* d) * (1.0 - betaM)
                } ) 
          val dphijab2 = bmState.m.phijab.zip(dphijab).map( dx => {
                        val m = dx._1
                        val d = dx._2
                        //betaM * m + (1.0 - betaM) * (d *:* d)
                        m * betaM + (d *:* d) * (1.0 - betaM)
                } ) 

          val newBMState = 
                  bmState.copy ( bmInteractions = BMInteractions(phia, phijab),
                                grad = grad,
                                v = BMInteractions(g2phia, g2phijab), m = BMInteractions(dphia2, dphijab2),
                                betaV = betaV, betaM = betaM,
                                learningRate = bmState.learningRate, 
                                learningRates = learningRates,
                                eps = eps, step = bmState.step + 1 )

        //val optionLearningRates: Option[BMInteractions] = None
        //(newBMState, optionLearningRates)     // 240428: ??? Adadelta canot be used for the soft thresholding function for L1.
          (newBMState, Option(learningRates) )     // 240428: ??? Adadelta canot be used for the soft thresholding function for L1.
    }

    def learningByAdam(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1, outDetails: Boolean = true ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {

          import breeze.linalg.InjectNumericOps

          assert(minOrMax.abs == 1)

        //val step = bmState.step + 1
          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 ) 

          val learningRate = bmState.learningRate.lr(step)
          val betaM = bmState.betaM
          val betaV = bmState.betaV
          val eps = bmState.eps
          assert(0.0 <= betaM && betaM < 1.0)
          assert(0.0 <= betaV && betaV < 1.0)

          val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val betaMt = 1.0 - math.pow(betaM, step)
          val betaVt = 1.0 - math.pow(betaV, step)

          val invSqrtVtPhia = vphia.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
          val invSqrtVtPhijab = vphijab.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 

          val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val invVt = invSqrtVtPhia(i)
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 
          val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val invVt = invSqrtVtPhijab(ij)
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 

          val learningRatePhia = invSqrtVtPhia.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )
          val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )

          if ( outDetails ) {
            val vtphia = vphia.map( v => breeze.linalg.max(v) ).max
            val vtphijab = vphijab.map( v =>  breeze.linalg.max(v) ).max
            val vtphi = math.max(vtphia, vtphijab) / betaVt

            val invSqrtVtPhi = 1.0 / (math.sqrt(vtphi) + eps)

            val learningR = invSqrtVtPhi * learningRate

            val vtphia_sum = vphia.map( v => breeze.linalg.sum(v) ).sum
            val vtphijab_sum = vphijab.map( v =>  breeze.linalg.sum(v) ).sum
            val vtphi_sum = (vtphia_sum + vtphijab_sum) / betaVt

            val sqrt_vtphi_sum = math.sqrt(vtphi_sum)
            val sqrt_vtphi = math.sqrt(vtphi)

            stderr.print("Adam: step= %d  %d  sqrt(max_i v^_i)=   %f  sqrt(sum_i v^_i)=   %f  ratio= %f  learning rate for ModAdam= %f  learning rate= %f\n".format(
                          step, bmState.step + 1, sqrt_vtphi, sqrt_vtphi_sum, sqrt_vtphi / sqrt_vtphi_sum, learningR, learningRate ) )
        //}
        //
        //if ( outDetails ) {
            import breeze.linalg.{max, sum}
            import breeze.numerics.abs

          //val mtphia = mphia.map( m => breeze.linalg.max(m *:* m) ).max
          //val mtphijab = mphijab.map( m =>  breeze.linalg.max(m *:* m) ).max
          //val mtphi = math.max(mtphia, mtphijab) / (betaMt * betaMt)

          //val sqrt_mtphi = math.sqrt(mtphi)

            val mtphia = mphia.map( m => max(abs(m)) ).max
            val mtphijab = mphijab.map( m =>  max(abs(m)) ).max
            val mtphi = math.max(mtphia, mtphijab) / betaMt

            val sqrt_mtphi = mtphi

            val mtphia_sum = mphia.map( m => breeze.linalg.sum(m *:* m) ).sum
            val mtphijab_sum = mphijab.map( m =>  breeze.linalg.sum(m *:* m) ).sum
            val mtphi_sum = (mtphia_sum + mtphijab_sum) / (betaMt * betaMt)

            val sqrt_mtphi_sum = math.sqrt(mtphi_sum)

            stderr.print("Adam: step= %d  %d  sqrt(max_i m^_i^2)= %f  sqrt(sum_i m^_i^2)= %f  ratio= %f  learning rate for ModAdam= %f  learning rate= %f\n".format(
                          step, bmState.step + 1, sqrt_mtphi, sqrt_mtphi_sum, sqrt_mtphi / sqrt_mtphi_sum, learningR, learningRate ) )
          }
        /**/

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
        //val optionLearningRates: Option[BMInteractions] = None

          ( bmState.copy (
                bmInteractions = BMInteractions(phia, phijab),
                grad = grad,
                v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
                betaV = betaV, betaM = betaM,
                learningRate = bmState.learningRate,
                learningRates = learningRates,
                eps = eps, step = bmState.step + 1),
              //optionLearningRates     )       //240428: ??? Adam cannot be used for the soft thresholding function for L1.
                Option(learningRates)     )
    }

    def learningByRAdam(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {
        //RAdam         arXiv:1908.03265

          import breeze.linalg.InjectNumericOps

          assert(minOrMax.abs == 1)
          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )
          val learningRate = bmState.learningRate.lr(step)
          val betaM = bmState.betaM
          val betaV = bmState.betaV
          val eps = bmState.eps
          assert(0.0 <= betaM && betaM < 1.0)
          assert(0.0 <= betaV && betaV < 1.0)

          val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val betaMt = 1.0 - math.pow(betaM, step)
          val betaV_power_t = math.pow(betaV, step)
          val betaVt = 1.0 - betaV_power_t

          val rhoInfty = 2.0 / (1.0 - betaV) - 1.0
          val rhot = rhoInfty - 2.0 * (step * betaV_power_t) / betaVt           // rectification term to Adam

          val ( phia, phijab, learningRatePhia, learningRatePhijab ) = if ( rhot > 4.0 ) {

            val rt = math.sqrt( ( (rhot - 4.0 ) * (rhot - 2.0) * rhoInfty ) / ( (rhoInfty - 4.0) * (rhoInfty - 2.0) * rhot ) )

            val invSqrtVtPhia = vphia.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
            val invSqrtVtPhijab = vphijab.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 

            val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val invVt = invSqrtVtPhia(i)
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * rt * ( m *:* invVt )
                        s + ( m *:* invVt ) * ( (minOrMax * learningRate / betaMt) * rt )
                } ) 
            val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val invVt = invSqrtVtPhijab(ij)
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * rt * ( m *:* invVt )
                        s + ( m *:* invVt ) * ( (minOrMax * learningRate / betaMt) * rt )
                } ) 

            val learningRatePhia = invSqrtVtPhia.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * rt * invVt
                                invVt * ( (learningRate * (1.0 - betaM) / betaMt) * rt )
                        } )
            val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * rt * invVt
                                invVt * ( (learningRate * (1.0 - betaM) / betaMt) * rt )
                        } )

            ( phia, phijab, learningRatePhia, learningRatePhijab ) 

          } else {

            val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * ( m )
                        s + ( m ) * (minOrMax * learningRate / betaMt)
                } ) 
            val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * ( m )
                        s + ( m ) * (minOrMax * learningRate / betaMt)
                } ) 

            val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt)

            val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
            val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

            ( phia, phijab, learningRatePhia, learningRatePhijab ) 

          }

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
        //val optionLearningRates: Option[BMInteractions] = None

          ( bmState.copy (
                bmInteractions = BMInteractions(phia, phijab),
                grad = grad,
                v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
                betaV = betaV, betaM = betaM,
                learningRate = bmState.learningRate,
                learningRates = learningRates,
                eps = eps, step = bmState.step + 1),
              //optionLearningRates     )       //240428: ??? RAdam cannot be used for the soft thresholding function for L1.
                Option(learningRates)     )       //240428: ??? RAdam cannot be used for the soft thresholding function for L1.
    }

    def learningByModAdamMax(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1, outDetails: Boolean = true ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {

          import breeze.linalg.{InjectNumericOps, max, sum}

          assert(minOrMax.abs == 1)
          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )
          val learningRate = bmState.learningRate.lr(step)
          val betaM = bmState.betaM
          val betaV = bmState.betaV
          val eps = bmState.eps
          assert(0.0 <= betaM && betaM < 1.0)
          assert(0.0 <= betaV && betaV < 1.0)

          val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val betaMt = 1.0 - math.pow(betaM, step)
          val betaVt = 1.0 - math.pow(betaV, step)

        /*
          val invSqrtVtPhia = vphia.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
          val invSqrtVtPhijab = vphijab.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
          val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val invVt = invSqrtVtPhia(i)
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 
          val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val invVt = invSqrtVtPhijab(ij)
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 

          val learningRatePhia = invSqrtVtPhia.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )
          val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )
        */
        //
          val vtphia = vphia.map( v => breeze.linalg.max(v) ).max 
          val vtphijab = vphijab.map( v =>  breeze.linalg.max(v) ).max
          val vtphi = math.max(vtphia, vtphijab) / betaVt

          val invSqrtVtPhi = 1.0 / (math.sqrt(vtphi) + eps)

          val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m 
                        s + m * (invSqrtVtPhi * (minOrMax * learningRate / betaMt) )
                } ) 
          val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m
                        s + m * (invSqrtVtPhi * (minOrMax * learningRate / betaMt) )
                } ) 

          //val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt) * invSqrtVtPhi
          val actualLearnRate = invSqrtVtPhi * ( learningRate * ((1.0 - betaM) / betaMt) )

          val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
          val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

        /**/
          //val learningR = learningRate * invSqrtVtPhi
          val learningR = invSqrtVtPhi * learningRate

          if ( outDetails ) {
            val vtphia_sum = vphia.map( v => breeze.linalg.sum(v) ).sum 
            val vtphijab_sum = vphijab.map( v =>  breeze.linalg.sum(v) ).sum
            val vtphi_sum = (vtphia_sum + vtphijab_sum) / betaVt

            val sqrt_vtphi_sum = math.sqrt(vtphi_sum)
            val sqrt_vtphi = math.sqrt(vtphi)
            stderr.print("ModAdam: step= %d  %d  sqrt(max_i v^_i)=   %f  sqrt(sum_i v^_i)=   %f  ratio= %f  learning rate for ModAdam= %f  learning rate= %f\n".format(
                          step, bmState.step + 1, sqrt_vtphi, sqrt_vtphi_sum, sqrt_vtphi / sqrt_vtphi_sum, learningR, learningRate ) )
          }

          if ( outDetails ) {
            import breeze.linalg.{max, sum}
            import breeze.numerics.abs

          //val mtphia = mphia.map( m => breeze.linalg.max(m *:* m) ).max 
          //val mtphijab = mphijab.map( m =>  breeze.linalg.max(m *:* m) ).max
          //val mtphi = math.max(mtphia, mtphijab) / (betaMt * betaMt)

          //val sqrt_mtphi = math.sqrt(mtphi)

            val mtphia = mphia.map( m => max(abs(m)) ).max 
            val mtphijab = mphijab.map( m =>  max(abs(m)) ).max
            val mtphi = math.max(mtphia, mtphijab) / betaMt

            val sqrt_mtphi = mtphi

            val mtphia_sum = mphia.map( m => breeze.linalg.sum(m *:* m) ).sum 
            val mtphijab_sum = mphijab.map( m =>  breeze.linalg.sum(m *:* m) ).sum
            val mtphi_sum = (mtphia_sum + mtphijab_sum) / (betaMt * betaMt)

            val sqrt_mtphi_sum = math.sqrt(mtphi_sum)

            stderr.print("ModAdam: step= %d  %d  sqrt(max_i m^_i^2)= %f  sqrt(sum_i m^_i^2)= %f  ratio= %f  learning rate for ModAdam= %f  learning rate= %f\n".format(
                          step, bmState.step + 1, sqrt_mtphi, sqrt_mtphi_sum, sqrt_mtphi / sqrt_mtphi_sum, learningR, learningRate ) )
          }
        /**/

        //

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)

          ( bmState.copy (
                        bmInteractions = BMInteractions(phia, phijab),
                        grad = grad,
                        v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
                        betaV = betaV, betaM = betaM,
                        learningRate = bmState.learningRate,
                        learningRates = learningRates,
                        eps = eps, step = bmState.step + 1),
                        Option(learningRates) )

    }

    def learningByModAdamSum(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1, outDetails: Boolean = true ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {

          import breeze.linalg.{InjectNumericOps, max, sum}

          assert(minOrMax.abs == 1)
          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )
          val learningRate = bmState.learningRate.lr( step )
          val betaM = bmState.betaM
          val betaV = bmState.betaV
          val eps = bmState.eps
          assert(0.0 <= betaM && betaM < 1.0)
          assert(0.0 <= betaV && betaV < 1.0)

          val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val betaMt = 1.0 - math.pow(betaM, step)
          val betaVt = 1.0 - math.pow(betaV, step)

        /*
          val invSqrtVtPhia = vphia.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
          val invSqrtVtPhijab = vphijab.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
          val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val invVt = invSqrtVtPhia(i)
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 
          val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val invVt = invSqrtVtPhijab(ij)
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
                        s + ( m *:* invVt ) * (minOrMax * learningRate / betaMt)
                } ) 

          val learningRatePhia = invSqrtVtPhia.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )
          val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * invVt
                                invVt * (learningRate * (1.0 - betaM) / betaMt)
                        } )
        */
        //
        /* for ModAdamMax
          val vtphia = vphia.map( v => breeze.linalg.max(v) ).max 
          val vtphijab = vphijab.map( v =>  breeze.linalg.max(v) ).max
          val vtphi = math.max(vtphia, vtphijab) / betaVt
        */
          val vtphia = vphia.map( v => breeze.linalg.sum(v) ).sum 
          val vtphijab = vphijab.map( v =>  breeze.linalg.sum(v) ).sum
          val vtphi = (vtphia + vtphijab) / betaVt

          val invSqrtVtPhi = 1.0 / (math.sqrt(vtphi) + eps)

          val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m 
                        s + m * ( invSqrtVtPhi * (minOrMax * learningRate / betaMt) )
                } ) 
          val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m
                        s + m * ( invSqrtVtPhi * (minOrMax * learningRate / betaMt) )
                } ) 

          //val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt) * invSqrtVtPhi
          val actualLearnRate = invSqrtVtPhi * ( learningRate * ((1.0 - betaM) / betaMt) )

          val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
          val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

        /**/
          //val learningR = learningRate * invSqrtVtPhi
          val learningR = invSqrtVtPhi * learningRate

          if ( outDetails ) {
            val vtphia_max = vphia.map( v => breeze.linalg.max(v) ).max 
            val vtphijab_max = vphijab.map( v =>  breeze.linalg.max(v) ).max
            val vtphi_max = math.max(vtphia_max, vtphijab_max) / betaVt

            val sqrt_vtphi_max = math.sqrt(vtphi_max)
            val sqrt_vtphi = math.sqrt(vtphi)
            stderr.print("step = %d  sqrt(max_i v^_i)  = %f  sqrt(sum_i v^_i)  = %f  ratio= %f  learnng rate= %f\n".format( step, sqrt_vtphi_max, sqrt_vtphi, sqrt_vtphi_max / sqrt_vtphi, learningR ) )
          }

          if ( outDetails ) {
            import breeze.linalg.{max, sum}
            import breeze.numerics.abs

          //val mtphia = mphia.map( m => breeze.linalg.max(m *:* m) ).max 
          //val mtphijab = mphijab.map( m =>  breeze.linalg.max(m *:* m) ).max
          //val mtphi = math.max(mtphia, mtphijab) / (betaMt * betaMt)

          //val sqrt_mtphi = math.sqrt(mtphi)

            val mtphia = mphia.map( m => max(abs(m)) ).max 
            val mtphijab = mphijab.map( m =>  max(abs(m)) ).max
            val mtphi = math.max(mtphia, mtphijab) / betaMt

            val sqrt_mtphi = mtphi

            val mtphia_sum = mphia.map( m => breeze.linalg.sum(m *:* m) ).sum 
            val mtphijab_sum = mphijab.map( m =>  breeze.linalg.sum(m *:* m) ).sum
            val mtphi_sum = (mtphia_sum + mtphijab_sum) / (betaMt * betaMt)

            val sqrt_mtphi_sum = math.sqrt(mtphi_sum)
            stderr.print("step= %d  sqrt(max_i m^_i^2)= %f  sqrt(sum_i m^_i^2)= %f  ratio= %f  learning rate= %f\n".format( step, sqrt_mtphi, sqrt_mtphi_sum, sqrt_mtphi / sqrt_mtphi_sum, learningR ) )
          }
        /**/

        //

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)

          ( bmState.copy (
                        bmInteractions = BMInteractions(phia, phijab),
                        grad = grad,
                        v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
                        betaV = betaV, betaM = betaM,
                        learningRate = bmState.learningRate,
                        learningRates = learningRates,
                        eps = eps, step = bmState.step + 1),
                        Option(learningRates) )

    }

    def learningByModRAdamMax(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1, outDetails: Boolean = true ):
                                                        Tuple2[BMState, Option[BMInteractions]] = {
        //RAdam         arXiv:1908.03265

          import breeze.linalg.{InjectNumericOps, max, sum}

          assert(minOrMax.abs == 1)
          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )
          val learningRate = bmState.learningRate.lr( step )
          val betaM = bmState.betaM
          val betaV = bmState.betaV
          val eps = bmState.eps
          assert(0.0 <= betaM && betaM < 1.0)
          assert(0.0 <= betaV && betaV < 1.0)

          val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
                        val m = dx._1
                        val g = dx._2
                        //betaM * m + (1.0 - betaM) * g
                        m * betaM + g * (1.0 - betaM)
                } ) 
          val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
                        val v = dx._1
                        val g = dx._2
                        //betaV * v + (1.0 - betaV) * (g *:* g)
                        v * betaV + (g *:* g) * (1.0 - betaV)
                } ) 
          val betaMt = 1.0 - math.pow(betaM, step)
          val betaV_power_t = math.pow(betaV, step)
          val betaVt = 1.0 - betaV_power_t

          val rhoInfty = 2.0 / (1.0 - betaV) - 1.0
          val rhot = rhoInfty - 2.0 * (step * betaV_power_t) / betaVt           // rectification term to Adam

          val ( phia, phijab, learningRatePhia, learningRatePhijab, learningR ) = if ( rhot > 4.0 ) {

            val rt = math.sqrt( ( (rhot - 4.0 ) * (rhot - 2.0) * rhoInfty ) / ( (rhoInfty - 4.0) * (rhoInfty - 2.0) * rhot ) )

         /*
            val invSqrtVtPhia = vphia.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 
            val invSqrtVtPhijab = vphijab.map( v => {
                                val vt = v / betaVt
                                1.0 / (breeze.numerics.sqrt(vt) + eps)
                        } ) 

            val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val invVt = invSqrtVtPhia(i)
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * rt * ( m *:* invVt )
                        s + ( m *:* invVt ) * ( (minOrMax * learningRate / betaMt) * rt )
                } ) 
            val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val invVt = invSqrtVtPhijab(ij)
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * rt * ( m *:* invVt )
                        s + ( m *:* invVt ) * ( (minOrMax * learningRate / betaMt) * rt )
                } ) 

            val learningRatePhia = invSqrtVtPhia.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * rt * invVt
                                invVt * ( (learningRate * (1.0 - betaM) / betaMt) * rt )
                        } )
            val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
                                //(learningRate * (1.0 - betaM) / betaMt) * rt * invVt
                                invVt * ( (learningRate * (1.0 - betaM) / betaMt) * rt )
                        } )
         */
        //
            val vtphia = vphia.map( v => breeze.linalg.max(v) ).max 
            val vtphijab = vphijab.map( v =>  breeze.linalg.max(v) ).max
            val vtphi = math.max(vtphia, vtphijab) / betaVt

            val invSqrtVtPhi = 1.0 / (math.sqrt(vtphi) + eps)

            val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi * rt ) * m 
                        s + m * ( invSqrtVtPhi * (minOrMax * learningRate / betaMt * rt ) )
                } ) 
            val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt * invSqrtVtPhi * rt ) * m
                        s + m * ( invSqrtVtPhi * (minOrMax * learningRate / betaMt * rt ) )
                } ) 

            //val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt) * invSqrtVtPhi * rt
            val actualLearnRate = invSqrtVtPhi * ( learningRate * ((1.0 - betaM) / betaMt) * rt )

            val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
            val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

            //val learningR = learningRate * invSqrtVtPhi * rt
            val learningR = invSqrtVtPhi * ( learningRate * rt )

        /**/
            if ( outDetails) {
              val vtphia_sum = vphia.map( v => breeze.linalg.sum(v) ).sum 
              val vtphijab_sum = vphijab.map( v =>  breeze.linalg.sum(v) ).sum
              val vtphi_sum = (vtphia_sum + vtphijab_sum) / betaVt

              val sqrt_vtphi_sum = math.sqrt(vtphi_sum)
              val sqrt_vtphi = math.sqrt(vtphi)
              stderr.print("step = %d  sqrt(max_i v^_i)  = %f  sqrt(sum_i v^_i)  = %f  ratio= %f  learning rate= %f\n".format( 
                        step, sqrt_vtphi, sqrt_vtphi_sum, sqrt_vtphi / sqrt_vtphi_sum, learningR ) )
            }
        /**/

        //

            ( phia, phijab, learningRatePhia, learningRatePhijab, learningR ) 

          } else {

            val phia = mphia.zipWithIndex.map( x => {
                        val (m, i) = x
                        val s = bmState.bmInteractions.phia(i)
                        //s + (minOrMax * learningRate / betaMt) * ( m )
                        s + ( m ) * (minOrMax * learningRate / betaMt)
                } ) 
            val phijab = mphijab.zipWithIndex.map( x => {
                        val (m, ij) = x
                        val s = bmState.bmInteractions.phijab(ij)
                        //s + (minOrMax * learningRate / betaMt) * ( m )
                        s + ( m ) * (minOrMax * learningRate / betaMt)
                } ) 

            val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt)

            val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
            val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

            ( phia, phijab, learningRatePhia, learningRatePhijab, learningRate ) 

          }

        /**/
          if ( outDetails ) {
            import breeze.linalg.{max, sum}
            import breeze.numerics.abs

          //val mtphia = mphia.map( m => breeze.linalg.max(m *:* m) ).max 
          //val mtphijab = mphijab.map( m =>  breeze.linalg.max(m *:* m) ).max
          //val mtphi = math.max(mtphia, mtphijab) / (betaMt * betaMt)

          //val sqrt_mtphi = math.sqrt(mtphi)

            val mtphia = mphia.map( m => max(abs(m)) ).max 
            val mtphijab = mphijab.map( m =>  max(abs(m)) ).max
            val mtphi = math.max(mtphia, mtphijab) / betaMt

            val sqrt_mtphi = mtphi

            val mtphia_sum = mphia.map( m => breeze.linalg.sum(m *:* m) ).sum 
            val mtphijab_sum = mphijab.map( m =>  breeze.linalg.sum(m *:* m) ).sum
            val mtphi_sum = (mtphia_sum + mtphijab_sum) / (betaMt * betaMt)

            val sqrt_mtphi_sum = math.sqrt(mtphi_sum)
            stderr.print("step= %d  sqrt(max_i m^_i^2)= %f  sqrt(sum_i m^_i^2)= %f  ratio= %f  learning rate= %f\n".format( step, sqrt_mtphi, sqrt_mtphi_sum, sqrt_mtphi / sqrt_mtphi_sum, learningR ) )
          }
        /**/

          val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
          val optionLearningRates: Option[BMInteractions] = None

          ( bmState.copy (
                bmInteractions = BMInteractions(phia, phijab),
                grad = grad,
                v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
                betaV = betaV, betaM = betaM,
                learningRate = bmState.learningRate,
                learningRates = learningRates,
                eps = eps, step = bmState.step + 1),
              //optionLearningRates     )       //240428: ??? RAdam cannot be used for the soft thresholding function for L1.
                Option(learningRates)     )       //240428: ??? RAdam cannot be used for the soft thresholding function for L1.
    }

    // J. P. Barton1, E. De Leonardis, A. Coucke, and S. Cocco
    // 16_Bioinformatics_32_3089-3097 and supplementary_data
    // By Matteo Figliuzzi, Pierre Barrat-Charlaix, and Martin Weigt
    // Mol. Biol. Evol. 35, 1018-1027, 2018
    def learningByRPROPLR(bmState: BMState, grad: BMInteractions, 
                        regTerm: String, propL1h: Double, propL1J: Double, minOrMax: Int = 1 ):
                                                Tuple2[BMState, Option[BMInteractions]] = {

          import breeze.linalg.InjectNumericOps

        /*
          val minLearnRatePhia = 0.001
          val minLearnRatePhijab = 0.00001
          val maxLearnRatePhia = 2.50
          val maxLearnRatePhijab = 2.50
        */
        /*
          val minLearnRatePhia = bmState.learningRate * 0.001
          val minLearnRatePhijab = bmState.learningRate * 0.001
          val maxLearnRatePhia = bmState.learningRate * 10.0
          val maxLearnRatePhijab = bmState.learningRate * 10.0
        */
        /*
          val minLearnRatePhia = bmState.learningRate * 0.01
          val minLearnRatePhijab = bmState.learningRate * 0.01
          val maxLearnRatePhia = bmState.learningRate * 100.0
          val maxLearnRatePhijab = bmState.learningRate * 100.0
        */
        /*
          val minLearnRatePhia = bmState.learningRate * 0.001
          val minLearnRatePhijab = bmState.learningRate * 0.001
          val maxLearnRatePhia = bmState.learningRate * 1000.0
          val maxLearnRatePhijab = bmState.learningRate * 1000.0
        */
        /*
          val minLearnRatePhia = bmState.learningRate * 0.0001
          val minLearnRatePhijab = bmState.learningRate * 0.0001
          val maxLearnRatePhia = bmState.learningRate * 1000.0
          val maxLearnRatePhijab = bmState.learningRate * 1000.0
        */
          val minLearnRatePhia = bmState.learningRateForRPROPLR.minLearningRate
          val minLearnRatePhijab = bmState.learningRateForRPROPLR.minLearningRate
          val maxLearnRatePhia = bmState.learningRateForRPROPLR.maxLearningRate
          val maxLearnRatePhijab = bmState.learningRateForRPROPLR.maxLearningRate
        /*
          val rateRange = (r: Double, minR: Double, maxR: Double) =>  
                                        if ( r > maxR )
                                                maxR
                                        else if ( r < minR )
                                                minR 
                                        else
                                                r
        */
          val rateRangePhia = bmState.learningRateForRPROPLR.lr(_: Double, minLearnRatePhia, maxLearnRatePhia )
          val rateRangePhijab = bmState.learningRateForRPROPLR.lr(_: Double, minLearnRatePhijab, maxLearnRatePhijab )

          assert(minOrMax.abs == 1)
          val rateDecrease = bmState.learningRateForRPROPLR.rateDecrease       // bmState.betaV
          val rateIncrease = bmState.learningRateForRPROPLR.rateIncrease       // bmState.betaM
          //assert( rateIncrease > rateDecrease )
          assert( rateIncrease > 1.0 && rateDecrease < 1.0 )
          val prevLearningRates = bmState.learningRates // bmState.m
          val prevGrad = bmState.grad

          val step = bmState.step + 1 - bmState.stepOffset
          assert( step >= 1 )

          val learningRatePhia = prevLearningRates.phia.zipWithIndex.map( x => {
                                val rate = x._1
                                val i = x._2
                                val change = (grad.phia(i) *:* prevGrad.phia(i)).map( 
                                        x => {if ( x > 0 ) rateIncrease else if ( x < 0 ) rateDecrease else 1.0} ) 
                                (rate *:* change).map( r => rateRangePhia(r) )
                        } )

          val learningRatePhijab = prevLearningRates.phijab.zipWithIndex.map( x => {
                                val rate = x._1
                                val ij = x._2
                                val change = (grad.phijab(ij) *:* prevGrad.phijab(ij)).map(
                                        x => {if ( x > 0 ) rateIncrease else if ( x < 0 ) rateDecrease else 1.0} ) 
                                (rate *:* change).map( r => rateRangePhijab(r) )
                        } )

          val phia = bmState.bmInteractions.phia.zipWithIndex.map( x => {
                                val phi = x._1
                                val i = x._2
                                phi + (minOrMax.toDouble * learningRatePhia(i)) *:* grad.phia(i)
                        } )

          val phijab = bmState.bmInteractions.phijab.zipWithIndex.map( x => {
                                val phij = x._1
                                val ij = x._2
                                phij + (minOrMax.toDouble * learningRatePhijab(ij)) *:* grad.phijab(ij)
                        } )

          val learningRates = 
                if ( regTerm == "GL1" || regTerm == "GL1L2" ) {
                        val modLearningRatePhia =
                                learningRatePhia
                                /*      does not converge.
                                if ( propL1h <= 0.0 ) {
                                        learningRatePhia
                                } else {
                                        learningRatePhia.map( phi => {
                                                //val minRate = breeze.linalg.min(phi)
                                                //phi.map(minRate)
                                                val aveRate = breeze.linalg.sum(phi) / phi.size
                                                phi.map(x => aveRate)
                                        } )
                                }
                                */

                        val modLearningRatePhijab =
                                learningRatePhijab
                                /*      does not converge.
                                if (propL1J <= 0.0) {
                                        learningRatePhijab
                                } else {
                                        learningRatePhijab.map ( phij => {
                                                //val minRate = breeze.linalg.min(phij)
                                                //phij.map(minRate)
                                                val aveRate = breeze.linalg.sum(phij) / phij.size
                                                phij.map(x => aveRate)
                                        } )
                                }
                                */

                        BMInteractions(modLearningRatePhia, modLearningRatePhijab)
                } else {
                        BMInteractions(learningRatePhia, learningRatePhijab)
                }

          val newBMState =
                  bmState.copy (
                        bmInteractions = BMInteractions(phia, phijab),
                        grad = grad,
                        learningRateForRPROPLR = bmState.learningRateForRPROPLR,
                      //minLearningRate = bmState.minLearningRate, maxLearningRate = bmState.maxLearningRate,
                      //rateDecrease = bmState.rateDecrease, rateIncrease = bmState.rateIncrease,
                        learningRates = learningRates,
                        eps = bmState.eps, step = bmState.step + 1 )

          (newBMState, Option(learningRates))   //240428: ??? RPROP-LR cannot be used for the soft thresholding function for L1.
    }

    def initializeBMStateForNAG(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.9, betaM: Double = 0.0, eps: Double = 1.0e-6, step: Int = 0 ) = {

        initializeBMStateForAdam(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForAdam(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {

        val nUnits = bmInteractions.phia.size 
        val nPairs = (nUnits * (nUnits - 1)) / 2
        val nStatesofUnit = bmInteractions.phia(0).size

        val grad = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
                 ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )
        val v = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
                 ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) ) 
        val m = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
                 ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )

        val step1LR = learningRate.lr(step = 1)
        val learningRates = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.fill(nStatesofUnit)(step1LR) ), 
            ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.fill(nStatesofUnit, nStatesofUnit)(step1LR) ) )

        val bmState = BM.createBMState(gradientDescentMethod = gradientDescentMethod, bmInteractions = bmInteractions, 
                grad = grad, v = v, m = m, betaV = betaV, betaM = betaM, eps = eps,
                learningRate = learningRate, learningRates = learningRates, step = step )

        bmState
    }

    def initializeBMStateForRAdam(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {

        initializeBMStateForAdam(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForModAdamMax(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {

        initializeBMStateForAdam(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForModAdamSum(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {

        initializeBMStateForModAdamMax(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForModRAdamMax(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {

        initializeBMStateForRAdam(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForAdadelta(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRate,
        betaV: Double = 0.9, betaM: Double = 0.9, eps: Double = 1.0e-5, step: Int = 0 ) = {

        assert( betaM == betaV )
        initializeBMStateForAdam(gradientDescentMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForRPROPLR(gradientDescentMethod: String, bmInteractions: BMInteractions, learningRate: LearningRateForRPROPLR,
      //minLearningRate: Double = 0.000001, maxLearningRate: Double = 10.0,
      //rateDecrease: Double = 0.5, rateIncrease: Double = 1.2,
        eps: Double = 1.0e-6, step: Int = 0 ) = {
        // betaV: rateDecrease , betaM: rateInccrease
        // bmState.m: prevLearningRates

        val nUnits = bmInteractions.phia.size 
        val nPairs = (nUnits * (nUnits - 1)) / 2
        val nStatesofUnit = bmInteractions.phia(0).size

        val minLearningRate = learningRate.minLearningRate

        val grad = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
            ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )

        val learningRates = BM.BMInteractions( ArraySeq.range(0, nUnits).map ( ve => DenseVector.fill(nStatesofUnit)(minLearningRate) ), 
            ArraySeq.range(0, nPairs).map ( ma => DenseMatrix.fill(nStatesofUnit, nStatesofUnit)(minLearningRate) ) )

        val bmState = BM.createBMState( gradientDescentMethod = gradientDescentMethod, bmInteractions = bmInteractions, 
                grad = grad, 
                eps = eps,
                learningRateForRPROPLR = learningRate,
              //minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
              //rateDecrease = rateDecrease, rateIncrease = rateIncrease,
                learningRates = learningRates, step = step )

        bmState
    }

    /* sigma_J: standard deviation for the elements of J; with mean = 0. */
    def initializehJ ( fia: ArraySeq[DenseVector[Double] ], observedN: Double, pseudoN: Double = 100.0, sigma_J: Double = 0.01) = {

      //import breeze.stats.distributions.Rand.FixedSeed.*            //for scala 3
        import SetSeedForRand.randBasis
        import breeze.stats.distributions.Gaussian
        import breeze.linalg.{InjectNumericOps, sum}
        import breeze.numerics.log

        require( pseudoN >= 0.0 )
        val nStatesOfUnit = fia(0).size
        val nUnits = fia.size 
        val nPairs = (nUnits * (nUnits - 1)) / 2

        val Jijab = ArraySeq.range(0, nPairs).map( _ => DenseMatrix.zeros[Double](nStatesOfUnit, nStatesOfUnit) )

      /*
        val ratio = observedN / (observedN + pseudoN)
        val pseudo = (1 - ratio) / nStatesOfUnit

        val fiaWithPseudo = if ( pseudo <= 0.0 ) {
                              fia
                            } else {
                              fia.map ( fi => fi * ratio + pseudo  )
                            }
      */

        val fiaWithPseudo = bayesianCorrection(fia, observedN, pseudoN ) 

        val hia = fiaWithPseudo.map ( fi => { 
                        val hi = breeze.numerics.log(fi)
                        hi - breeze.linalg.sum(hi) / hi.size    // Ising gauge
                } )

        val jijab = if ( sigma_J <= 0.0 ) {
                      Jijab
                    } else {
                      val g = Gaussian(0.0, sigma_J) 
                      Jijab.map{ mat => mat.map{ ab => g.draw() } }
                    } 
        (hia, jijab)
      //val (hia_Ising, jijab_Ising) = toIsingGauge( hia, Jijab )
      //(hia_Ising, jijab_Ising)
    }

    def initializeBMState(
        optMethod: String,
        bmInteractions: BMInteractions,
        learningRate: LearningRate,
        learningRateForRPROPLR: LearningRateForRPROPLR,
      //minLearningRate: Double, maxLearningRate: Double,
      //rateDecrease: Double, rateIncrease: Double,
        betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6, step: Int = 0 ) = {   // for Adam

        assert( step >= 0 )
        val bmState = {
              if ( optMethod == "NAG" )
                initializeBMStateForNAG(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "ModRAdamMax" || optMethod == "ModRAdam" )
                initializeBMStateForModRAdamMax(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "ModAdamMax" || optMethod == "ModAdam" )
                initializeBMStateForModAdamMax(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "ModAdamSum" )
                initializeBMStateForModAdamSum(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "Adam" )
                initializeBMStateForAdam(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "RAdam" )
                initializeBMStateForRAdam(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if ( optMethod == "Adadelta" )
                initializeBMStateForAdadelta(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              else if  ( optMethod == "RPROP-LR" || optMethod == "MF" )
                      initializeBMStateForRPROPLR(optMethod, bmInteractions = bmInteractions, learningRate = learningRateForRPROPLR, 
                        //minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
                        //rateDecrease = rateDecrease, rateIncrease = rateIncrease, 
                          eps = eps, step = step )
              else {
                sys.error("Not supported: " + optMethod) 
                initializeBMStateForAdam(optMethod, bmInteractions, learningRate, betaV, betaM, eps, step )
              }
        }

        bmState.copy( learningRate = learningRate, learningRateForRPROPLR = learningRateForRPROPLR )
    }

    def initialize(fia: ArraySeq[DenseVector[Double]],
        initialInteractions: MCMC.Interactions,
        optMethod: String,
        learningRate: LearningRate,
        learningRateForRPROPLR: LearningRateForRPROPLR,
      //minLearningRate: Double, maxLearningRate: Double,
      //rateDecrease: Double, rateIncrease: Double,
        betaV: Double, betaM: Double, eps: Double = 1.0e-6, step: Int = 0,
        betaLAP: Double = 0.0, betaLAPForKL: Double = 0.0 ) = {

      //val nUnits = fia.size 
      //val nPairs = (nUnits * (nUnits - 1)) / 2
      //val nStatesofUnit = fia(0).size
        
        val interactions = initialInteractions

        val (phia, phijab) = BM.hJToPhi (interactions.hia, interactions.jijab, fia )
        val bmInteractions = BM.BMInteractions( phia, phijab, interactions.gauge )

        val bmState = initializeBMState(
                          optMethod,
                          bmInteractions,
                          learningRate,
                          learningRateForRPROPLR,
                        //minLearningRate: Double, maxLearningRate: Double,
                        //rateDecrease: Double, rateIncrease: Double,
                          betaV, betaM, eps, step )

        (interactions, bmState.copy(betaLAP = betaLAP, betaLAPForKL = betaLAPForKL) )
    }

    def runBM1Step( ioFiles: HashMap[String, File],
                stateOrderString: String,
                interactions: MCMC.Interactions,
              //initialConfigurations: ArraySeq[IArray[Byte]],
              //initialMCStates: ArraySeq[MCMC.State],
                nInitialIterationsPerUnit: Int = 0, 
                everyNIterationsPerUnit: Int = 10,
                nSamples: Int = 1,
                initialT: Double = 1.0, 
                finalT: Double = 1.0, 
                annealingRate: Double = 0.99,
                maxExtendedIterations: Int = 1,
              //kernel: (MCMC.State, RandBasis) => Rand[MCMC.State],
                mcmcKernel: String = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"
                proposedPia: ArraySeq[DenseVector[Double]],       // fia corrected with pseudocounts

                fia: ArraySeq[DenseVector[Double]],
                fijab: ArraySeq[DenseMatrix[Double]],

                regTerm: String,
                propL1h: Double, propL1J: Double,

                lambdaPhi: Double, lambdaPhij: Double,
                bmState: BMState,
                optMethod: String,

                noSoftThresholding: Boolean = false
                ) = {

        require( mcmcKernel == "MH" || 
                        mcmcKernel == "GibbsWithMHStep" || mcmcKernel == "MultiBlockMH" || 
                        mcmcKernel == "Gibbs", 
                        "Not supported: mcmcKernel == MH, GibbsWithMHStep or Gibbs" )

        require( lambdaPhi >= 0.0 && lambdaPhij >= 0.0 )

        val mcmc = new MCMC(ioFiles, stateOrderString, interactions, proposedPia = proposedPia )
        val kernelForMCMC: (MCMC.State, RandBasis) => Rand[MCMC.State] =
          if ( mcmcKernel == "MH" ) {
              //mcmc.kernelMH(_)(_)
                mcmc.kernelMH(_)(using _: RandBasis)
          } else if (  mcmcKernel == "GibbsWithMHStep" || mcmcKernel == "MultiBlockMH" ) {
              //mcmc.kernelGibbsWithMHStep(_)(_)
                mcmc.kernelGibbsWithMHStep(_)(using _: RandBasis)
          } else {      // if ( mcmcKernel == "Gibbs" ) {
              //mcmc.kernelGibbs(_)(_)
                mcmc.kernelGibbs(_)(using _: RandBasis)
          }
        val ( revInitialMCStates, nNonEquil, nExtendedIterations, independentSamplings) = 
            mcmc.runMC( bmState.initialMCStates,            //initialConfigurations,
                nInitialIterationsPerUnit,
                everyNIterationsPerUnit,
                nSamples,
                initialT,
                annealingRate,
                finalT = finalT,
                maxExtendedIterations = maxExtendedIterations,
                kernel = kernelForMCMC)

        val currentPia = mcmc.frequenciesAtUnitInSamples( independentSamplings )
        val currentPijab = mcmc.pairwiseFrequenciesInSamples( independentSamplings)
        val ensembleAverages = MCMC.EnsembleAverages( pia = currentPia , pijab = currentPijab )

        def leakyAvePiaPijab ( step: Int, beta: Double,
                               leakyAveP: MCMC.EnsembleAverages,
                               currentPia: ArraySeq[DenseVector[Double]], currentPijab: ArraySeq[DenseMatrix[Double]] ) = { 

            val (lAPia, lAPijab) = (leakyAveP.pia, leakyAveP.pijab)
            if ( beta == 0.0 ) 
                (currentPia, currentPijab)
            else if ( lAPia.size == 0 || lAPia(0).size == 0 ) { 
                assert( lAPijab.size == 0 || lAPijab(0).size == 0 )
                assert( step == 1 )
                (currentPia, currentPijab)
            } else {
                (
                  leakyAveVec( beta, lAPia, currentPia) ,
                  leakyAveMat( beta, lAPijab, currentPijab)
                )
            }

        }
        val (pia, pijab)           = leakyAvePiaPijab ( bmState.step + 1 - bmState.stepOffset, 
                                         bmState.betaLAP, bmState.leakyAveP, currentPia, currentPijab ) 

        val (piaForKL, pijabForKL) = leakyAvePiaPijab ( bmState.step + 1 - bmState.stepOffset, 
                                         bmState.betaLAPForKL, bmState.leakyAvePForKL, currentPia, currentPijab ) 

        val (dldphia, dldphijab) =
                if(regTerm == "GL1L2")
                        BM.dLdphiWithGL1L2(propL1h, propL1J, fia, fijab, pia, pijab, 
                                bmState.bmInteractions.phia, bmState.bmInteractions.phijab, lambdaPhi, lambdaPhij )
                else if(regTerm == "L1L2")
                        BM.dLdphiWithL1L2(propL1h, propL1J, fia, fijab, pia, pijab, 
                                bmState.bmInteractions.phia, bmState.bmInteractions.phijab, lambdaPhi, lambdaPhij )
                else if(regTerm == "L1")
                        BM.dLdphiWithL1(fia, fijab, pia, pijab, 
                                bmState.bmInteractions.phia, bmState.bmInteractions.phijab, lambdaPhi, lambdaPhij )
                else if(regTerm == "L2")
                        BM.dLdphiWithL2(fia, fijab, pia, pijab, 
                                bmState.bmInteractions.phia, bmState.bmInteractions.phijab, lambdaPhi, lambdaPhij )
                else {
                        sys.error("Not supported: " + regTerm) 
                        BM.dLdphiWithL2(fia, fijab, pia, pijab, 
                                bmState.bmInteractions.phia, bmState.bmInteractions.phijab, lambdaPhi, lambdaPhij )
                }

        val grad = BM.BMInteractions( dldphia, dldphijab)

        val bmStateWithMCSamples = bmState.copy( grad = grad, initialMCStates = revInitialMCStates,
              independentSamplings = independentSamplings, ensembleAverages = ensembleAverages,
              leakyAveP = if ( bmState.betaLAP <= 0.0 ) MCMC.EnsembleAverages() else MCMC.EnsembleAverages(pia, pijab) ,
              leakyAvePForKL = if ( bmState.betaLAPForKL <= 0.0 ) MCMC.EnsembleAverages() else MCMC.EnsembleAverages(piaForKL, pijabForKL) )

        val oldBMInteractions = bmState.bmInteractions 

        //val (nextBMStateWithoutSoft: BMState, learningRates: BMInteractions) =
        val nextBMStateWithoutSoftPlus =
                if ( optMethod == "NAG" )
                        BM.learningByNAG(bmStateWithMCSamples, grad )
                else if ( optMethod == "ModRAdamMax" || optMethod == "ModRAdam" ) 
                        BM.learningByModRAdamMax(bmStateWithMCSamples, grad )
                else if ( optMethod == "ModAdamMax" || optMethod == "ModAdam" ) 
                        BM.learningByModAdamMax(bmStateWithMCSamples, grad )
                else if ( optMethod == "ModAdamSum" ) 
                        BM.learningByModAdamSum(bmStateWithMCSamples, grad )
                else if ( optMethod == "Adam" ) 
                        BM.learningByAdam(bmStateWithMCSamples, grad )
                else if ( optMethod == "RAdam" ) 
                        BM.learningByRAdam(bmStateWithMCSamples, grad )
                else if ( optMethod == "Adadelta" ) 
                        BM.learningByAdadelta(bmStateWithMCSamples, grad )
                else if ( optMethod == "RPROP-LR" || optMethod == "MF" ) 
                        // The last grad, bmState.grad, is needed.
                        BM.learningByRPROPLR(bmState, grad, regTerm, propL1h, propL1J )
                else {
                        sys.error("Not supported: " + optMethod) 
                        BM.learningByAdam(bmStateWithMCSamples, grad )
                }
        val nextBMStateWithoutSoft = nextBMStateWithoutSoftPlus._1
        val optionLearningRates = nextBMStateWithoutSoftPlus._2

        val nextBMState = if ( noSoftThresholding ) {
                nextBMStateWithoutSoft
            } else {
                if( regTerm == "L2" || ( (regTerm == "L1L2" || regTerm == "GL1L2") && (propL1h == 0.0 && propL1J == 0.0) ) || optionLearningRates == None )
                        nextBMStateWithoutSoft
                else if ( regTerm == "GL1" || regTerm == "GL1L2" )
                        softThresholdingForGL1(oldBMInteractions, optMethod, regTerm, 
                                propL1h, propL1J,
                                lambdaPhi, lambdaPhij,
                                nextBMStateWithoutSoft, optionLearningRates.get )
                else if ( regTerm == "L1" || regTerm == "L1L2" )
                        softThresholdingForL1(oldBMInteractions, optMethod, regTerm,
                                propL1h, propL1J,
                                lambdaPhi, lambdaPhij,
                                nextBMStateWithoutSoft, optionLearningRates.get )
                else {
                        sys.error(s"Not supported: $regTerm")
                        nextBMStateWithoutSoft
                }
            }

      //(mcmc, initialMCStates, nExtendedIterations, bmStateWithMCSamples, nextBMState)
        (mcmc, nNonEquil, nExtendedIterations, bmStateWithMCSamples,
               nextBMState.copy( 
                  leakyAveP = if ( bmState.betaLAP <= 0.0 ) MCMC.EnsembleAverages() else MCMC.EnsembleAverages(pia, pijab) ,
                  leakyAvePForKL = if ( bmState.betaLAPForKL <= 0.0 ) MCMC.EnsembleAverages() else MCMC.EnsembleAverages(piaForKL, pijabForKL) )
        )
    }

    def bayesianCorrection(fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]], observedN: Double, pseudoN: Double ) = {
      require( observedN >= 0.0 && pseudoN >= 0.0 )
      assert( fia(0).size * fia(0).size == fijab(0).size )
      assert( (fia.size * (fia.size -1)) / 2  == fijab.size )

      if ( pseudoN <= 0 ) {
        (fia, fijab)
      } else if ( fia.size <= 0 || fia(0).size <= 0 ) {
        (fia, fijab)
      } else {
        val ratioO = observedN / (observedN + pseudoN)
        val ratioP = 1.0 - ratioO
        val pseudoA = 1.0 / fia(0).size * ratioP
        val pseudoAB = 1.0 / fijab(0).size * ratioP

        ( fia.map  ( fa  => fa  * ratioO + pseudoA  ),
          fijab.map( fab => fab * ratioO + pseudoAB ) )
      }
    }

    def bayesianCorrection(fia: ArraySeq[DenseVector[Double]], observedN: Double, pseudoN: Double ) = {
      require( observedN >= 0.0 && pseudoN >= 0.0 )

      if ( pseudoN <= 0 ) {
        fia
      } else if ( fia.size <= 0 || fia(0).size <= 0 ) {
        fia
      } else {
        val ratioO = observedN / (observedN + pseudoN)
        val ratioP = 1.0 - ratioO
        val pseudoA = 1.0 / fia(0).size * ratioP
        fia.map{ fa => fa * ratioO + pseudoA }
      }
    }

    def bayesianCorrection(fa: DenseVector[Double], observedN: Double, pseudoN: Double ) = {
      require( observedN >= 0.0 && pseudoN >= 0.0 )

      if ( pseudoN <= 0 ) {
        fa
      } else if ( fa.size <= 0 ) {
        fa
      } else {
        val ratioO = observedN / (observedN + pseudoN)
        val ratioP = 1.0 - ratioO
        val pseudoA = 1.0 / fa.size * ratioP
        fa * ratioO + pseudoA
      }
    }

    def replace0By1ForKL (fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]] ) = {
        import breeze.numerics.log
        val pia   = fia.map  ( fi  => fi.map ( a  => if (a  <= 0.0 ) 1.0 else a  ) )
        val pijab = fijab.map( fij => fij.map( ab => if (ab <= 0.0 ) 1.0 else ab ) )
        ( pia, pijab ) 
    }

    def freqsForBMandKL( fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]] , fa: DenseVector[Double],
                effectiveNConfigs: Double,
                pseudoNCountsForBM: Double = 0.0,        // If the regularization is not used, this may be used to avoid 0 counts in this.fia and this.fijab
                pseudoNCountsForKL: Double = 1.0         // used to calculate KL valunes from ensemble averages of Pij(a,b)
        ) = {

        require( pseudoNCountsForBM >= 0.0 && pseudoNCountsForKL > 0.0 )

        val (fiaForBM, fijabForBM, faForBM ) = 
            if ( pseudoNCountsForBM > 0.0 ) {
              val (fiaForBM, fijabForBM) = bayesianCorrection(fia, fijab, effectiveNConfigs, pseudoNCountsForBM)
              val faForBM = bayesianCorrection( fa, effectiveNConfigs * fia.size, pseudoNCountsForBM) 
              (fiaForBM, fijabForBM, faForBM )
            } else {
              (fia, fijab, fa )
            }

        val (fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1) = {
            if ( pseudoNCountsForBM > 0.0 ) {
              (fiaForBM, fijabForBM, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]]() )
            } else {
            // Be carefull
            /*
              if ( pseudoNCountsForKL > 0.0 ) {
                val (fiaForKL, fijabForKL) = bayesianCorrection(fia, fijab, effectiveNConfigs, pseudoNCountsForKL)
                (fiaForKL, fijabForKL, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]]() )
              } else
            */
              {
                val (fiaWith0ReplacedBy1, fijabWith0ReplacedBy1) = replace0By1ForKL(fia, fijab)
                (fia, fijab, fiaWith0ReplacedBy1, fijabWith0ReplacedBy1 )
              }
            }
        }
        ( fiaForBM, fijabForBM, faForBM,  fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1 )
    }

    def calcKL(fiaP: ArraySeq[DenseVector[Double]], fijabP: ArraySeq[DenseMatrix[Double]], observedNP: Double, pseudoNP: Double,
               fiaQ: ArraySeq[DenseVector[Double]], fijabQ: ArraySeq[DenseMatrix[Double]], observedNQ: Double, pseudoNQ: Double,
               fiaPWith0ReplacedBy1: ArraySeq[DenseVector[Double]] = ArraySeq[DenseVector[Double]](),       
               fijabPWith0ReplacedBy1: ArraySeq[DenseMatrix[Double]] = ArraySeq[DenseMatrix[Double]]()      
              ) = {

        import breeze.linalg.{InjectNumericOps, sum}
        import breeze.numerics.log

        val ((piaP_0, piaP_1), (pijabP_0, pijabP_1)) = if (observedNP == scala.Double.PositiveInfinity ) { 
                                ( (fiaP, fiaP), (fijabP, fijabP) )          // means fiaP > 0 and fijabP > 0 
                             } else if ( pseudoNP <= 0.0 ) {
                                val piaP   = if( fiaPWith0ReplacedBy1.size > 0) {
                                     fiaPWith0ReplacedBy1
                                  } else {
                                     fiaP.map ( fi  => fi.map ( a  => if (a  <= 0.0 ) 1.0 else a  ) )
                                  }
                                val pijabP   = if( fijabPWith0ReplacedBy1.size > 0) {
                                     fijabPWith0ReplacedBy1
                                  } else {
                                     fijabP.map( fij => fij.map( ab => if (ab <= 0.0 ) 1.0 else ab ) )
                                  }
                                ( (fiaP, piaP), (fijabP, pijabP) ) 
                             } else {
                                val (piaP, pijabP) =
                                    bayesianCorrection(fiaP, fijabP, observedNP, pseudoNP)
                                ( (piaP, piaP), (pijabP, pijabP) )
                             }
        val ((piaQ_0, piaQ_1), (pijabQ_0, pijabQ_1)) =  if (observedNQ == scala.Double.PositiveInfinity ) {
                                ( (fiaQ, fiaQ), (fijabQ, fijabQ) )          // means fiaQ > 0 and fijabQ > 0
                             } else if ( pseudoNQ <= 0.0 ) {
                                val piaQ   = fiaQ.zipWithIndex.map   {(fi, i)    => fi.pairs.map { (a, f)  => if (piaP_0(i)(a) <= 0.0 && f <= 0.0) 1.0 else f } } 
                                val pijabQ = fijabQ.zipWithIndex.map { (fij, ij) => fij.pairs.map{ (ab, f) => if (pijabP_0(ij)(ab) <= 0.0 && f <= 0.0) 1.0 else f } }  
                                ( (fiaQ, piaQ), (fijabQ, pijabQ) ) 
                             } else {
                                val (piaQ, pijabQ) =
                                    bayesianCorrection(fiaQ, fijabQ, observedNQ, pseudoNQ)
                                ( (piaQ, piaQ), (pijabQ, pijabQ) )
                             }

        val KL1 = piaP_0.zip  (piaP_1.zip(piaQ_1)).map     (pi  => breeze.linalg.sum(pi._1  *:* breeze.numerics.log(pi._2._1  /:/ pi._2._2) ) ).sum
        val KL2 = pijabP_0.zip(pijabP_1.zip(pijabQ_1)).map (pij => breeze.linalg.sum(pij._1 *:* breeze.numerics.log(pij._2._1 /:/ pij._2._2) ) ).sum
        
        val nUnits = fiaP.size
        val nPairs = fijabP.size

        ( KL1 / nUnits, KL2 / nPairs )
    }

    def keepBestKLs(bestKLs: Queue[Option[ Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ], 
                caseKL: Tuple5[Double, Double, Int, Option[File], Option[File] ]) = {
        val nBestKLs = bestKLs.size
        @annotation.tailrec
        def scanQ(i: Int): (Queue[ Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ], Boolean) = {
                if(i < 0) {
                  val files = bestKLs.takeRight(1).head
                  if ( files != None && files.get._4 != None ) files.get._4.get.delete()
                  if ( files != None && files.get._5 != None ) files.get._5.get.delete()
                  (Option(caseKL) +: bestKLs.dropRight(1), true)
                } else {
                        if ( bestKLs(i) == None || caseKL._2 < bestKLs(i).get._2 ) {
                                scanQ(i - 1)
                        } else {
                                if(i >= nBestKLs -1) {
                                        (bestKLs, false)
                                } else {
                                        val files = bestKLs.takeRight(1).head
                                        if ( files != None && files.get._4 != None ) files.get._4.get.delete()
                                        if ( files != None && files.get._5 != None ) files.get._5.get.delete()
                                        val q1 = bestKLs.take(i + 1) 
                                        if ( nBestKLs - i - 1 > 1) {
                                                val q2 = Option(caseKL) +: bestKLs.drop(i +1)
                                                (q1 ++: q2.dropRight(1), true)
                                        } else {
                                                (q1 :+ Option(caseKL), true)
                                        }
                                }
                        }
                }
        }
        scanQ(nBestKLs - 1)
    }

    def deleteBestKLs(bestKLs: Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ] ) = {
        val nBestKLs = bestKLs.size

        val noneQueue = bestKLs.map( i => {
                  if ( i != None && i.get._4 != None ) i.get._4.get.delete()
                  if ( i != None && i.get._5 != None ) i.get._5.get.delete()
                  val none: Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] = None
                  none
                } )
        noneQueue
    }

    def printBMInteractions(outphi: PrintStream, stateOrderString: String, 
                bmState: BMState ) = {

        val phia = bmState.bmInteractions.phia
        val phijab = bmState.bmInteractions.phijab
        MCMC.printhOrPi(outphi, stateOrderString, phia, "phi(i)" )
        MCMC.printJOrPij(outphi, stateOrderString, phijab, "phi(ij)")
    }

    def miniBatchSizeCorrected (miniBatchSize: Int, fullBatchSize: Int) = {
        if ( miniBatchSize <= 0 || miniBatchSize > fullBatchSize )
           fullBatchSize
        else
           miniBatchSize
    }

    def nMiniBatchesInFullBatch(miniBatchSize: Int, fullBatchSize: Int) = {
        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)
        (fullBatchSize.toDouble / miniBSize).floor.toInt
    }

    def startPosInFB(startMBPos: Int, miniBatchSize: Int, fullBatchSize: Int ) = {    //nMBsInFB: Int) = {
        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)
        assert( startMBPos < nMiniBatchesInFullBatch(miniBatchSize, fullBatchSize) )
      //assert( startMBPos < nMBsInFB )
        startMBPos * miniBSize
    }

    def startMBPosInFB(startPos: Int, miniBatchSize: Int, fullBatchSize: Int) = {
        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)
        assert( startPos < endPosOfMiniBatch(startPos, miniBatchSize, fullBatchSize)._3 )
        assert( startPos % miniBSize == 0 )
        startPos / miniBSize
    }

    def endPosOfMiniBatch(startPos: Int, miniBatchSize: Int, fullBatchSize: Int) = {
        // Return (startPos, endPos, if(endPos == 0) fullBatchSize else endPos )

        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)
      //if ( miniBatchSize <= 0 || miniBatchSize > fullBatchSize ) {
      //      assert( startPos == 0 )
      //      (startPos, 0, fullBatchSize)
      //} else
        {
            //require( startPos >= 0 && startPos < fullBatchSize )
            //assert( miniBatchSize <= fullBatchSize )
              require( startPos >= 0 && startPos % miniBSize == 0 )
              val s = startPos
              if ( s + miniBSize * 2 > fullBatchSize) {
                  (s, 0, fullBatchSize)          // last miniBatch size > miniBatchSize
              } else if ( s + miniBSize > fullBatchSize ) {
                  assert( s + miniBSize <= fullBatchSize, "Must be: startPos + miniBatchSize <= fullBatchSize\n" )
                  (s, 0, fullBatchSize)          // dummy
              } else if ( s + miniBSize == fullBatchSize ) {
                  (s, 0, fullBatchSize)          // last miniBatch size <= miniBatchSize
              } else {
                  (s, s + miniBSize, s + miniBSize)
              }
        }
    }

    def sliceMiniBatch[A](rangeExcl: Tuple2[Int, Int], array: ArraySeq[A]) (using classTagA: ClassTag[A]): ArraySeq[A] = {

        val (s, e) = rangeExcl    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= array.size )
          array.slice(s, e) 
        } else { //  ( e - s <= 0 )
          assert( s < array.size )
          val head = array.slice(s, array.size)
          if ( e <= 0 ) {
            head
          } else {                                // This case does not occur in the present endPosOfMiniBatch.
            assert( e < array.size )
            val tail = array.slice(0, e)
          //ArraySeq.concat(head, tail)
            head ++ tail
          }
        }
    }

    def sliceMiniBatch[A](rangeExcl: Tuple2[Int, Int], vector: Vector[A]) (using classTagA: ClassTag[A]): Vector[A] = {
        val array = vector
        val (s, e) = rangeExcl    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= array.size )
          array.slice(s, e) 
        } else { //  ( e - s <= 0 )
          assert( s < array.size )
          val head = array.slice(s, array.size)
          if ( e <= 0 ) {
            head
          } else {                                // This case does not occur in the present endPosOfMiniBatch.
            assert( e < array.size )
            val tail = array.slice(0, e)
          //ArraySeq.concat(head, tail)
            head ++ tail
          }
        }
    }

    def sliceMiniBatch[A](rangeExcl: Tuple2[Int, Int], vector: ParVector[A]) (using classTagA: ClassTag[A]): ParVector[A] = {
        val array = vector
        val (s, e) = rangeExcl    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= array.size )
          array.slice(s, e) 
        } else { //  ( e - s <= 0 )
          assert( s < array.size )
          val head = array.slice(s, array.size)
          if ( e <= 0 ) {
            head
          } else {                                // This case does not occur in the present endPosOfMiniBatch.
            assert( e < array.size )
            val tail = array.slice(0, e)
          //ArraySeq.concat(head, tail)
            head ++ tail
          }
        }
    }

    def sliceMiniBatch[A](range: Tuple2[Int, Int], array: ArraySeq[A], byIndex: IndexedSeq[Int]) (using classTagA: ClassTag[A]): ArraySeq[A] = {

        assert( array.size == byIndex.size )
        val (s, e) = range    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= byIndex.size )
          ArraySeq.range(s, e).map( i => array( byIndex(i) ) )
        } else { //  ( e - s <= 0 )
          assert( s < byIndex.size )
          val head = ArraySeq.range(s, array.size).map( i => array( byIndex(i) ) )
          if ( e <= 0 ) {
            head
          } else {
            assert( e < byIndex.size )
            val tail = ArraySeq.range(0, e).map( i => array( byIndex(i) ) )
          //ArraySeq.concat(head, tail)
            head ++ tail
          }
        }
    }

    def sliceMiniBatch[A](range: Tuple2[Int, Int], vector: Vector[A], byIndex: IndexedSeq[Int]) (using classTagA: ClassTag[A]): Vector[A] = {

        val array = vector
        assert( array.size == byIndex.size )
        val (s, e) = range    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= byIndex.size )
          Vector.range(s, e).map( i => array( byIndex(i) ) )
        } else { //  ( e - s <= 0 )
          assert( s < byIndex.size )
          val head = Vector.range(s, array.size).map( i => array( byIndex(i) ) )
          if ( e <= 0 ) {
            head
          } else {
            assert( e < byIndex.size )
            val tail = Vector.range(0, e).map( i => array( byIndex(i) ) )
          //ArraySeq.concat(head, tail)
            head ++ tail
          }
        }
    }

    def sliceMiniBatch[A](range: Tuple2[Int, Int], vector: ParVector[A], byIndex: IndexedSeq[Int]) (using classTagA: ClassTag[A]): ParVector[A] = {

        sliceMiniBatch(range, vector.toVector, byIndex).par
      /*
        val array = vector
        assert( array.size == byIndex.size )
        val (s, e) = range    // e == 0 means e == array.size
        if ( e - s > 0 ) {
          assert( e <= byIndex.size )
          Vector.range(s, e).map( i => array( byIndex(i) ) ).par
        } else { //  ( e - s <= 0 )
          assert( s < byIndex.size )
          val head = Vector.range(s, array.size).map( i => array( byIndex(i) ) )
          if ( e <= 0 ) {
            head.par
          } else {
            assert( e < byIndex.size )
            val tail = Vector.range(0, e).map( i => array( byIndex(i) ) )
          //ArraySeq.concat(head, tail)
            (head ++ tail).par
          }
        }
      */
    }

  //def copyToSlice[A](rangeExcl: Tuple2[Int, Int], to_vector: Vector[A], vector: Vector[A]) (using classTagA: ClassTag[A]): Vector[A] = {
    def copyToSlice[A](rangeExcl: Tuple2[Int, Int], to_vector: Vector[A], vector: Vector[A]) : Vector[A] = {
        // 0 <= endPos <= to_vector.size 

          val toArraySeq = to_vector
          val array = vector

          val (startPos, endPos) = rangeExcl
          require(endPos >= 0 && endPos <= toArraySeq.size)
          require(startPos >= 0 && startPos < toArraySeq.size)

          if( startPos == endPos && array.size == 0 ) {
            toArraySeq
          } else {
            val newArraySeq =
              if ( endPos - startPos > 0 ) {
                assert( endPos - startPos == array.size )
                val part0 = toArraySeq.slice(0, startPos) 
                val part1 = array
                if ( endPos  == toArraySeq.size ) {
                    part0 ++ part1
                } else {
                    val part2 = toArraySeq.slice(endPos, toArraySeq.size) 
                    part0 ++ part1 ++ part2
                }
              } else { //  ( endPos - startPos <=  0 )
                if ( startPos == 0 && endPos == 0 ) {
                    array
                } else if ( endPos == 0 ) {
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array
                    part1 ++ part2
                } else {                                             // This case does not occur in the present endPosOfMiniBatch.
                    assert( endPos < startPos )
                    val part0 = array.slice(toArraySeq.size - startPos, array.size)
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array.slice(0, toArraySeq.size - startPos)
                    part0 ++ part1 ++ part2
                }
              }
            assert ( newArraySeq.size == toArraySeq.size)
            newArraySeq
          }
    }

  //def copyToSlice[A](rangeExcl: Tuple2[Int, Int], to_vector: Vector[A], vector: Vector[A]) (using classTagA: ClassTag[A]): Vector[A] = {
    def copyToSlice[A](rangeExcl: Tuple2[Int, Int], to_vector: ParVector[A], vector: ParVector[A]) : ParVector[A] = {
        // 0 <= endPos <= to_vector.size 

          val toArraySeq = to_vector
          val array = vector

          val (startPos, endPos) = rangeExcl
          require(endPos >= 0 && endPos <= toArraySeq.size)
          require(startPos >= 0 && startPos < toArraySeq.size)

          if( startPos == endPos && array.size == 0 ) {
            toArraySeq
          } else {
            val newArraySeq =
              if ( endPos - startPos > 0 ) {
                assert( endPos - startPos == array.size )
                val part0 = toArraySeq.slice(0, startPos) 
                val part1 = array
                if ( endPos  == toArraySeq.size ) {
                    part0 ++ part1
                } else {
                    val part2 = toArraySeq.slice(endPos, toArraySeq.size) 
                    part0 ++ part1 ++ part2
                }
              } else { //  ( endPos - startPos <=  0 )
                if ( startPos == 0 && endPos == 0 ) {
                    array
                } else if ( endPos == 0 ) {
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array
                    part1 ++ part2
                } else {                                             // This case does not occur in the present endPosOfMiniBatch.
                    assert( endPos < startPos )
                    val part0 = array.slice(toArraySeq.size - startPos, array.size)
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array.slice(0, toArraySeq.size - startPos)
                    part0 ++ part1 ++ part2
                }
              }
            assert ( newArraySeq.size == toArraySeq.size)
            newArraySeq
          }
    }

    def copyToSlice[A](rangeExcl: Tuple2[Int, Int], to_array: ArraySeq[A], array: ArraySeq[A]) (using classTagA: ClassTag[A]): ArraySeq[A] = {
        // 0 <= endPos <= toArraySeq.size 

          val toArraySeq = to_array

          val (startPos, endPos) = rangeExcl
          require(endPos >= 0 && endPos <= toArraySeq.size)
          require(startPos >= 0 && startPos < toArraySeq.size)

          if( startPos == endPos && array.size == 0 ) {
            toArraySeq
          } else {
            val newArraySeq =
              if ( endPos - startPos > 0 ) {
                assert( endPos - startPos == array.size )
                val part0 = toArraySeq.slice(0, startPos) 
                val part1 = array
                if ( endPos  == toArraySeq.size ) {
                  //ArraySeq.concat(part0, part1)
                    part0 ++ part1
                } else {
                    val part2 = toArraySeq.slice(endPos, toArraySeq.size) 
                  //ArraySeq.concat(part0, part1, part2)
                    part0 ++ part1 ++ part2
                }
              } else { //  ( endPos - startPos <=  0 )
                if ( startPos == 0 && endPos == 0 ) {
                    array
                } else if ( endPos == 0 ) {
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array
                  //ArraySeq.concat(part1, part2)
                    part1 ++ part2
                } else {                                             // This case does not occur in the present endPosOfMiniBatch.
                    assert( endPos < startPos )
                    val part0 = array.slice(toArraySeq.size - startPos, array.size)
                    val part1 = toArraySeq.slice(endPos, startPos)
                    val part2 = array.slice(0, toArraySeq.size - startPos)
                  //ArraySeq.concat(part0, part1, part2)
                    part0 ++ part1 ++ part2
                }
              }
            assert ( newArraySeq.size == toArraySeq.size)
            newArraySeq
          }
    }

//def rep[T, C <: Seq[T]](l: C & scala.collection.SeqOps[T, scala.collection.Seq, C], n: Int): C = if (n <=1) l else (l ++ rep(l,n-1)).asInstanceOf[C]

  //def randomizeConfigurations[T,C <: immutable.IndexedSeq[T]](
  //    configurations: C & scala.collection.SeqOps[T, scala.collection.Seq, C]): C = {
  //def randomizeConfigurations[C <: immutable.IndexedSeq[MCMC.Configuration]](
  //    configurations: C & scala.collection.SeqOps[MCMC.Configuration, scala.collection.Seq, C]): C = {
  //def randomizeConfigurations[C <: immutable.IndexedSeq[MCMC.Configuration]](
  //    configurations: C )(using classTagC: ClassTag[C]): C = {
    def randomizeConfigurations( configurations: Vector[MCMC.Configuration] ): Vector[MCMC.Configuration] = {
          configurations.map{ s =>
            //val randomOrder = Random.shuffle( IArray.range(0, s.size) )
            //IArray(randomOrder.map( i => s(i) )*) 
              IArray( Random.shuffle( s )* )
            }
    }

    def randomizeConfigurations( configurations: ParVector[MCMC.Configuration] ): ParVector[MCMC.Configuration] = {
          // To reproduce the same result
          val seeds = Vector.range(0, configurations.size).map{ _ => Random.nextInt() }
          configurations.zip(seeds).map( (s, r) => 
            {
              val random = new Random( r )
            //val randomOrder = random.shuffle( IArray.range(0, s.size) )
            //IArray(randomOrder.map( i => s(i) )*) 
              IArray( random.shuffle( s )* )
            } )
    }

    def nMiniBatchesForEnsembleAverage(miniBatchSize: Int, maxSize: Int, nMiniBatchesForEnsembleAve: Int ) = { 

        assert(maxSize > 0)
        val miniBS = miniBatchSizeCorrected (miniBatchSize, maxSize)
        val nMBForAve = //if ( miniBatchSize <= 0 || miniBatchSize >= maxSize ) {
                        //  1
                        //} else 
                        if ( nMiniBatchesForEnsembleAve <= 0 || miniBS * nMiniBatchesForEnsembleAve > maxSize ) {
                            math.max( (maxSize / miniBS.toDouble ).floor.toInt, 1)
                        } else {
                            nMiniBatchesForEnsembleAve
                        }
        nMBForAve
    }

  /*
    def rangeExclForEnsembleAverage(independentSamplingsSize: Int, fullBatchSize: Int, miniBatchSize: Int, 
                miniBatchRangeExcl: Tuple2[Int, Int], nMiniBatchesForEnsembleAve: Int ) = {

        require(independentSamplingsSize >= fullBatchSize)
        val (startPos, untilEndPos) = miniBatchRangeExcl
        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)  
        assert ( startPos >= 0 && startPos % miniBSize == 0 && untilEndPos <= fullBatchSize ) 
        assert ( untilEndPos - startPos >= miniBSize ) 

        val nMBperFB = nMiniBatchesInFullBatch(miniBSize, fullBatchSize)
      //val lastMBSize = (fullBatchSize - nMBperFB * miniBSize) + miniBSize
        
        val nFB = (independentSamplingsSize.toDouble / fullBatchSize).ceil.toInt

        val miniBPosInFullB = Vector.range(0, nMBperFB).map( i => i * miniBSize )

        val miniBPosInNFB = Vector.range(0, nFB + 1).map( i => {
                               if ( i < nFB ) {
                                    miniBPosInFullB.map(j => i * fullBatchSize + j )
                               } else {
                                    Vector[Int](i * fullBatchSize)
                               }
                            } ).flatten

        assert( miniBPosInNFB.size == nFB * nMBperFB + 1 )

        val unSuppliedSize = (nFB * fullBatchSize - independentSamplingsSize)
        assert( unSuppliedSize % miniBSize == 0 )
        val unSuppliedNMB = unSuppliedSize / miniBSize

        val unUsedNMBAtTail = nMBperFB - ((startPos / miniBSize) + 1)
        val usedNMBAtTail = math.min((startPos / miniBSize) + 1 , nMiniBatchesForEnsembleAve )
        val unUsedHeadNMBAtTail = nMBperFB - usedNMBAtTail - unUsedNMBAtTail

        assert(nMiniBatchesForEnsembleAve >= 0)
        val usedNMBAtNonTail = math.min(nFB * nMBperFB - unSuppliedNMB - nMBperFB, nMiniBatchesForEnsembleAve - usedNMBAtTail)

        val totalNMBAtNonTail = nFB * nMBperFB - unSuppliedNMB - nMBperFB
        assert(totalNMBAtNonTail >= 0)

        val unUsedNMBAtNonTail = ( totalNMBAtNonTail - usedNMBAtNonTail )

        assert( unSuppliedNMB + (unUsedNMBAtNonTail + usedNMBAtNonTail) +
                (unUsedHeadNMBAtTail + usedNMBAtTail) + unUsedNMBAtTail == nFB * nMBperFB  )
         
        val fromPosNMBForEnsembleAve = if ( usedNMBAtNonTail > 0 )
                             unSuppliedNMB + unUsedNMBAtNonTail
                         else 
                             unSuppliedNMB + unUsedNMBAtNonTail + unUsedHeadNMBAtTail
        val untilPosNMBForEnsembleAve = unSuppliedNMB + (unUsedNMBAtNonTail + usedNMBAtNonTail) + (unUsedHeadNMBAtTail + usedNMBAtTail)

        val fromPosForEnsembleAve = miniBPosInNFB(fromPosNMBForEnsembleAve) - miniBPosInNFB(unSuppliedNMB)
        val untilPosForEnsembleAve = miniBPosInNFB(untilPosNMBForEnsembleAve) - miniBPosInNFB(unSuppliedNMB)

        assert( untilPosForEnsembleAve == (independentSamplingsSize - fullBatchSize) + untilEndPos )

        val miniBPosInRange = miniBPosInNFB.slice(fromPosNMBForEnsembleAve, untilPosNMBForEnsembleAve + 1).map( p => p - miniBPosInNFB(fromPosNMBForEnsembleAve) ) 

        assert( untilPosNMBForEnsembleAve - fromPosNMBForEnsembleAve + 1 == miniBPosInRange.size )
        assert( miniBPosInRange(0) == 0 && miniBPosInRange(miniBPosInRange.size - 1) == untilPosForEnsembleAve - fromPosForEnsembleAve ) 
        assert( miniBPosInRange(miniBPosInRange.size - 1) - miniBPosInRange(miniBPosInRange.size - 2 ) == untilEndPos - startPos ) 

        (fromPosForEnsembleAve, untilPosForEnsembleAve, miniBPosInRange )
    }
  */

    def miniBatchRangeExclForEnsembleAverage(nMBinPrevMiniBatches: Int, nMBinFullB: Int,
        startMBPosInFull: Int, nMiniBatchesForEnsembleAveCorrected: Int ) = {    // startMBPos == startPos / miniBatchSize  

        assert( nMiniBatchesForEnsembleAveCorrected > 0 )
        assert( nMBinFullB >= startMBPosInFull + 1 )
        val untilEndMBPos = nMBinPrevMiniBatches + startMBPosInFull + 1
        val startMBPos = math.max (untilEndMBPos - nMiniBatchesForEnsembleAveCorrected, 0)

        (startMBPos, untilEndMBPos)
    }

    def printKL(outKL: PrintStream, bmStep: Int, nNonEquil: Int, nExtendedIterations: Int, 
                ensembleKL: Tuple2[Double, Double], aveKL: Tuple2[Double, Double], aveKL2: Tuple2[Double, Double] = (0.0, 0.0) ) = {
            //if ( bmStep < 0 ) { 
            //  outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and all samples with pseudocounts = %g\n".format(pseudoNCounts))
            //} else {
                if ( aveKL2._1 == 0.0 && aveKL2._2 == 0.0 ) 
                  outKL.print("%d\t%d %d\t%g %g\t%g %g\n".format(bmStep, nNonEquil, nExtendedIterations, ensembleKL._1, ensembleKL._2, aveKL._1, aveKL._2))
                else
                  outKL.print("%d\t%d %d\t%g %g\t%g %g\t%g %g\n".format(bmStep, nNonEquil, nExtendedIterations,
                                                                        ensembleKL._1, ensembleKL._2, aveKL._1, aveKL._2, aveKL2._1, aveKL2._2))
            //}
    }

    def printLog( bmStep: Int ,
                  logFiles: Option[Tuple2[Option[java.io.File], Option[java.io.File]]],
                  stateOrderString: String,
                  interactions: MCMC.Interactions, miniBacthesForEnsembleAve: Vector[BMState],
                  lastAndCurrentFullBatch: Vector[BMState], outTE: PrintStream, fa: DenseVector[Double],
                  startPos: Int, miniBatchSize: Int, nMiniBatchesForEnsembleAve: Int, printRawInteractions: Boolean = true ) = {

      assert( bmStep == miniBacthesForEnsembleAve( miniBacthesForEnsembleAve.size - 1).step )
      val fullBatchSize = lastAndCurrentFullBatch.map{bms => math.max(bms.initialMCStates.size, bms.independentSamplings.size) }.sum
      val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)  

      val nMBsInFullB = nMiniBatchesInFullBatch(miniBSize, fullBatchSize )
      assert( nMBsInFullB == lastAndCurrentFullBatch.size)

      val (_, endPos, untilEndPos) = endPosOfMiniBatch( startPos, miniBSize, fullBatchSize )
      val startMBPos = startMBPosInFB(startPos, miniBSize, fullBatchSize) 

      assert( bmStep == lastAndCurrentFullBatch(startMBPos).step )

      val (interactions_Ising) = hJ_in_IsingGauge(interactions)

      val vectIndependentSamplings: Vector[MCMC.IndependentSamplings] =
            miniBacthesForEnsembleAve.map{ bms => bms.independentSamplings }

      val vectIndependentSamplings_Ising = vectIndependentSamplings.map{ independentSamplings =>
         val ( _ , independentSamplings_Ising) = hJ_and_TE_in_IsingGauge(interactions_Ising, independentSamplings)
         independentSamplings_Ising
      }

      if (logFiles != None) {

      //val logFileInteractions = new File(ioFiles("outputDir"), logFiles._1)
      //val logFileIndpendentMCsamplings = new File(ioFiles("outputDir"), logFiles._2)
        val logFileInteractions = logFiles.get._1
        val logFileIndpendentMCsamplings = logFiles.get._2

        if( logFileInteractions != None ) {
          val outhJ = new PrintStream( gzipOutputStream( logFileInteractions.get ) )

          outhJ.println("# Step: " + bmStep.toString )

          if ( printRawInteractions ) 
              MCMC.printInteractions(outhJ, stateOrderString, interactions )
          else
              MCMC.printInteractions(outhJ, stateOrderString, interactions_Ising )

          outhJ.close()
        }

        if( logFileIndpendentMCsamplings != None ) {
          val outSamples = new PrintStream( gzipOutputStream( logFileIndpendentMCsamplings.get ) )

        //if ( ( (bmStep + 1) % nMBsInFullB == 0 || bmStep == 0 )  ) {
          if ( ( (bmStep + 1) % nMBsInFullB == 0 || bmStep == 0 ) && miniBacthesForEnsembleAve.size < lastAndCurrentFullBatch.size ) {
              /* The fullBatch is printed. */
              val fullBatchSamples_Ising =  
                if ( bmStep == 0 ) {
                    vectIndependentSamplings_Ising
                } else if ( miniBacthesForEnsembleAve.size < lastAndCurrentFullBatch.size ) {
                    assert( startMBPos + 1 == nMBsInFullB )
                    val partSize = lastAndCurrentFullBatch.size - miniBacthesForEnsembleAve.size
                    val part: Vector[MCMC.IndependentSamplings] =
                               lastAndCurrentFullBatch.slice(0, partSize).map{ bms => bms.independentSamplings }
                    val part_Ising = part.map{ independentSamplings =>
                      val ( _ , samples_Ising) = hJ_and_TE_in_IsingGauge(interactions_Ising, independentSamplings)
                      samples_Ising
                    }
                    val fullBSamplesIsing = part_Ising ++ vectIndependentSamplings_Ising
                    fullBSamplesIsing
                } else {
                    assert( startMBPos + 1 == nMBsInFullB )
                    val len = vectIndependentSamplings_Ising.size
                    vectIndependentSamplings_Ising.slice(len - nMBsInFullB, len)
                }
              
              cfor(0)( m => m < startMBPos, m => m + 1)( m => {
                    val step = bmStep - ( startMBPos - m )
                  //assert( step >= 0 )
                    outSamples.print(s"# Step: ${step}  The total energies are calculated with the interactions of step= ${bmStep} .\n")
                    MCMC.printIndependentMCsamplings(outSamples, stateOrderString, fullBatchSamples_Ising(m) )
              } )
              outSamples.print(s"# Step: ${bmStep}\n")
              MCMC.printIndependentMCsamplings(outSamples, stateOrderString, fullBatchSamples_Ising(startMBPos) )
          } else {
              /* In this case, the configuratios employed for ensemble averages are printed, but will not be used. */ 
              cfor(0)( m => m < vectIndependentSamplings_Ising.size - 1, m => m + 1)( m => {
                    val step = bmStep - ( vectIndependentSamplings_Ising.size - 1 - m )
                  //assert( step >= 0 )
                    assert( step == miniBacthesForEnsembleAve(m).step )
                    outSamples.print(s"# Step: ${step}  The total energies are calculated with the interactions of step= ${bmStep} .\n")
                    MCMC.printIndependentMCsamplings(outSamples, stateOrderString, vectIndependentSamplings_Ising(m) )
              } )
              outSamples.print(s"# Step: ${bmStep}\n")
              val m = vectIndependentSamplings_Ising.size - 1
              MCMC.printIndependentMCsamplings(outSamples, stateOrderString, vectIndependentSamplings_Ising(m) )
          }

          outSamples.close()
        }
      }

        /**************************************/

      def TE_of_sampleConfigurations( interactions: MCMC.Interactions, 
        configurations: ParVector[IArray[Byte]] ): Tuple2[ArraySeq[Double], ArraySeq[Double]]  = {

          val interactions_Ising = hJ_in_IsingGauge(interactions)
          val confsTE_Ising = MCMC.totalE(configurations, interactions_Ising)

          val randomizedConfs = randomizeConfigurations( configurations )
          val randomizedConfsTE_Ising = MCMC.totalE(randomizedConfs, interactions_Ising)

          (confsTE_Ising, randomizedConfsTE_Ising)
      }
  /*
    //def TE_of_sampleConfigurations( interactions: MCMC.Interactions, 
    //    configurations: ArraySeq[IArray[Byte]] ): Tuple2[ArraySeq[Double], ArraySeq[Double]]  = {
      def TE_of_sampleConfigurations[C <: immutable.Seq[IArray[Byte]]]( interactions: MCMC.Interactions,
        configurations: C )(using classTagA: ClassTag[C]): Tuple2[ArraySeq[Double], ArraySeq[Double]]  = {

          val confs = ArraySeq(configurations* )

          val (interactions_Ising, confsTE_Ising) =
            hJ_and_TE_in_IsingGauge(interactions, confs)

          val randomizedConfs = confs.map( s => 
            {
            //val randomOrder = Random.shuffle( IArray.range(0, s.size) )
            //randomOrder.map( i => s(i) ) 
              IArray( Random.shuffle( s )* )
            } )
          val (interactions_Ising_1, randomizedConfsTE_Ising) =
            hJ_and_TE_in_IsingGauge(interactions_Ising, randomizedConfs)

          (confsTE_Ising, randomizedConfsTE_Ising)
      }
  */

        /**************************************/

    /*
      if ( outTE == null ) {

          (interactions_Ising, independentSamplings_Ising, 
           meanAndVariance.MeanAndVariance(999.9, 999.9, 0), 
           meanAndVariance.MeanAndVariance(999.9, 999.9, 0),
           meanAndVariance.MeanAndVariance(999.9, 999.9, 0),
           meanAndVariance.MeanAndVariance(999.9, 999.9, 0),
           meanAndVariance.MeanAndVariance(999.9, 999.9, 0),
           (999.9, 999.9, 999.9) )

      } else {
      */
      {
          val allEnsembleTE_Ising: Vector[ParVector[ArraySeq[Double]]] = 
                vectIndependentSamplings_Ising.map( bmstate => bmstate.map(mc => mc.map(s => s.energy)) )
          val allEnsembleTE_MeanAndVariance = meanAndVariance(allEnsembleTE_Ising.flatten.flatten )

          val (ensembleTE_Ising, ensembleTE_MeanAndVariance) =
              {
                //val ensembleTE_Ising = vectIndependentSamplings_Ising( vectIndependentSamplings_Ising.size - 1 ).map( s => s.energy ).flatten.toVector
                  val ensembleTE_Ising: ParVector[ArraySeq[Double]] = allEnsembleTE_Ising( allEnsembleTE_Ising.size - 1 )
                  val ensembleTE_MeanAndVariance = meanAndVariance(ensembleTE_Ising.flatten.toVector )
                  ( ensembleTE_Ising, ensembleTE_MeanAndVariance )
              }

        //stderr.print(s"allEnsembleTE_Ising.flatten.flatten.size=${allEnsembleTE_Ising.flatten.flatten.size}\n")
        //stderr.print(s"ensembleTE_Ising.flatten.size=${ensembleTE_Ising.flatten.size}\n")
        //if (allEnsembleTE_Ising.size == nMiniBatchesForEnsembleAve)
        //    sys.error(s"nMiniBatchesForEnsembleAve=${nMiniBatchesForEnsembleAve}\n")

          val (sampleTE_MeanAndVariance, randomizedSampleTE_MeanAndVariance, 
               allRandomizedSampleTE_MeanAndVariance,
               allSampleTE_MeanAndVariance, randomSampleTE_MeanAndVariance) = {

          //val allSampleConfigurations: Vector[ParVector[MCMC.Configuration]] = vectIndependentSamplings.map{ bms => bms.map(mc => mc(mc.size -1).nativeConfig.sample ) }
          //val allSampleConfigurations: Vector[ParVector[MCMC.Configuration]] = lastAndCurrentFullBatch.map{ bms => bms.map(mc => mc(mc.size -1).nativeConfig.sample ) }
          //val allSampleConfigurations: Vector[ParVector[ArraySeq[MCMC.Configuration]]] = lastAndCurrentFullBatch.map{ bms =>
          //                             bms.initialMCStates.map{s => ArraySeq.fill(nSamples)(s.nativeConfig.sample) }}

            val allSampleConfigurations: Vector[ParVector[MCMC.Configuration]] = lastAndCurrentFullBatch.map{ bms =>
                                             if ( bms.initialMCStates.size > 0 )
                                                 bms.initialMCStates.map{s => s.nativeConfig.sample }
                                             else
                                                 bms.independentSamplings.map{mc => mc(0).nativeConfig.sample }
                                         }

            assert(allSampleConfigurations.flatten.size == fullBatchSize )

          //val startPosInAllSamples = startPos * nSamples
          //val untilEndPosInAllSamples = untilEndPos * nSamples

            if ( allSampleConfigurations.size == 0 || allSampleConfigurations(0).size == 0) {
              ( meanAndVariance.MeanAndVariance(999.9, 999.9, 0),   // Dummy
                meanAndVariance.MeanAndVariance(999.9, 999.9, 0),   // Dummy
                meanAndVariance.MeanAndVariance(999.9, 999.9, 0),   // Dummy
                meanAndVariance.MeanAndVariance(999.9, 999.9, 0),   // Dummy
                energyOfRandomSeq( fa , interactions_Ising.hia, interactions_Ising.jijab )  )
            } else {

            //val (allSampleConfsTE_Ising, allRandomizedSampleConfsTE_Ising) =
            //            TE_of_sampleConfigurations( interactions_Ising, allSampleConfigurations.flatten.par )
              val allSampleConfsTE_Ising = MCMC.totalE(allSampleConfigurations.flatten , interactions_Ising) 
              val allSampleTE_M_V = meanAndVariance(allSampleConfsTE_Ising)

              val nSamples = lastAndCurrentFullBatch(0).independentSamplings(0).size

              val totalNSamples = (nSamples.toDouble * nMiniBatchesForEnsembleAve / nMBsInFullB ).ceil.toInt    
              val allSampleConfsForRandomized: Vector[MCMC.Configuration] =
                      randomizeConfigurations(
                        allSampleConfigurations.map{ mc => mc.map{ s => ArraySeq.fill(totalNSamples)(s) } }.flatten.flatten
                      )
              val allRandomizedSampleConfsTE_Ising = MCMC.totalE(allSampleConfsForRandomized , interactions_Ising)
              val allRandomizedSampleTE_M_V = meanAndVariance(allRandomizedSampleConfsTE_Ising)

            //stderr.print(s"allSampleConfsTE_Ising.size=${allSampleConfsTE_Ising.size}\n")
            //stderr.print(s"allRandomizedSampleConfsTE_Ising.size=${allRandomizedSampleConfsTE_Ising.size}\n")

              val (sampleConfsTE_Ising, randomizedSampleConfsTE_Ising, sampleTE_M_V, randomizedSampleTE_M_V) = {
                  {
                      val startMBPos = startMBPosInFB(startPos, miniBSize, fullBatchSize) 
                      val sampleConfigurations: ParVector[MCMC.Configuration] =
                        //lastAndCurrentFullBatch(startMBPos).independentSamplings.map{mc => mc(mc.size -1).nativeConfig.sample }
                          lastAndCurrentFullBatch(startMBPos).independentSamplings.map{mc => mc(0).nativeConfig.sample }

                      assert( lastAndCurrentFullBatch(startMBPos).step >= 0 && 
                              ( startMBPos + 1 >= nMBsInFullB  || lastAndCurrentFullBatch(startMBPos + 1).step < 0 ) ) 

                      assert( sampleConfigurations.size == untilEndPos - startPos )

                      val sampleConfsTE_Ising = MCMC.totalE(sampleConfigurations , interactions_Ising)

                      val sampleConfsForRandomized: ParVector[MCMC.Configuration] = 
                              randomizeConfigurations(  sampleConfigurations.map{ s => ArraySeq.fill(nSamples)(s) }.flatten  )
                      val randomizedSampleConfsTE_Ising = MCMC.totalE(sampleConfsForRandomized , interactions_Ising)

                    //val (sampleConfsTE_Ising, randomizedSampleConfsTE_Ising) =
                    //    TE_of_sampleConfigurations( interactions_Ising, sampleConfigurations )

                    //stderr.print(s"sampleConfsTE_Ising.size=${sampleConfsTE_Ising.size}\n")
                    //stderr.print(s"randomizedSampleConfsTE_Ising.size=${randomizedSampleConfsTE_Ising.size}\n")

                      val sampleTE_M_V = meanAndVariance(sampleConfsTE_Ising)
                      val randomizedSampleTE_M_V = meanAndVariance(randomizedSampleConfsTE_Ising)

                      (sampleConfsTE_Ising, randomizedSampleConfsTE_Ising, sampleTE_M_V, randomizedSampleTE_M_V)
                  }
              }

              val randomSampleTE_M_V = energyOfRandomSeq( fa , interactions_Ising.hia, interactions_Ising.jijab )

              ( sampleTE_M_V, randomizedSampleTE_M_V, allRandomizedSampleTE_M_V, allSampleTE_M_V, randomSampleTE_M_V )
            }
          }

          val nUnits = vectIndependentSamplings_Ising(0)(0)(0).configuration.size

          if ( outTE != null && bmStep <= 0 ) 
            outTE.print("# step\tL=%d\t(all randomized sample mean - all randomized sample variance/L)  TE_/L  <TE>/L\tall randomized sample mean\tall randomized sample variance/L  TE sample variance/L  <(TE-<TE>)^2>/L\n".format(nUnits) )

          if ( outTE != null ) 
              outTE.print( "%d\t%f %f %f\t\t\t%f\t%f %f %f         \tA\n".format( bmStep,
                      (allRandomizedSampleTE_MeanAndVariance.mean - allRandomizedSampleTE_MeanAndVariance.variance) / nUnits,
                      sampleTE_MeanAndVariance.mean / nUnits, ensembleTE_MeanAndVariance.mean / nUnits,
                      allRandomizedSampleTE_MeanAndVariance.mean / nUnits,
                      allRandomizedSampleTE_MeanAndVariance.variance / nUnits,
                      sampleTE_MeanAndVariance.variance / nUnits, ensembleTE_MeanAndVariance.variance / nUnits ) )

          if ( outTE != null && bmStep <= 0 ) 
            outTE.print("# step\tL=%d\t(random sample mean - random sample variance)/L  TE_all/L  <TE>/L  <TE>_m/L\trandom sample mean\trandom sample variance/L  TE all sample variance/ L  <(TE-<TE>)^2>/L  <(TE-<TE>)^2>_m/L ; <...>_m over max %d minibatch ensembles\n".format(nUnits, nMiniBatchesForEnsembleAve) )

          if ( outTE != null ) 
              outTE.print( "%d\t%f %f %f %f\t\t%f\t%f %f %f %f\tB\n".format( bmStep,
                      (randomSampleTE_MeanAndVariance._2 - randomSampleTE_MeanAndVariance._3) / nUnits,
                      allSampleTE_MeanAndVariance.mean / nUnits, 
                      ensembleTE_MeanAndVariance.mean / nUnits, allEnsembleTE_MeanAndVariance.mean / nUnits,
                      randomSampleTE_MeanAndVariance._2 / nUnits,
                      randomSampleTE_MeanAndVariance._3 / nUnits,
                      allSampleTE_MeanAndVariance.variance / nUnits,
                      ensembleTE_MeanAndVariance.variance / nUnits, allEnsembleTE_MeanAndVariance.variance / nUnits ) )

          (interactions_Ising, vectIndependentSamplings_Ising,
           ensembleTE_MeanAndVariance,       sampleTE_MeanAndVariance, allRandomizedSampleTE_MeanAndVariance,
           allEnsembleTE_MeanAndVariance, allSampleTE_MeanAndVariance,        randomSampleTE_MeanAndVariance )
      }
    }

  } // End of object BM


  class BM ( 
    val ioFiles: HashMap[String, File],
    val stateOrderString: String,
    val effectiveNConfigs: Double,
    val fia: ArraySeq[DenseVector[Double]], val fijab: ArraySeq[DenseMatrix[Double]], 
    val regTerm: String,
    val propL1h: Double, val propL1J: Double
  //val lambdaPhi: Double, val lambdaPhij: Double,
  //val learningrate: Double,
  //val betaV: Double = 0.999, val betaM: Double = 0.9, val eps: Double = 1.0e-6,
  //val betaLAP: Double = 0.0
                                     ) {

    //val betaV = 0.999, val betaM = 0.9, val learningRate = 0.000001 to 10, val eps = 1.0e-6 ) = {

    import miyazawa.potts.BM
    import miyazawa.potts.BM.*

    if ( regTerm == "L1L2" || regTerm == "GL1L2" ) {
                assert( propL1h >= 0.0 && propL1h <= 1.0)       // propL1h means the proportion of L1 even for GL1L2.
                assert( propL1J >= 0.0 && propL1J <= 1.0)       // propL1J means the proportion of L1 for L1L2 and that of GL1 for GL1L2.
    } else if ( regTerm == "L2" ) {
       assert( propL1h == 0.0 && propL1J == 0.0 )
    } else {
      require( regTerm == "L2" || regTerm == "GL1L2" || regTerm == "L1L2" )
    }
    val nUnits = fia.size
    val nStatesOfUnit = fia(0).size
    val nPairs = (nUnits * (nUnits - 1)) / 2
    require( fijab.size == nPairs )
    require( fijab(0).size == nStatesOfUnit * nStatesOfUnit)

    val fa = toProb(fia)
  //require( lambdaPhi >= 0.0 && lambdaPhij >= 0.0 )

  /*
    def this( ioFiles: HashMap[String, File],
        stateOrderString: String,
        effectiveNConfigs: Double, fia: ArraySeq[ArraySeq[Double]], fijab: ArraySeq[ArraySeq[ArraySeq[Double]]], 
        regTerm: String, 
        propL1h: Double, propL1J: Double ) = {
      //lambdaPhi: Double, lambdaPhij: Double ) = {
      //betaV: Double, betaM: Double, eps: Double ) = {
                this( ioFiles, stateOrderString, effectiveNConfigs, fia.map(fi => DenseVector(fi)),
                        //fijab.map( fij => DenseMatrix(fij:_*) ),                      // for scala 2.11
                          fijab.map( fij => DenseMatrix(fij*) ),                        // for scala 2.13
                        //fijab.map( fij => new DenseMatrix(fij(0).size, fij.size, fij.flatten.toArray).t ),    // for scala 2.13
                                regTerm, propL1h, propL1J )    //, lambdaPhi, lambdaPhij, betaV, betaM, eps )
    }
  */

  /*
    def this( ioFiles: HashMap[String, File],
        stateOrderString: String,
        effectiveNConfigs: Double, fia: ArraySeq[DenseVector[Double]], fijab: ArraySeq[DenseMatrix[Double]], 
        regTerm: String, 
        propL1h: Double, propL1J: Double ) = {
        betaV: Double, betaM: Double, eps: Double ) = {
                this( ioFiles, stateOrderString, effectiveNConfigs, fia, fijab,
                                regTerm, propL1h, propL1J, 0.0, 0.0, betaV, betaM, eps )
    }
  */

  /*
    def this( ioFiles: HashMap[String, File],
        stateOrderString: String,
        effectiveNConfigs: Double, fia: ArraySeq[ArraySeq[Double]], fijab: ArraySeq[ArraySeq[ArraySeq[Double]]], 
        regTerm: String, 
        propL1h: Double, propL1J: Double ) = {
      //betaV: Double, betaM: Double, eps: Double ) = {
                this( ioFiles, stateOrderString, effectiveNConfigs, fia.map(fi => DenseVector(fi)),
                        //fijab.map( fij => DenseMatrix(fij:_*) ),                      // for scala 2.11
                          fijab.map( fij => DenseMatrix(fij*) ),                        // for scala 2.13
                        //fijab.map( fij => new DenseMatrix(fij(0).size, fij.size, fij.flatten.toArray).t ),    // for scala 2.13
                                regTerm, propL1h, propL1J)    //, 0.0, 0.0, betaV, betaM, eps )
    }
  */

    def runBM(
                optionInitialInteractions: Option[MCMC.Interactions],   // If optionInitialInteractions == None/null, fia will be used to generate hia.
                initialConfigurations: ArraySeq[IArray[Byte]],      // For batch; If initialConfigurations(i).size == 0, randomConfiguration will be called.
                sampleConfigurations: ArraySeq[IArray[Byte]] = ArraySeq[IArray[Byte]](),       // Samples: the energy distribution of which is compared with its equilibrium distribution. 
                sampleIDs: ArraySeq[String] = ArraySeq[String](),
                initialConfigurationIndex: ArraySeq[Int] = ArraySeq[Int](),    //Vector[Int] = Vector[Int](),

                sigma_initial_J: Double = 0.001,   // Jijab is initialize by the Gaussian of this standard deviation and mean 0.0 for optionInitialInteractions == None/null
                                          // If this is equal to 0.0, then all elements will be initialized to be equal to 0.0.

                nInitialIterationsPerUnit: Int = 0, 
                everyNIterationsPerUnit: Int = 10,

                miniBatchSize: Int = 100,        // if miniBatchSize != 0, then initialConfigurations.size >= miniBatchSize

                // nSamples or initialNSamples and maxNSamples must be provided.
                // nSamples != 0 means initialNSamples=maxNSamples=nSamples.
                nSamples: Int = 1,                      // This should be non zero; 1 is appropriate.
                initialNSamples: Int = 0,               // obsolete
                incrementSamples: Double = 0.0,         // obsolete
                maxNSamples: Int = 0,                   // obsolete

                initialT: Double = 1.0,
                finalT: Double = 1.0,
                annealingRate: Double = 0.99,
                maxExtendedIterations: Int = 1,
                mcmcKernel: String = "GibbsWithMHStep", //"GibbsWithMHStep"=="MultiBlockMH", // or "MH" // or "Gibbs"

              //fia: ArraySeq[DenseVector[Double]],
              //fijab: ArraySeq[DenseMatrix[Double]],

              //regTerm: String,
              //propL1h: Double, propL1J: Double,
                lambdaPhi: Double, lambdaPhij: Double,
              ////bmState: BMState, 

                gauge: String = "unused",  // "phi_zeroSum" for L2 and L2GL1, "ungauged"=="unused" for all
                gradientDescentMethod: String = "ModAdamMax",
                betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-6,

                betaLAP: Double = 0.0,       // It is safe to set betaLAP = 0.0.

                learningRate: LearningRate = LearningRate(),
                learningRateForRPROPLR: LearningRateForRPROPLR = LearningRateForRPROPLR(),

              //minLearningRate: Double = 1.0e-6,
              //maxLearningRate: Double = 10.0,
              //rateDecrease: Double = 0.5,
              //rateIncrease: Double = 1.2,

                pseudoNCountsForBM: Double = 0.0,        // If the regularization is not used, this may be used to avoid 0 counts in this.fia and this.fijab
                pseudoNCountsForhia: Double = 10.0,      // used to generate hia from this.fia
                pseudoNCountsForProposedPia: Double = 10.0,    // used to generate proposedPia from this.fia 
                pseudoNCountsForKL: Double = 1.0,        // used to calculate KL valunes from ensemble averages of Pij(a,b)

                nBestKLs: Int = 1,
                minLearningStepsForBestKL: Int = 1710,
                minLearningSteps: Int = 1800,
                maxNoLearnings: Int = 90,

                nMiniBatchesForEnsembleAve: Int = 0,     // The ensemble average of TE will be done over max nMiniBatchesForEnsembleAve minibatch ensembles; 0 means over-epoch.
                betaLAPForKL: Double = 0.0,              // pseudo-counts for this leaky average are not clear so that simple average would be better.

                logInterval: Int = 100,
                logLevel: String = "default"     // Detail or Default
                ) = {

        require( pseudoNCountsForhia > 0.0 && pseudoNCountsForProposedPia > 0.0 )
        require( pseudoNCountsForBM >= 0.0 && pseudoNCountsForKL > 0.0 )

        require( pseudoNCountsForBM == 0.0, "*** Strongly recommended: pseudoNCountsForBM = 0.0, if regularization is employed.")
        require( betaLAPForKL == 0.0, "*** Strongly recommended: betaLAPForKL = 0")
        require( betaLAP == 0.0, "*** Strongly recommended: betaLAP = 0")

      //require(minLearningStepsForBestKL <= minLearningSteps )  // Now minLearningStepsForBestKL <= or >= minLearningSteps 

        require( mcmcKernel == "MH" || 
                        mcmcKernel == "GibbsWithMHStep" || mcmcKernel == "MultiBlockMH" || 
                        mcmcKernel == "Gibbs", 
                        "Not supported: mcmcKernel == MH, GibbsWithMHStep or Gibbs" )

        require( lambdaPhi >= 0.0 && lambdaPhij >= 0.0 )

        // If the Ising gauge is used, convergence will be not good.
        //require( gauge == "unused" )

        if ( this.regTerm == "L1L2" ) {
                require( gauge == "unused" || gauge == "ungauged" )
        } else {
                require( gauge == "unused" || gauge == "ungauged" || 
                        gauge == "phi_zeroSum" || gauge == "phi_ZeroSum" )      // gauge for phi and phij
        }

        val optMethod = 
            if ( gradientDescentMethod == "NAG" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "ModRAdamMax" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "ModRAdam" )
                gradientDescentMethod   // "ModRAdamMax"
            else if ( gradientDescentMethod == "ModAdamMax" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "ModAdamSum" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "ModAdam" )
                gradientDescentMethod   // "ModAdamMax"
            else if ( gradientDescentMethod == "Adam" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "RAdam" )
                gradientDescentMethod
            else if ( gradientDescentMethod == "Adadelta" )
                gradientDescentMethod
            else if  ( gradientDescentMethod == "RPROP-LR" || gradientDescentMethod == "MF" )
                "RPROP-LR"
            else {
                sys.error("Not supported: " + gradientDescentMethod) 
                ""
            }

        val nIndependentMC = initialConfigurations.size

        if ( miniBatchSize > 0 ) {
        //print(" %d  %d %d\n".format(nSamples, initialConfigurations.size, miniBatchSize) )
          require( initialConfigurations.size >= miniBatchSize) 
          require( nSamples >= 1 )
          require( sampleConfigurations.size == 0 || sampleConfigurations.size == initialConfigurations.size )
        } else {
          require( nSamples >= 0 )
          require( initialConfigurations.size > 0 ) 
        }

        val fullBatchSize = initialConfigurations.size
        val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)

        if ( nSamples == 0 ) {
          require( initialNSamples != 0 && maxNSamples != 0 )
        }

        val stepsPerEpoch = nMiniBatchesInFullBatch(miniBSize, fullBatchSize).toDouble

        val logInterval_rev = if (stepsPerEpoch <= 1.0) {
                                 if ( logInterval > 0 ) logInterval else 1
                              } else if ( stepsPerEpoch <= logInterval ) { 
                                 (stepsPerEpoch * (logInterval / stepsPerEpoch).floor ).toInt
                              } else {
                                 stepsPerEpoch.toInt
                              }

        val minLearningSteps_rev = (minLearningSteps.toDouble / logInterval_rev).ceil.toInt * logInterval_rev
        val minLearningStepsForBestKL_rev =  minLearningSteps_rev + (minLearningStepsForBestKL - minLearningSteps)
        val learningRate_rev = learningRate.newLearningRate( maxLREpochs = minLearningSteps_rev - learningRate.warmupEpochs )

        def printParameters(out: PrintStream): Unit = {
            out.print("# %s  propL1h= %g propL1J= %g  lambdaPhi= %g  lambdaPhij= %g  %s  %s  learningRate= %g  betaV= %g  betaM= %g  eps= %g  batchSize=%d  miniBatchSize= %d  stepsPerEpoch= %g".format(
                this.regTerm, this.propL1h, this.propL1J, lambdaPhi, lambdaPhij, gauge, optMethod, 
                learningRate_rev.maxLR, betaV, betaM, eps,
                initialConfigurations.size, miniBatchSize, stepsPerEpoch) )
            if ( betaLAP != 0.0 )
              out.print("  betaLAP= %g\n".format(betaLAP) )
            else
              out.print("\n")
            if ( optMethod == "RPROP-LR" ) {
                learningRateForRPROPLR.printParameters(out)
            } else { 
                learningRate_rev.printParameters(out)
            }

            if (optionInitialInteractions == null || optionInitialInteractions == None ) {
              out.print(s"# initial_h: to reproduce observed amino acid frequencies at each site; pseudoNCountsForhia= ${pseudoNCountsForhia}\n")
              out.print(s"# sigma_initial_J for Gaussian: ${sigma_initial_J}\n")
            } else {
              out.print(s"# hJ initialized: by optionInitialInteractions.get")
            }
            out.print(s"# MCMC Kernel: ${mcmcKernel}  with pseudoNCountsForProposedPia= ${pseudoNCountsForProposedPia}\n")
            if ( miniBatchSize > 0 ) {
                out.print("# nIndependentMC= %d".format(nIndependentMC))
                out.print("  miniBatchSize= %d  for MCMC \n".format(miniBatchSize))
                out.print("# initial configuration: observed sequence data\n")
            } else {
                out.print("# nIndependentMC= %d\n".format(nIndependentMC))
              //if ( initialConfigurations(0) == null || initialConfigurations(0).size == 0 ) {
                if ( initialConfigurations(0).size == 0 ) {
                    out.print("# initial configurations: random sequences generated with initial h\n")
                }
            }
            if ( nSamples != 0 )
                out.print("# nSamples for each MC run: %d\n".format(nSamples))
            else
                out.print("# initialNSamples = %d  incrementSamples = %g  maxNSamples = %d  for each MC run\n".format(
                                        initialNSamples, incrementSamples, maxNSamples))
            out.print(s"# nUnits= ${this.fia.size}  nInitialIterationsPerUnit= ${nInitialIterationsPerUnit}  everyNIterationsPerUnit= ${everyNIterationsPerUnit}\n")
            out.print("# initialT= %g  finalT= %g  annealingRate= %g ; annealingRate may be adjusted.\n".format(initialT, finalT, annealingRate))
            out.print("# maxExtendedIterations= %d\n".format(maxExtendedIterations))
            out.print("# minLearningStepsForBestKL= %d  minLearningSteps= %d  maxNoLearnings= %d\n".format(minLearningStepsForBestKL_rev, minLearningSteps_rev, maxNoLearnings) )
            out.print("#\n")
            out.print(s"# Frequencies, fia and fijab, are calculated from the effective number of sequences, effectiveNConfigs= ${this.effectiveNConfigs} ,\n")
            out.print(s"# with the corrections, pseudoNCountsForBM= ${pseudoNCountsForBM} , to avoid 0-counts.\n" )
            out.print("#\n")
        }

        val fileTE = new File(ioFiles("outputDir"), 
                //(s"Energy_distribution_${this.regTerm}_${this.propL1h}_${this.propL1J}_%.1e_%.1e_${gauge}_${optMethod}_${learningRate}_${betaV}_${betaM}.txt").format(lambdaPhi, lambdaPhij))
                (f"Energy_distribution_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate_rev.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
        val outTE  = new PrintStream(fileTE)

        printParameters(outTE)

        outTE.print("#First line:\n")
        outTE.print("# allRandomizedSampleTE_MeanAndVariance.mean - allRandomizedSampleTE_MeanAndVariance.variance: mean of TE in the Gaussian approx. for TE density\n")
        outTE.print("# sampleTE_MeanAndVariance.mean: mean of TE of native sequences for a mini-batch\n")
        outTE.print("# ensembleTE_MeanAndVariance.mean: <TE> : mean of TE over the ensemble generated in the MCMC for a mini-batch\n")
        outTE.print("#\n")
        outTE.print("# allRandomizedSampleTE_MeanAndVariance.mean: mean of TE of all randomized native sequences for a full-batch (allRandomizedSamples)\n")
        outTE.print("#\n")
        outTE.print("# allRandomizedSampleTE_MeanAndVariance.varience: variance of TE of all randomized native sequences for a full-batch (allRandomizedSamples)\n")
        outTE.print("# sampleTE_MeanAndVariance.varience: variance of TE of native sequences for a mini-batch\n")
        outTE.print("# ensembleTE_MeanAndVariance.varience: <(TE -<TE>)^2> : variance of TE over the ensemble generated in the MCMC for a mini-batch\n")
        outTE.print("#\n")
        outTE.print("#Second line:\n")
        outTE.print("# randomSampleTE_MeanAndVariance._2 - randomSampleTE_MeanAndVariance._3: expected mean of TE in the Gaussian approx. for TE density\n")
        outTE.print("# allSampleTE_MeanAndVariance.mean: mean of TE of all native sequences (sampleConfigurations)\n")
        outTE.print("# ensembleTE_MeanAndVariance.mean: <TE> : mean of TE over the ensemble generated in the MCMC for a mini-batch\n")
        outTE.print("# allEnsembleTE_MeanAndVariance.mean: <TE>_m : mean of TE over multiple minibatch ensembles generated in the MCMC\n" )
        outTE.print("#\n")
        outTE.print("# randomSampleTE_MeanAndVariance._2: mean of TE of random sequences with amino acid frequencies of native sequences (randomSamples)\n")
        outTE.print("#\n")
        outTE.print("# randomSampleTE_MeanAndVariance._3: expected variance of TE of random sequences with amino acid frequencies of native sequences (randomSamples)\n")
        outTE.print("# allSampleTE_MeanAndVariance.varience: variance of TE of all native sequences (sampleConfigurations)\n")
        outTE.print("# ensembleTE_MeanAndVariance.varience: <(TE -<TE>)^2> : variance of TE over the ensemble generated in the MCMC for a mini-batch\n")
        outTE.print("# allEnsembleTE_MeanAndVariance.varience: <(TE -<TE>)^2>_m: variance of TE over multiple minibatch ensembles generated in the MCMC for a full-batch\n")
        outTE.print("#\n")
        outTE.print("# The second line should be used;\n")
        outTE.print("#   The large variances for randomized sequences result from the variation of amino acid frequencies among native sequences.\n")
        outTE.print("#\n")
        
        val fileKL = new File(ioFiles("outputDir"), 
                (f"KL_of_each_step_${this.regTerm}_${this.propL1h}_${this.propL1J}_${lambdaPhi}%.1e_${lambdaPhij}%.1e_${gauge}_${optMethod}_${learningRate_rev.maxLR}_${betaV}_${betaM}_${initialConfigurations.size}_${miniBatchSize}_${sigma_initial_J}.txt"))
        val outKL  = new PrintStream(fileKL)

        printParameters(outKL)
      //printKL(outKL, bmStep = -1, nNonEquil = 0, nExtendedIterations = 0, (0.0, 0.0), (0.0, 0.0) )     // print a title

      //(valueKL1, valueKL2, step, fileInteractions, fileIndpendentMCsamplings)
        @annotation.tailrec
        def fillNone(i: Int, queue: Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ]):
                Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ] = {
                assert( i > 0 )
                val none: Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] = None
                if(i == 1)
                  none +: queue
                else
                  fillNone(i - 1, none +: queue) 
        }
        val bestKLs = fillNone(nBestKLs, Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ]() )
        val bestKL = fillNone(1, Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ]() )


        def preFirstBatch(miniBatchSize: Int,
                                  initialConfigurationIndex: ArraySeq[Int],
                                  initialConfigurations: ArraySeq[MCMC.Configuration],
                                  sampleConfigurations: ArraySeq[MCMC.Configuration],
                                  sampleIDs: ArraySeq[String],
                                  bmState: BMState ) = { 

            require( initialConfigurations.size > 0 )
            val fullBatchSize = initialConfigurations.size

            val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize )           
            val batchSize = miniBSize

            if ( sampleConfigurations.size != 0 ) {
              require( sampleConfigurations.size == initialConfigurations.size )
            }

            val (revisedSampleConfigurations, revisedSampleIDs) = {
                  val revSampleIDs = if ( sampleIDs.size == 0 ) {
                        ArraySeq.fill(initialConfigurations.size)("")
                  } else {
                        require( sampleIDs.size == initialConfigurations.size)
                        sampleIDs
                  }
                  if ( initialConfigurationIndex.size == 0 ) {
                       if ( sampleConfigurations.size == 0 ) {
                              (initialConfigurations, revSampleIDs)
                       } else {
                              (sampleConfigurations, revSampleIDs)
                       }
                  } else {
                       require( initialConfigurationIndex.size == initialConfigurations.size )

                       if ( sampleConfigurations.size == 0 ) {
                           val sortByIndex = initialConfigurationIndex.zipWithIndex
                           val sortedByIndex = stableSort( sortByIndex, (x, y) => x._1 <= y._1 )
                           val inverseIndex = ArraySeq.range(0, initialConfigurationIndex.size).map(i => sortedByIndex(i)._2)
                           ( inverseIndex.map(i => initialConfigurations(i)), inverseIndex.map(i => revSampleIDs(i)) )
                       } else {
                           (sampleConfigurations, revSampleIDs)
                       }
                  }
            }
            val revInitialConfigurationIndex = if ( initialConfigurationIndex.size == 0 ) {
                Vector.range(0, initialConfigurations.size)
            } else {
                initialConfigurationIndex
            }
 
            val nativeConfigs: ArraySeq[MCMC.NativeConfig] = ArraySeq.range(0, revisedSampleConfigurations.size).map{ i =>
                val id = revisedSampleIDs(i)
                val sample = revisedSampleConfigurations(i)
                MCMC.NativeConfig(index = i, id = id, sample = sample)
            }

          /*
            val stateOrderString = "ARNDCQEGHILKMFPSTWYV-"
            cfor(0)(i => i < initialConfigurations.size, i => i + initialConfigurations.size - 1 ) { i =>
              val confIndex = initialConfigurationIndex(i)
              stderr.print("%d  %d  %s  ".format(i, nativeConfigs(confIndex).index, nativeConfigs(confIndex).id) )
              cfor(0)(j => j < nativeConfigs(confIndex).sample.size, j => j + 1 ) { j =>
                stderr.print("%c".format( stateOrderString(nativeConfigs(confIndex).sample(j) ) ) )
              }
              stderr.print("\n")
            }
            sys.error("Debug ends.")
          */
          
          // Shuffling of configuration order
          // If previous ensembles are used, random shuffling must be done even for miniBSize == fullBatchSize

            val shuffledInitialConfigs = {
                val randomOrder = Random.shuffle( Vector.range(0, initialConfigurations.size) )
                randomOrder.map( i => (revInitialConfigurationIndex(i), initialConfigurations(i)) )
            }

            val initialEachInFullBatch: ParVector[MCMC.State] = shuffledInitialConfigs.par.map{ c =>
                  val index  = c._1
                  val config = c._2
                  MCMC.State( nativeConfig = nativeConfigs(index) , configuration = config )
            }

          /*
            val stateOrderString = "ARNDCQEGHILKMFPSTWYV-"
            cfor(0)(i => i < initialEachInFullBatch.size, i => i + initialEachInFullBatch.size - 1 ) { i =>
              val nativeConfig = initialEachInFullBatch(i).nativeConfig
              val config = initialEachInFullBatch(i).configuration
              stderr.print("%d  %d  %s  ".format(i, nativeConfig.index, nativeConfig.id) )
              cfor(0)(j => j < config.size, j => j + 1 ) { j =>
                stderr.print("%c".format( stateOrderString(config(j) ) ) )
              }
              stderr.print("\n")
            }
            sys.error("Debug ends.")
          */

            val nMBinFull = nMiniBatchesInFullBatch(miniBSize, fullBatchSize)
 
            val nextLastAndCurrentFullBatch: Vector[BMState] =
                Vector.range(0, nMBinFull).map{ iMB => 
                    val startPosInFull  = startPosInFB(iMB, miniBSize, fullBatchSize )
                    val (s, e, untilE) = endPosOfMiniBatch(startPosInFull, miniBSize, fullBatchSize)
                    val initialMCStates: MCMC.InitialMCStates = initialEachInFullBatch.slice(s, untilE)
                    if ( iMB == 0 ) {
                        bmState.copy( initialMCStates = initialMCStates,
                          independentSamplings = ParVector[ArraySeq[MCMC.State]](), ensembleAverages = MCMC.EnsembleAverages() )
                    } else {
                        bmState.copy( step = scala.Int.MinValue, initialMCStates = initialMCStates,
                          independentSamplings = ParVector[ArraySeq[MCMC.State]](), ensembleAverages = MCMC.EnsembleAverages() )
                    }
                }

            val nextBMState = nextLastAndCurrentFullBatch(0)

            val nextPrevMiniBatches: Vector[BMState] = Vector[BMState]()

            val nextStartPos = 0
 
            (
             nativeConfigs: ArraySeq[MCMC.NativeConfig],
             nextLastAndCurrentFullBatch: Vector[BMState],
             nextPrevMiniBatches: Vector[BMState],
             nextStartPos: Int,
             nextBMState: BMState
            )
        }

        def postMiniBatch(miniBatchSize: Int,
            lastAndCurrentFullBatch: Vector[BMState],
            prevMiniBatches: Vector[BMState],
            startPos: Int,
            bmState: BMState,
            nMiniBatchesForEnsembleAve: Int,
            nextBMState: BMState ) = {

          val fullBatchSize = lastAndCurrentFullBatch.map{ bms => math.max(bms.initialMCStates.size, bms.independentSamplings.size) }.sum

          val miniBSize = miniBatchSizeCorrected(miniBatchSize, fullBatchSize)

          assert( startPos % miniBSize == 0 )

          val nMBinFull = nMiniBatchesInFullBatch(miniBSize, fullBatchSize)

          val (_, endPos, untilEndPos) = endPosOfMiniBatch( startPos, miniBSize, fullBatchSize )

          assert( startPos + bmState.initialMCStates.size == untilEndPos )

        //val miniBatchPos = startPos / miniBSize
          val miniBatchPos = startMBPosInFB(startPos, miniBSize, fullBatchSize)

          val cleanBMS = cleanBMState(bmState)
          val revisedLastAndCurrentFullBatch = {
              // Release unnecessary memories except memories for future uses.
              copyToSlice( (miniBatchPos, miniBatchPos + 1), lastAndCurrentFullBatch, Vector( cleanBMS ) )
          }

          assert( revisedLastAndCurrentFullBatch(miniBatchPos).step >= 0,
                 "revicedLastAndCurrentFullBatch(miniBatchPos).step: %d\n".format(revisedLastAndCurrentFullBatch(miniBatchPos).step) )

          val (nextStartPos, nextEndPos, nextUntilEndPos ) = endPosOfMiniBatch(startPos = endPos, miniBSize, fullBatchSize)

          assert( nMBinFull == revisedLastAndCurrentFullBatch.size ) 

          val (nextLastAndCurrentFullBatch, nextPrevMiniBatches ) = {
              val nextPrev = 
                         if ( nextStartPos > 0 ) {
                              val tooManyMBs = prevMiniBatches.size + (miniBatchPos + 2) - nMiniBatchesForEnsembleAve 
                              if ( tooManyMBs <= 0 ) {
                                       prevMiniBatches
                              } else {
                                       prevMiniBatches.slice(tooManyMBs, prevMiniBatches.size)
                              }
                         } else if (nextStartPos == 0 ) {
                              val total = prevMiniBatches ++ revisedLastAndCurrentFullBatch
                              val s = math.max(total.size - (nMiniBatchesForEnsembleAve - 1), 0)
                              if ( s == 0 )
                                 total
                              else
                                 total.slice(s, total.size)
                         } else {
                            prevMiniBatches
                         } 

            //if ( nextStartPos == 0 && nextEndPos != 0 ) { // If previous ensembles are used, random shuffling must be done even for miniBSize == fullBatchSize
              if ( nextStartPos == 0 ) { 
                  // Pick the last configuration for each MC sampling
                  val eachInFull: Vector[MCMC.State] = revisedLastAndCurrentFullBatch.map{ bms => 
                                                           assert(bms.step >= 0)
                                                           bms.independentSamplings.map{ states => states(states.size - 1) } 
                                                    }.flatten

                  assert( eachInFull.size == fullBatchSize)
      
                  val nextEachInFull = {
                    //val randomOrder = Random.shuffle( Vector.range(0, fullBatchSize) )
                    //randomOrder.par.map( i => eachInFull(i) )
                      Random.shuffle( eachInFull )
                  }

                  val nextCurrentFullBatch: Vector[BMState] = Vector.range(0, nMBinFull).map{ iMB => 
                        val startPosInFull  = startPosInFB(iMB, miniBSize, fullBatchSize )
                        val (s, e, untilE) = endPosOfMiniBatch(startPos = startPosInFull, miniBSize, fullBatchSize)
                        val initialMCStates = nextEachInFull.slice(s, untilE).map{ s =>
                            s.copy( interactions = MCMC.Interactions(), energy = 0.0, kT = 0.0, annealingRate = 0.0, step = 0, finalkT = 0.0 ) 
                        }
                        cleanBMS.copy( step = scala.Int.MinValue, initialMCStates = initialMCStates.par,
                           independentSamplings = ParVector[ArraySeq[MCMC.State]](), ensembleAverages = MCMC.EnsembleAverages() )
                  }

                  ( nextCurrentFullBatch, nextPrev )
              } else {
                  ( revisedLastAndCurrentFullBatch, nextPrev )
              } 
          }

          val revisedNextBMState = {
                     val nextStartMBPos = startMBPosInFB(nextStartPos, miniBSize, fullBatchSize)
                     val bms = nextLastAndCurrentFullBatch( nextStartMBPos )
                     assert( bms.step < 0 )             //It confirms that this bmState is initialized.
                     nextBMState.copy( initialMCStates = bms.initialMCStates,
                                       independentSamplings = bms.independentSamplings,
                                       ensembleAverages = bms.ensembleAverages )
          }
          assert( revisedNextBMState.step >= 0 )
         
          (nextLastAndCurrentFullBatch, nextPrevMiniBatches, revisedLastAndCurrentFullBatch, nextStartPos, revisedNextBMState)
        }

      /*
        require( pseudoNCountsForhia > 0 && pseudoNCountsForProposedPia > 0 )
        require( pseudoNCountsForBM >= 0 && pseudoNCountsForKL > 0 )

        val (fiaForBM, fijabForBM, faForBM ) = 
            if ( pseudoNCountsForBM > 0.0 ) {
              val (fiaForBM, fijabForBM) = bayesianCorrection(this.fia, this.fijab, this.effectiveNConfigs, pseudoNCountsForBM)
              val faForBM = bayesianCorrection( this.fa, this.effectiveNConfigs * this.fia.size, pseudoNCountsForBM) 
              (fiaForBM, fijabForBM, faForBM )
            } else {
              (this.fia, this.fijab, this.fa )
            }

        val (fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1) = {
            if ( pseudoNCountsForBM > 0.0 ) {
              (fiaForBM, fijabForBM, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]]() )
            } else {
            // Be carefull
            /*
              if ( pseudoNCountsForKL > 0.0 ) {
                val (fiaForKL, fijabForKL) = bayesianCorrection(this.fia, this.fijab, this.effectiveNConfigs, pseudoNCountsForKL)
                (fiaForKL, fijabForKL, ArraySeq[DenseVector[Double]](), ArraySeq[DenseMatrix[Double]]() )
              } else
            */
              {
                val (fiaWith0ReplacedBy1, fijabWith0ReplacedBy1) = replace0By1ForKL(this.fia, this.fijab)
                (this.fia, this.fijab, fiaWith0ReplacedBy1, fijabWith0ReplacedBy1 )
              }
            }
        }
      */

        val (fiaForBM, fijabForBM, faForBM,  fiaForKL, fijabForKL, fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1) =
            freqsForBMandKL( this.fia, this.fijab, this.fa,
                this.effectiveNConfigs,
                pseudoNCountsForBM,        // If the regularization is not used, this may be used to avoid 0 counts in this.fia and this.fijab
                pseudoNCountsForKL         // used to calculate KL valunes from ensemble averages of Pij(a,b)
            )

        val calcKL = BM.calcKL(fiaForKL, fijabForKL, observedNP = this.effectiveNConfigs, pseudoNP = 0.0, 
                               _: ArraySeq[DenseVector[Double]], _: ArraySeq[DenseMatrix[Double]], _: Double, _: Double,
                               fiaForKLWith0ReplacedBy1, fijabForKLWith0ReplacedBy1 )

        val proposedPia = bayesianCorrection( this.fia, this.effectiveNConfigs, pseudoNCountsForProposedPia ) 

        val runBM1Step = BM.runBM1Step(
                this.ioFiles,
                this.stateOrderString,
                _: MCMC.Interactions,
             // _: ArraySeq[State],         // _: ArraySeq[IArray[Byte]],
                nInitialIterationsPerUnit,
                everyNIterationsPerUnit,
                _: Int,
              //nIndependentMC,
                initialT,
                finalT,
                annealingRate,
                maxExtendedIterations,
                mcmcKernel,
                proposedPia,

              // Please specify ...
              //_: ArraySeq[DenseVector[Double]],
              //_: ArraySeq[DenseMatrix[Double]],

              //this.fia, this.fijab, // replaced by fiaForBM, fijabForBM
                fiaForBM, fijabForBM, // In the case of no regularization fiaForBM and fijabForBM can be employed to avoid the zero value in fia and fijab.

                this.regTerm,
                this.propL1h, this.propL1J,
                _: Double, _: Double, //this.lambdaPhi, this.lambdaPhij,
                _: BMState,
                _: String             //optMethod
                )

        @annotation.tailrec
        def iterateBMSteps(
              //runBM1Step: (MCMC.Interactions, ArraySeq[IArray[Byte]], Int, Double, Double, BMState, String ) =>
              //runBM1Step: (MCMC.Interactions, ArraySeq[MCMC.State], Int, Double, Double, BMState, String ) =>
                runBM1Step: (MCMC.Interactions, Int, Double, Double, BMState, String ) =>
                            (MCMC, Int, Int, BMState, BMState) ,
                interactions: MCMC.Interactions, 

              //fia, fijab, fa,                               // replaced by f.*withPseudo to avoid too-large interactions.
                fiaForBM: ArraySeq[DenseVector[Double]],
                fijabForBM: ArraySeq[DenseMatrix[Double]],
                faForBM: DenseVector[Double],

                miniBatchSize: Int,
              //initialConfigurations: ArraySeq[IArray[Byte]],   // for minibatch
              //initialMCStates: ArraySeq[MCMC.State],             // for minibatch

                lastAndCurrentFullBatch:        Vector[BMState],
                prevMiniBatches:                Vector[BMState],

                startPos: Int,
                nSamples: Int,
              //initialNSamples: Int,          // obsolete
              //incrementSamples: Double,      // obsolete
              //maxNSamples: Int,              // obsolete

                lambdaPhi: Double = lambdaPhi, lambdaPhij: Double = lambdaPhij,

                bmState: BMState,

                bestKLs: Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ],
                bestKL:  Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ],

                minLearningStepsForBestKL: Int,
                minLearningSteps: Int,
                maxNoLearnings: Int,

                nMiniBatchesForEnsembleAve: Int,

                logInterval: Int,

                optimizer: String = optMethod
                ): (MCMC.Interactions, Vector[BM.BMState]) = {

        //val initialConfigIndex = lastConfigIndex

        //val nSamples = { val n = initialNSamples + (incrementSamples * bmState.step).floor.toInt
        //                   if ( n > maxNSamples ) maxNSamples else n
        //               }
        //val initialNSamples = nSamples
        //val incrementSamples = 0.0
        //val maxNSamples = nSamples

          val fullBatchSize = lastAndCurrentFullBatch.map{ bms => math.max(bms.initialMCStates.size, bms.independentSamplings.size) }.sum

          val miniBSize = miniBatchSizeCorrected (miniBatchSize, fullBatchSize)  
          val nMBsForEnsembleAve = nMiniBatchesForEnsembleAverage(miniBSize, math.max(fullBatchSize, 50000), nMiniBatchesForEnsembleAve )

          val (mcmc, nNonEquil, nExtendedIterations, bmStateWithMCSamples, nextBMState_0) = 
                runBM1Step( interactions, nSamples, lambdaPhi, lambdaPhij, bmState, optimizer )
              //runBM1Step( interactions, initialConfigurations, nSamples, lambdaPhi, lambdaPhij, bmState, optimizer )

          val independentSamplings = bmStateWithMCSamples.independentSamplings
          val lastMiniBatchSamplings = independentSamplings
        //val revBMState = bmStateWithMCSamples

          val ensembleAverages = bmStateWithMCSamples.ensembleAverages

          val (nextLastAndCurrentFullBatch, nextPrevMiniBatches, revisedLastAndCurrentFullBatch, nextStartPos, nextBMState)
                   = postMiniBatch(miniBSize,                                                 // initialMCStatesOfBMStep,
                       lastAndCurrentFullBatch, prevMiniBatches, startPos, bmStateWithMCSamples, nMBsForEnsembleAve, nextBMState_0)

          val (_, endPos, untilEndPos) = endPosOfMiniBatch( startPos, miniBSize, fullBatchSize)

          val nMBsInFullB = nMiniBatchesInFullBatch(miniBSize, fullBatchSize ) 
          assert( revisedLastAndCurrentFullBatch.size == nMBsInFullB )
          val stepsPerAveBlock = nMBsForEnsembleAve

        //val (fromPosForEnsembleAve, untilPosForEnsembleAve, _ ) = BM.rangeExclForEnsembleAverage(
        //       prevSamplings_of_IndependentMC.size + revisedLastSamplings_of_IndependentMC.size, 
        //       revisedLastSamplings_of_IndependentMC.size,
        //       miniBatchSize, (startPos, untilEndPos), nMBsForEnsembleAve )
 
          val startMBPosInFull = startMBPosInFB(startPos, miniBSize, fullBatchSize)
          val (startMBPos, untilEndMBPos) = miniBatchRangeExclForEnsembleAverage(prevMiniBatches.size , revisedLastAndCurrentFullBatch.size,
                                  startMBPosInFull = startMBPosInFull, nMiniBatchesForEnsembleAveCorrected = nMBsForEnsembleAve )

          val stepsPerEpoch = nMBsInFullB.toDouble  //(fullBatchsize / miniBSize.toDouble ).floor 

          val totalNSamples = (independentSamplings.size * independentSamplings(0).size).toDouble 
      
        //val ensembleKL = BM.calcKL(this.fia, this.fijab, observedNP = this.effectiveNConfigs, pseudoNP = pseudoNCounts, 
        //                        ensembleAverages.pia, ensembleAverages.pijab, observedNQ = totalNSamples, pseudoNQ = pseudoNCounts)
        //val ensembleKL = BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0, 
        //                        ensembleAverages.pia, ensembleAverages.pijab, observedNQ = totalNSamples, pseudoNQ = pseudoNCountsForKL)
          val ensembleKL = calcKL( ensembleAverages.pia, ensembleAverages.pijab, totalNSamples, pseudoNCountsForKL )
          if(bmStateWithMCSamples.step == 0) {
                val pia = fiaForBM    // mcmc.proposedDistributions.proposedDistributionsAtAllUnits
                val nUnits = pia.size
                val pijab = ArraySeq.range(0, (nUnits * (nUnits - 1)) / 2).map { ij =>
                        val (i, j) = inversePairIndex(ij)
                        val (ij_ , a, b) = pairIndex(i, j, 0, 1)
                        assert( ij == ij_ )
                      //pijab(ij) =
                        if ( a == 0 ) 
                                pia(i) * pia(j).t
                        else {
                                assert(a == 0, "Please check pairIndex.\n")
                                pia(j) * pia(i).t
                        }
                }
              //val aveKL = BM.calcKL(this.fia, this.fijab, observedNP = this.effectiveNConfigs, pseudoNP = pseudoNCounts, 
              //val aveKL = BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0, 
              //                   pia, pijab, observedNQ = Double.PositiveInfinity, pseudoNQ = 0.0)
                val aveKL = calcKL( pia, pijab, this.effectiveNConfigs, 0.0 )

                val aveKL2 = if(bmStateWithMCSamples.betaLAPForKL == 0.0 ) { 
                                 (0.0, 0.0)
                             } else {
                               //BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0, 
                               //BM.calcKL(fiaForBM, fijabForBM, observedNP = this.effectiveNConfigs, pseudoNP = 0.0, 
                               //  bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                               //  observedNQ = totalNSamples, pseudoNQ = pseudoNCountsForKL )
                                 calcKL( bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                                         totalNSamples, pseudoNCountsForKL )
                             }

                outKL.print("#")
                printKL(outKL, 0, nNonEquil = 0, nExtendedIterations = 0, ensembleKL, aveKL, aveKL2 )
          }

          val miniBacthesForEnsembleAve = 
                    (prevMiniBatches ++ revisedLastAndCurrentFullBatch).slice(startMBPos, untilEndMBPos)

          val (aveKL, aveKL2) = {
            if ( startMBPos == prevMiniBatches.size + startMBPosInFull ) {
              assert( untilEndMBPos == startMBPos + 1 )
              ( ensembleKL , 
                if(bmStateWithMCSamples.betaLAPForKL == 0.0 ) {
                    (0.0, 0.0)
                } else {
                  //BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0,
                  //BM.calcKL(fiaForBM, fijabForBM, observedNP = this.effectiveNConfigs, pseudoNP = 0.0,
                  //   bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                  //   observedNQ = totalNSamples, pseudoNQ = pseudoNCountsForKL )
                    calcKL( bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                            totalNSamples, pseudoNCountsForKL )
                }
              )

            } else {
            //val rangeExclMBForEnsembleAve = 
            //      (prevMiniBatches ++ revisedLastAndCurrentFullBatch).slice(startMBPos, untilEndMBPos)
              val independentS: ParVector[ArraySeq[MCMC.State]] =
                  miniBacthesForEnsembleAve.map{bms => bms.independentSamplings}.flatten.par

            //val nStatesOfUnit = this.fia(0).size
              val nStatesOfUnit = fiaForBM(0).size
              val pia = MCMC.frequenciesAtUnitInSamples( independentS, nStatesOfUnit )
              val pijab = MCMC.pairwiseFrequenciesInSamples( independentS, nStatesOfUnit )
            //calcKL(this.fia, this.fijab, effectiveNConfigs, pseudoNCounts, pia, pijab, independentS.size.toDouble, pseudoNCounts)
              (
              //BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0, 
              //BM.calcKL(fiaForBM, fijabForBM, observedNP = this.effectiveNConfigs, pseudoNP = 0.0, 
              //           pia, pijab, observedNQ = (independentS.size * independentS(0).size).toDouble, pseudoNQ = pseudoNCountsForKL) ,
                calcKL( pia, pijab, (independentS.size * independentS(0).size).toDouble, pseudoNCountsForKL) ,
                if(bmStateWithMCSamples.betaLAPForKL == 0.0 ) {
                    (0.0, 0.0)
                } else {
                  //BM.calcKL(fiaForBM, fijabForBM, observedNP = Double.PositiveInfinity, pseudoNP = 0.0, 
                  //BM.calcKL(fiaForBM, fijabForBM, observedNP = this.effectiveNConfigs, pseudoNP = 0.0, 
                  //   bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                  //   observedNQ = (independentS.size * independentS(0).size).toDouble, pseudoNQ = pseudoNCountsForKL )
                    calcKL( bmStateWithMCSamples.leakyAvePForKL.pia, bmStateWithMCSamples.leakyAvePForKL.pijab,
                            (independentS.size * independentS(0).size).toDouble, pseudoNCountsForKL )
                }
              )
            }
          }
          printKL(outKL, bmStateWithMCSamples.step, nNonEquil, nExtendedIterations, ensembleKL, aveKL, aveKL2 )  // for the step-th interactions

          val bestKLsRevised =
            if ( (bmStateWithMCSamples.step + 1) % logInterval == 0 || ( (bmStateWithMCSamples.step + 1) % logInterval ) % stepsPerAveBlock == 0 || bmStateWithMCSamples.step == 0) {

                val logFiles: Option[Tuple2[Option[java.io.File], Option[java.io.File]]] =
                  if ( (bmStateWithMCSamples.step + 1) % logInterval == 0 || bmStateWithMCSamples.step == 0 ) {
                    val logFileInteractions = new File(ioFiles("outputDir"), s"log_Interactions_hJ_${bmStateWithMCSamples.step}.txt.gz")
                    val logFileIndpendentMCsamplings = new File(ioFiles("outputDir"), s"log_MC_samples_${bmStateWithMCSamples.step}.fasta.gz")
                  //val logFileIndpendentEnergies = new File(ioFiles("outputDir"), s"log_Config_energy_distributions_${bmStateWithMCSamples.step}.txt")
                    Option( (Option(logFileInteractions), Option(logFileIndpendentMCsamplings)) )
                  } else {
                    None
                  }

                val (interactions_Ising, vecIndependentSamplings_Ising_ForEnsembleAve,
                     ensembleTE_MeanAndVariance, sampleTE_MeanAndVariance, allRandomizedSampleTE_MeanAndVariance,
                     allEnsembleTE_MeanAndVariance, allSampleTE_MeanAndVariance, randomSampleTE_MeanAndVariance) =
                          BM.printLog( bmStateWithMCSamples.step, logFiles, this.stateOrderString,
                          interactions, miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch,
                          outTE, faForBM, startPos, miniBSize, nMiniBatchesForEnsembleAve = nMBsForEnsembleAve )

              //val nMinusR_varOfTE = allSampleTE_MeanAndVariance.variance - randomSampleTE_MeanAndVariance._3
              //val sT = randomSampleTE_MeanAndVariance._3 / (randomSampleTE_MeanAndVariance._2 - allSampleTE_MeanAndVariance.mean)

              //if ( bestKLs.size < (logInterval / 2.0) ) deleteBestKLs(bestKLs) else bestKLs
                if ( bestKLs.size < (logInterval / 2.0) )
                        fillNone(bestKLs.size, Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ]() )
                else
                        bestKLs
            } else {
                bestKLs
            }

          val fileInteractions = new File(ioFiles("outputDir"), s"Interactions_hJ_${bmStateWithMCSamples.step}.txt.gz")
          val fileIndpendentMCsamplings = new File(ioFiles("outputDir"), s"MC_samples_${bmStateWithMCSamples.step}.fasta.gz")

          val (newBestKLs, oneOfBestKLs) = keepBestKLs(bestKLsRevised, Tuple5(aveKL._1, aveKL._2, bmStateWithMCSamples.step, Option(fileInteractions), Option(fileIndpendentMCsamplings) ) )
          if( oneOfBestKLs && logLevel != "default" ) {
              val out = if ( (bmStateWithMCSamples.step + 1) % logInterval == 0 ) {
                        null
                      } else {
                        outTE
                      }
              val (interactions_Ising, vectIndependentSamplings_Ising_ForEnsembleAve, ensembleTE_MeanAndVariance, sampleTE_MeanAndVariance,
                    allRandomizedSampleTE_MeanAndVariance, allEnsembleTE_MeanAndVariance, allSampleTE_MeanAndVariance, randomSampleTE_MeanAndVariance) =
                          BM.printLog( bmStateWithMCSamples.step, Option( (Option(fileInteractions), Option(fileIndpendentMCsamplings)) ),
                          this.stateOrderString, interactions, miniBacthesForEnsembleAve, revisedLastAndCurrentFullBatch,
                          out, faForBM, startPos, miniBSize, nMiniBatchesForEnsembleAve = nMBsForEnsembleAve)
          }

 //Error: val (newBestKL, oneOfBestKL) = if ( nSamples >= maxNSamples && bmStateWithMCSamples.step - bmStateWithMCSamples.stepOffset >= minLearningStepsForBestKL ) {
          val (newBestKL, oneOfBestKL) = if ( nSamples >= maxNSamples && bmStateWithMCSamples.step + 1 >= minLearningStepsForBestKL ) {
                                                val none: Option[File] = None
                                                keepBestKLs(bestKL, Tuple5(aveKL._1, aveKL._2, bmStateWithMCSamples.step, none, none)) 
                                         } else {
                                                (bestKL, false)
                                         }
          if ( oneOfBestKL ) {
              outKL.print( s"# bestKL= ${newBestKL(0).get._1} ${newBestKL(0).get._2}  step= ${newBestKL(0).get._3}\n" )
          }

        //val (minLearningSteps_new, newOptimizer, param) = if(
          val (newStepOffset, minLearningSteps_rev, minLearningStepsForBestKL_rev, newOptimizer, param ) = if( 
              //  (bmStateWithMCSamples.step + 1) % logInterval == 0 &&
         //Error: bmStateWithMCSamples.step - bmStateWithMCSamples.stepOffset + 1 >= minLearningSteps ) {
         //Error: bmStateWithMCSamples.step - bmStateWithMCSamples.stepOffset + 1 == 
                  bmStateWithMCSamples.step + 1 == 
                      ( (minLearningSteps.toDouble / logInterval).ceil.toInt * logInterval )  ) {

              if ( ioFiles.contains("jobControlParamFile") && ioFiles("jobControlParamFile").exists() ) {
                  val jobControlParamFile = ioFiles("jobControlParamFile")
                
                  val param = MCMC.paramFileToHashMap(jobControlParamFile)

                  val newOpt = if( param.contains("gradientDescentMethod") ) {     
                                     param("gradientDescentMethod")
                               } else {
                                     optimizer
                               }

                  val minSteps = if( param.contains("minLearningSteps") ) {
                                     val min_rev = param("minLearningSteps").toInt

                                     val minLearning = (min_rev.toDouble / logInterval).ceil.toInt * logInterval
                                     require( minLearning > 0 )
                                     minLearning
                                 } else {
                                     minLearningSteps - nextBMState.stepOffset
                                 }

                  if ( newOpt != optimizer ) {
                      val newOffset = nextBMState.step
                    //val offset = nextBMState.stepOffset
                    //val minStepsForBestKL = (minSteps * ( (minLearningStepsForBestKL - offset).toDouble / (minLearningSteps - offset) ) ).floor.toInt
                    //assert( minStepsForBestKL <=  minSteps && minStepsForBestKL > 0 )
                      val minStepsForBestKL = minSteps + (minLearningStepsForBestKL - minLearningSteps)  
                      assert( minStepsForBestKL > 0 )

                      (newOffset, newOffset + minSteps, newOffset + minStepsForBestKL, newOpt, param )
                //} else if ( minSteps != minLearningSteps - nextBMState.stepOffset ) {
                  } else if ( minSteps > minLearningSteps - nextBMState.stepOffset ) {
                      val offset = nextBMState.stepOffset
                    //val minStepsForBestKL = (minSteps * ( (minLearningStepsForBestKL - offset).toDouble / (minLearningSteps - offset) )  ).floor.toInt
                    //assert( minStepsForBestKL <=  minSteps && minStepsForBestKL > 0 )
                      val minStepsForBestKL = minSteps + (minLearningStepsForBestKL - minLearningSteps)  
                      assert( minStepsForBestKL > 0 )

                    //error: 250413 (nextBMState.stepOffset, minLearningSteps, offset + minStepsForBestKL, optimizer, param )
                      (nextBMState.stepOffset, offset + minSteps, offset + minStepsForBestKL, optimizer, param )
                  } else {
                      (nextBMState.stepOffset, minLearningSteps, minLearningStepsForBestKL, optimizer, param )
                  }

              } else {
                /* 
                  if ( optimizer == "ModAdamMax" ) {
                  //(minLearningSteps * 2, "Adam", HashMap[String, String]()) 
                    (minLearningSteps, "Adam", HashMap[String, String]()) 
                  } else {
                    (minLearningSteps, optimizer, HashMap[String, String]())
                  }
                */
                (nextBMState.stepOffset, minLearningSteps, minLearningStepsForBestKL, optimizer, HashMap[String, String]())
              }

          } else {
                (nextBMState.stepOffset, minLearningSteps, minLearningStepsForBestKL, optimizer, HashMap[String, String]())
          }

          if( newOptimizer == optimizer &&
               (bmStateWithMCSamples.step + 1) % logInterval == 0 && 
      //Error: bmStateWithMCSamples.step - bmStateWithMCSamples.stepOffset + 1 >= minLearningSteps_rev &&
               bmStateWithMCSamples.step + 1 >= math.max( minLearningStepsForBestKL_rev, minLearningSteps_rev ) &&
               bmStateWithMCSamples.step - newBestKL(0).get._3 >= maxNoLearnings ) {
              
              ( BM.hJ_in_IsingGauge(interactions), revisedLastAndCurrentFullBatch)

          } else {

              val nextBMState_rev = {

                  val newBetaLAP = if( param.contains("betaLAP") ) {
                                     param("betaLAP").toDouble
                               } else {
                                     nextBMState.betaLAP
                               }

                /* this block was moved forward.
                  val newStepOffset = if ( newOptimizer != optimizer ) {
                     nextBMState.step
                  } else {
                     nextBMState.stepOffset
                  }
                */

                /* 241011 Error
                  val newMaxLREpochs = if( param.contains("maxLREpochs") ) {
                      param("maxLREpochs").toDouble
                //} else if ( nextBMState.learningRate.maxLREpochs < minLearningSteps_rev - newStepOffset ) {
                //    minLearningSteps_rev.toDouble - newStepOffset
                  } else {
                     nextBMState.learningRate.maxLREpochs
                  }
                */
                //val newMaxLREpochs = minLearningSteps_rev - minLearningSteps - nextBMState.learningRate.warmupEpochs 
                  val newMaxLREpochs = minLearningSteps_rev - newStepOffset - nextBMState.learningRate.warmupEpochs 
                  assert( newMaxLREpochs > 0.0 )

                  val newMaxLR = if ( param.contains("maxLR") ) {
                     param("maxLR").toDouble
                  } else {
                     nextBMState.learningRate.maxLR
                  }

                  if ( newOptimizer != optimizer ) {

                      // Default values are employed for betaV, and betaM.
                      val initialState: BMState = {
                        //initializeBMState(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate, nextBMState.learningRateForRPROPLR)
                          if (newOptimizer == "NAG" ) {
                            initializeBMStateForNAG(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "ModRAdamMax" || newOptimizer == "ModRAdam" ) {
                            initializeBMStateForModRAdamMax("ModRAdamMax", nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "ModAdamMax" || newOptimizer == "ModAdam" ) {
                            initializeBMStateForModAdamMax("ModRAdamMax", nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "ModAdamSum" ) {
                            initializeBMStateForModAdamSum(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "Adam" ) {
                            initializeBMStateForAdam(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "RAdam" ) {
                            initializeBMStateForRAdam(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if ( newOptimizer == "Adadelta" ) {
                            initializeBMStateForAdadelta(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          } else if  ( newOptimizer == "RPROP-LR" || newOptimizer == "MF" ) {
                            initializeBMStateForRPROPLR("RPROP-LR", nextBMState.bmInteractions, learningRate = nextBMState.learningRateForRPROPLR )
                          } else {
                            sys.error("Not supported: " + newOptimizer) 
                            initializeBMStateForAdam(newOptimizer, nextBMState.bmInteractions, nextBMState.learningRate )
                          }

                      }
                      if  ( newOptimizer == "RPROP-LR" || newOptimizer == "MF" ) {
                            stderr.print("If necessary, please modify your program to specify the parameters for the new optimizer.\n" )
                      } else {
                            assert( initialState.learningRate.maxLR > 0.0, "Please modify the program to specify the parameters for the new optimizer.\n" )
                      }
                      val newBetaV = if( param.contains("betaV") ) {
                                     param("betaV").toDouble
                                 } else {
                                     initialState.betaV
                                 }
                      val newBetaM = if( param.contains("betaM") ) {
                                     param("betaM").toDouble
                                 } else {
                                     initialState.betaM
                                 }

                      outTE.print("# New optimizer= %s  betaV= %g  betaM= %g  maxLR= %g  maxLREpochs= %g  stepOffset= %d  Total minLearningSteps= %d  minLearningStepsForBestKL= %d".format(newOptimizer, newBetaV, newBetaM, newMaxLR, newMaxLREpochs, newStepOffset, minLearningSteps_rev, minLearningStepsForBestKL_rev ) )
                      if ( newBetaLAP != nextBMState.betaLAP ) outTE.print("  betaLAP: %g".format(newBetaLAP) )
                      outTE.print(" ;  The step number for new gradientDescentMethod are virtually reset to 0.\n" )


                      outKL.print("# New optimizer= %s  betaV= %g  betaM= %g  maxLR= %g  maxLREpochs= %g  stepOffset= %d  Total minLearningSteps= %d  minLearningStepsForBestKL= %d".format(newOptimizer, newBetaV, newBetaM, newMaxLR, newMaxLREpochs, newStepOffset, minLearningSteps_rev, minLearningStepsForBestKL_rev ) )
                      if ( newBetaLAP != nextBMState.betaLAP ) outKL.print("  betaLAP: %g".format(newBetaLAP) )
                      outKL.print(" ;  The step number for new gradientDescentMethod are virtually reset to 0.\n" )

                      nextBMState.copy(gradientDescentMethod = initialState.gradientDescentMethod,
                                       grad = initialState.grad,
                                       v = initialState.v, m = initialState.m,  
                                       betaV = newBetaV, betaM = newBetaM,
                                       learningRate = initialState.learningRate.newLearningRate( maxLR = newMaxLR, maxLREpochs = newMaxLREpochs ) ,     
                                       learningRateForRPROPLR = initialState.learningRateForRPROPLR,
                                       learningRates = initialState.learningRates,
                                     //if stepOffset = 0, learningRate = nextBMState.learningRate.newLearningRate( stepOffset = nextBMState.step.toDouble ),
                                       stepOffset = newStepOffset, betaLAP = newBetaLAP )

         //Error: } else if((bmStateWithMCSamples.step + 1) % logInterval == 0 && bmStateWithMCSamples.step - bmStateWithMCSamples.stepOffset + 1 >= minLearningSteps ) {
         //       } else if((bmStateWithMCSamples.step + 1) % logInterval == 0 && bmStateWithMCSamples.step + 1 >= minLearningSteps ) {
         //       } else if((bmStateWithMCSamples.step + 1) % logInterval == 0 && bmStateWithMCSamples.step + 1 == minLearningSteps ) {
                  } else if ( param.size > 0 ) {

                    if ( minLearningSteps_rev != minLearningSteps || minLearningStepsForBestKL_rev != minLearningStepsForBestKL || newMaxLR != nextBMState.learningRate.maxLR || newMaxLREpochs != nextBMState.learningRate.maxLREpochs || newStepOffset != nextBMState.stepOffset || newBetaLAP != nextBMState.betaLAP ) {
                        outTE.print("# New minLearningSteps= %d  minLearningStepsForBestKL= %d  maxLR= %g  maxLREpochs= %g  stepOffset= %d".format(
                                     minLearningSteps_rev, minLearningStepsForBestKL_rev, newMaxLR, newMaxLREpochs, newStepOffset) )
                        outKL.print("# New minLearningSteps= %d  minLearningStepsForBestKL= %d  maxLR= %g  maxLREpochs= %g  stepOffset= %d".format(
                                     minLearningSteps_rev, minLearningStepsForBestKL_rev, newMaxLR, newMaxLREpochs, newStepOffset) )
                      
                        if ( newBetaLAP != nextBMState.betaLAP ) {
                          outTE.print("  betaLAP: %g\n".format(newBetaLAP) )
                          outKL.print("  betaLAP: %g\n".format(newBetaLAP) )
                        } else {
                          outTE.print("\n")
                          outKL.print("\n")
                        }

                        nextBMState.copy( learningRate = nextBMState.learningRate.newLearningRate( maxLR = newMaxLR, maxLREpochs = newMaxLREpochs ),
                                   stepOffset = newStepOffset, betaLAP = newBetaLAP
                                    )
                    } else {
                        nextBMState
                    }
                  } else {
                    nextBMState
                  }
              }

            //if(bmStateWithMCSamples.step + 1 >= math.max(minLearningStepsForBestKL, minLearningSteps_rev) && nSamples >= maxNSamples ) 
              val (newInteractions, newBMState) = 
                  BM.newPhiThroughhJ(nextBMState_rev, fiaForBM, gauge) 
                //BM.newPhiThroughhJ(nextBMState_rev, this.fia, gauge) 
          
              iterateBMSteps(runBM1Step,
                newInteractions,
                fiaForBM, fijabForBM, faForBM,

                miniBSize,

                nextLastAndCurrentFullBatch,
                nextPrevMiniBatches,

                nextStartPos,
                nSamples,
                lambdaPhi, lambdaPhij,

                newBMState, newBestKLs, newBestKL,

                minLearningStepsForBestKL_rev, minLearningSteps_rev, maxNoLearnings,

                nMiniBatchesForEnsembleAve = nMBsForEnsembleAve,

                logInterval,
                newOptimizer)
          }

        } // End of iterateBMSteps

      /* moved before printParameters:
        val logInterval_rev = if (stepsPerEpoch <= 1.0) {
                                 if ( logInterval > 0 ) logInterval else 1
                              } else if ( stepsPerEpoch <= logInterval ) { 
                                 (stepsPerEpoch * (logInterval / stepsPerEpoch).floor ).toInt
                              } else {
                                 stepsPerEpoch.toInt
                              }
      */

        val initialInteractions = { 
                val hJ =
                  if (optionInitialInteractions == null || optionInitialInteractions == None ) {
                        initializehJ (this.fia, this.effectiveNConfigs, pseudoNCountsForhia, sigma_J = sigma_initial_J)
                      //initializehJ (fiaForBM, this.effectiveNConfigs, 0.0, sigma_J = sigma_initial_J)
                      //initializehJ (proposedPia, this.effectiveNConfigs, 0.0, sigma_J = sigma_initial_J)
                  } else {
                        val interactions = optionInitialInteractions.get
                        assert(interactions.hia.size == this.fia.size && interactions.hia(0).size == this.fia(0).size )
                        assert(interactions.jijab.size == this.fijab.size && interactions.jijab(0).size == this.fijab(0).size )
                        ( interactions.hia, interactions.jijab )
                  }
                if ( gauge == "unused" || gauge == "ungauged" ) {
                  MCMC.Interactions(hJ._1,  hJ._2, gauge)
                } else {
                //val phiPhij = BM.hJToPhi (hJ._1,  hJ._2, this.fia)
                  val phiPhij = BM.hJToPhi (hJ._1,  hJ._2, fiaForBM)
                  val bmInteractions = BMInteractions(phiPhij._1, phiPhij._2)

                //val (newInteractions, newBMInteractions) = BM.newPhiThroughhJ(bmInteractions, this.fia, gauge)
                  val (newInteractions, newBMInteractions) = BM.newPhiThroughhJ(bmInteractions, fiaForBM, gauge)
                  newInteractions
                } 
        }

        val (interactions, bmState) = {
                         //val (inter, bms) = initialize(this.fia,
                           val (inter, bms) = BM.initialize(
                                  fiaForBM,
                                  initialInteractions,
                                  optMethod, 
                                  learningRate_rev, 
                                  learningRateForRPROPLR, 
                                //minLearningRate, maxLearningRate, 
                                //rateDecrease, rateIncrease, 
                                  betaV, betaM, eps, 0,
                                  betaLAP, betaLAPForKL )

                           if(gauge == inter.gauge ) {
                                        (inter, bms)
                           } else {
                                        assert( gauge == inter.gauge )
                                      //BM.newPhiThroughhJ(bms, this.fia, gauge) 
                                        BM.newPhiThroughhJ(bms, fiaForBM, gauge) 
                           }
        }

        val newInteractions = interactions
        
        val  ( 
               nativeConfigs: ArraySeq[MCMC.NativeConfig],
               lastAndCurrentFullBatch: Vector[BMState],
               prevMiniBatches: Vector[BMState],
               startPos: Int,
               newBMState: BMState
             ) = preFirstBatch(miniBatchSize: Int,
                                  initialConfigurationIndex: ArraySeq[Int],
                                  initialConfigurations: ArraySeq[MCMC.Configuration],
                                  sampleConfigurations: ArraySeq[MCMC.Configuration],
                                  sampleIDs: ArraySeq[String],
                                  bmState: BMState )

        assert( fullBatchSize == nativeConfigs.size )  
      //assert( fullBatchSize == lastAndCurrentFullBatch.map(bms => bms.initialMCStates.size).sum )
        assert( fullBatchSize == lastAndCurrentFullBatch.map{bms => math.max(bms.initialMCStates.size, bms.independentSamplings.size) }.sum )

      //val miniBSize = miniBatchSizeCorrected(miniBatchSize, fullBatchSize)

        val nSamples_rev = if ( nSamples <= 0 ) {
                               val n = initialNSamples + (incrementSamples * bmState.step).floor.toInt
                               if ( n > maxNSamples ) maxNSamples else n
                           } else {
                               nSamples
                           }

        val nMiniBatchesForAve = nMiniBatchesForEnsembleAverage(miniBSize, 
                                     math.max(fullBatchSize, (50000.0 / nSamples_rev).floor.toInt ), nMiniBatchesForEnsembleAve )

        if ( betaLAPForKL == 0.0 )  { 
            outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and those over %d minibatch ensembles with pseudoNCountsForKL= %g\n".format(
                    nMiniBatchesForEnsembleAve, pseudoNCountsForKL ) )
        } else {
            outKL.print("# step\tnNonEquil  nExtendedIterations of MC\t<KLpi>\t<KLpij> over ensemble and those over %d minibatch ensembles with pseudoNCountsForKL= %g  KL with the leaky rate %g for P\n".format(
                    nMiniBatchesForEnsembleAve, pseudoNCountsForKL, betaLAPForKL ) )
        }

        outTE.print("# TE is calculated in the Ising gauge for comparison.\n#\n") 
        outTE.print("# For <TE>_m and <(TE -<TE>)^2>_m, %d minibatch ensembles are employed for averaging.\n#\n".format(nMiniBatchesForEnsembleAve)) 

        val results  = iterateBMSteps(runBM1Step,
                          newInteractions,
                          fiaForBM, fijabForBM, faForBM,

                          miniBSize,      //if ( miniBatchSize > 0 ) miniBatchSize else initialConfigurations_for_miniBatch.size,
                       // initialMCStates: Vector[ArraySeq[MCMC.State]],

                          lastAndCurrentFullBatch,     //Vector[BMState],
                          prevMiniBatches,             //Vector[BMState],

                          startPos,
                          nSamples_rev,
                          lambdaPhi, lambdaPhij,
                          newBMState, bestKLs, bestKL,
                          minLearningStepsForBestKL_rev, minLearningSteps_rev, maxNoLearnings,
                          nMiniBatchesForAve,

                          logInterval_rev,
                          optMethod  )

        outKL.close()
        outTE.close()

        results

    } // End of RunBM


  } // End of class BM

} // End of org.sanzo.potts
