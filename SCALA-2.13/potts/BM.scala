
import org.{sanzo => miyazawa}

package org.sanzo.potts {

  import scala.collection.parallel.CollectionConverters._	// for par in scala 2.13
  import scala.math.Ordering.Double.TotalOrdering		// for scala 2.13

  import breeze.linalg
  import breeze.linalg.DenseMatrix
  import breeze.linalg.DenseVector
  import breeze.linalg.{max, sum}
  import breeze.numerics.{exp, log}
  import breeze.stats.meanAndVariance
  import breeze.stats.distributions.Uniform

  import breeze.stats.distributions.Process
  import breeze.stats.distributions.Rand

  //import spire.math.Integral
  import spire.syntax.cfor._

  import java.io.File 
  import java.io.PrintStream

  import scala.sys.process.stdin
  import scala.sys.process.stdout
  import scala.sys.process.stderr
  import scala.io.Source
  import scala.annotation.tailrec
  import scala.collection.parallel.mutable.ParArray
  import scala.collection.mutable.HashMap
  import scala.collection.immutable.Queue
  import scala.math
  import scala.util.matching.Regex

  //import miyazawa.sequence.SequenceArray.pairIndex
  import miyazawa.potts.MCMC.pairIndex
  import miyazawa.potts.MCMC.{inversePairIndex, fromPairIndex}

  object MCMC {

    case class State(val configuration: Array[Byte], val energy: Double, val kT: Double, val annealingRate: Double, val step: Int, val finalkT: Double ) {
	require ( kT >= finalkT )
    }
    case class Interactions(val hia: Array[DenseVector[Double]], val Jijab: Array[DenseMatrix[Double]], gauge: String = "ungauged")

    case class ProposedDistributions(val proposedDistributionsAtAllUnits: Array[DenseVector[Double]],
	val logProposedDistributionsAtAllUnits: Array[DenseVector[Double]] ,
	val log1_ProposedDistributionsAtAllUnits: Array[DenseVector[Double]] )

    case class EnsembleAverages(val pia: Array[DenseVector[Double]], val pijab: Array[DenseMatrix[Double]])

    type IndependentMC        = ParArray[Iterator[State]]	//type   IndependentMC[T] = ParArray[T]
    type IndependentSamplings = ParArray[Array[State]]		//type   IndependentSamplings[T] = ParArray[T]

    def pairIndex(i: Int,j: Int) = {
                if(i > j)
                        (i * (i - 1) / 2 + j)
                else if ( j > i )
                        (j * (j - 1) / 2 + i)
                else {
                        sys.error("*** Error: Pairwise index must be for i != j\n")
                        -1
                }
    }

    def pairIndex(i: Int,j: Int, ai: Int, aj: Int) = {
                if(i > j)
                        ((i * (i - 1) / 2 + j), ai, aj)
                else if ( j > i )
                        ((j * (j - 1) / 2 + i), aj, ai)
                else {
                        sys.error("*** Error: Pair index must be for i != j\n")
                        (-1, -1, -1)
                }
    }

    def inversePairIndex(ij: Int) = {
		val i = {
			  val i = math.sqrt(ij * 2.0).floor.toInt
			  if ( ij - (i * ( i - 1)) / 2 >= i ) 
				i + 1
			  else
				i
			}
		val j = ij - (i * ( i - 1)) / 2  
		require( pairIndex(i, j) == ij )

		(i, j)
    }

    def fromPairIndex(ij: Int) = inversePairIndex(ij)


    def readNSamples(lines: collection.Iterator[String] ) = {
	val (lines0, lines1) = lines.duplicate

	val nHeaderLine = ("""^#[ \t]*Effective_number_of_samples:[ \t]*([0-9.eE+]+).*$""").r
	val lines11 = lines1.dropWhile(line => ! nHeaderLine.findFirstMatchIn(line).nonEmpty )
	val effectiveNumberOfSamples =
	  if ( lines11.hasNext ) {
	    val line = lines11.next
	    val nHeaderLine(effectiveNumberOfSamplesString) = line
	    val effectiveNumberOfSamples = effectiveNumberOfSamplesString.toDouble
	    effectiveNumberOfSamples
	  } else {
	    0.0
	  }

	(effectiveNumberOfSamples, lines0)
    }

    def readN(lines: collection.Iterator[String]) = {
	val (lines0, lines1) = lines.duplicate
	val nHeaderLine = ("""^#[ \t]*n_units: *([0-9]+)[ \t]*n_states_of_unit:[ \t]*([0-9]+).*$""").r
	val lines11 = lines1.dropWhile(line => ! nHeaderLine.findFirstMatchIn(line).nonEmpty )
	val (nUnits, nStatesOfUnit, stateOrderString) =
	  if ( lines11.hasNext ) {
	    val line = lines11.next
	    val nHeaderLine(nUnitsString, nStatesOfUnitString) = line
	    val nUnits = nUnitsString.toInt
	    val nStatesOfUnit = nStatesOfUnitString.toInt

	    val sHeaderLine = ("""^.*[ \t]+State_Order_String:[ \t]*([^ \t]+).*$""").r
	    val stateOrderString: String = if(sHeaderLine.findFirstMatchIn(line).nonEmpty) {
			  val sHeaderLine(s) = line
			  s
			} else {
			  ""	// null
			}
	    (nUnits, nStatesOfUnit, stateOrderString)
	  } else {
	    (0, 0, "")	// (0, 0, null)
	  }

	(nUnits, nStatesOfUnit, stateOrderString, lines0)
    }

    def readhOrPi(lines: collection.Iterator[String], nUnits: Int, nStatesOfUnit: Int, whichhOrPi: String ) = {
	val lines1 = lines
	
	val commentline = """^#.*""".r
        val hPiHeaderLine = ("""^""" + whichhOrPi + """[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9eE.+-]+).*$""").r
	val hPiLine = hPiHeaderLine
	val lines11 = lines1.dropWhile(line => ! hPiLine.findFirstMatchIn(line).nonEmpty )

	val h = (new Array[DenseVector[Double]](nUnits)).map( hi => new DenseVector[Double](nStatesOfUnit) )

	@annotation.tailrec
	def processhPiLine(lines: collection.Iterator[String], h: Array[DenseVector[Double]], n: Int): Int = {
		  if(lines.hasNext) {
		    val line = lines.next
		    if ( (hPiLine.findFirstMatchIn(line) ).nonEmpty ) {
			val field = line.split("""[ \t]+""")
			val ir = field(1).toInt + 1
			val ia_orig = field(2).toInt

			if ( ir <= nUnits ) {
			//cfor(0)(i => i < nAADelTypes, i => i + 1) ( i => {
                        //val i = ia_orig
			//	val ia = aaDelOrderMap.compoundOrderNumber( aaOrderString(i), -1 )	/* error generator */
			//	h(ir - 1)(ia) = field(i + 1).toDouble
			  val ia = ia_orig
			  h(ir - 1)(ia) = field(3).toDouble
			//} )
			  if( n + 1 >= nUnits * nStatesOfUnit) { 
				n + 1
			  } else {
				processhPiLine(lines, h, n + 1)
			  }
			} else {
			  sys.error("*** Error: the number of residues is larger than %d\n".format(nUnits))
			  - n
			}
		    } else if ( (commentline.findFirstMatchIn(line) ).nonEmpty ) {
			processhPiLine(lines, h, n)
		    } else {
			if ( n > 0 )
			  sys.error("*** Error in %s lines: %s".format(whichhOrPi, line))
			//- n
			processhPiLine(lines, h, n)
		    }
		  } else {
			sys.error("*** Error: no enough number of %s lines.\n".format(whichhOrPi))
			- n
		  }
	}
	val nh = processhPiLine(lines, h, 0)
	if( nh < nUnits )
		sys.error("*** Error: no enough number of %s lines: %d\n".format(whichhOrPi, nh))
	
	(h, lines)
    }

    def readJOrPij(lines: collection.Iterator[String], nUnits: Int, nStatesOfUnit: Int, whichJOrPij: String ) = {
	val lines1 = lines
	
	val commentline = """^#.*""".r
        val JPijHeaderLine = ("""^""" + whichJOrPij + 
		"""[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]([0-9]+)[ \t]+([0-9]+)[ \t]+([0-9eE.+-]+).*$""").r
	val JPijLine = JPijHeaderLine
	val lines11 = lines1.dropWhile(line => ! JPijLine.findFirstMatchIn(line).nonEmpty )

	val nPairs = (nUnits * (nUnits - 1)) / 2
	val J = (new Array[DenseMatrix[Double]](nPairs)).map( Jij => new DenseMatrix[Double](nStatesOfUnit, nStatesOfUnit) )

	@annotation.tailrec
	def processJPijLine(lines: collection.Iterator[String], J: Array[DenseMatrix[Double]], n: Int): Int = {
		  if(lines.hasNext) {
		    val line = lines.next
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
			//	val ia = aaDelOrderMap.compoundOrderNumber( aaOrderString(i), -1 )	/* error generator */
			val (k, a, b) = pairIndex(jr -1, ir -1, ja, ia)
			//	J(k)(a) = field(i + 3).toDouble
			J(k)(a, b) = field(5).toDouble
			//} )
			  val nn = (nUnits * ( nUnits - 1)) / 2 * nStatesOfUnit * nStatesOfUnit
			  if( n + 1 >= nn ) { 
				//(ir == nUnits && jr == nUnits -1) ||
				//(ir == nUnits - 1 && jr == nUnits) ) {
				n + 1
			  } else {
				processJPijLine(lines, J, n + 1)
			  }
			} else {
			  sys.error("*** Error: the number of residues is larger than %d\n".format(nUnits))
			  - n
			}
		    } else if ( (commentline.findFirstMatchIn(line) ).nonEmpty ) {
			processJPijLine(lines, J, n)
		    } else {
			if ( n > 0 )
			  sys.error("*** Error in %s lines: %s".format(whichJOrPij, line))
			//- n
			processJPijLine(lines, J, n)
		    }
		  } else {
			sys.error("*** Error: no enough number of %s lines.\n".format(whichJOrPij))
			- n
		  }
	}
	val nJ = processJPijLine(lines, J, 0)
	if( nJ < nUnits * (nUnits - 1) / 2 * nStatesOfUnit)
		sys.error("*** Error: no enough number of %s lines: %d\n".format(whichJOrPij, nJ))
	
	(J, lines)
    }

    def readhJ(lines: collection.Iterator[String] ) = {

	val (nUnits, nStatesOfUnit, stateOrderString, lines1) = readN(lines)
    	val ( h, lines2 ) = readhOrPi(lines1, nUnits, nStatesOfUnit, "h" )
    	val ( j, lines3 ) = readJOrPij(lines2, nUnits, nStatesOfUnit, "J" )

	(h, j, stateOrderString)
    }

    def readPiPij(lines: collection.Iterator[String] ) = {

	val (nUnits, nStatesOfUnit, stateOrderString, lines1) = readN(lines)
    	val ( pia, lines2 ) = readhOrPi(lines1, nUnits, nStatesOfUnit, "Pi" )
    	val ( pijab, lines3 ) = readJOrPij(lines2, nUnits, nStatesOfUnit, "Pij" )

	(pia, pijab, stateOrderString)
    }

    def printNSamples(out: PrintStream, effectiveNumberOfSamples: Double) = {
	out.print("# Effective_number_of_samples: %g\n".format( effectiveNumberOfSamples ) )
    }

    def printhOrPi(outhOrPi: PrintStream, stateOrderString: String, 
		hOrPi: Array[DenseVector[Double]], whichhOrPi: String ) = {

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
		  //	stateOrderString(ia), hOrPi(i)(ia) ) )
		  //} else
		  {
		    outhOrPi.print("%s %d %d %g\n".format(whichhOrPi, i,
			ia, hOrPi(i)(ia) ) )
		  }
		} )
	} ) 
    }

    def printhOrPi(outhOrPi: PrintStream, stateOrderString: String, 
		hOrPi: Array[Array[Double]], whichhOrPi: String ): Unit = {

	printhOrPi(outhOrPi, stateOrderString,
		hOrPi.map( x => DenseVector(x) ), whichhOrPi)
     }

    def printJOrPij(outJOrPij: PrintStream, stateOrderString: String, 
		JOrPij: Array[DenseMatrix[Double]], whichJOrPij: String ) = {

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
		  //	stateOrderString(ia), stateOrderString(jb), interactions.Jijab(ij)(a, b) ) )
		  //} else {
		    outJOrPij.print("%s %d %d %d %d %g\n".format(whichJOrPij, i, j, 
			ia, jb, JOrPij(ij)(a, b)   ))
		  //}
		} )
		} )
	  } ) 
	} ) 
    }

    def printJOrPij(outJOrPij: PrintStream, stateOrderString: String, 
		JOrPij: Array[Array[Array[Double]]], whichJOrPij: String ): Unit = {
    	printJOrPij(outJOrPij, stateOrderString, 
		JOrPij.map(pij => {
			//DenseMatrix( pij: _*)		// for 2.11
			  new breeze.linalg.DenseMatrix( pij(0).size, pij.size, pij.flatten).t	// for 2.13
		} ), 
		whichJOrPij )
    }

    def printInteractions(outhJ: PrintStream, stateOrderString: String, 
		interactions: Interactions ) = {

	val h = interactions.hia
	val J = interactions.Jijab

	outhJ.println("# Gauge: " + interactions.gauge )
    	printhOrPi(outhJ, stateOrderString, h, "h" )
    	printJOrPij(outhJ, stateOrderString, J, "J")
    }

    def printhJ(outhJ: PrintStream, stateOrderString: String, 
		hia: Array[DenseVector[Double]],
		Jijab: Array[DenseMatrix[Double]] ) = {

    	printhOrPi(outhJ, stateOrderString, hia, "h" )
    	printJOrPij(outhJ, stateOrderString, Jijab, "J")
    }

    def printhJ(outhJ: PrintStream, stateOrderString: String, 
		hia: Array[Array[Double]],
		Jijab: Array[Array[Array[Double]]] ) = {

    	printhOrPi(outhJ, stateOrderString, hia, "h" )
    	printJOrPij(outhJ, stateOrderString, Jijab, "J")
    }

    def printPiPij(outPiPij: PrintStream, stateOrderString: String, 
		pia: Array[DenseVector[Double]],
		pijab: Array[DenseMatrix[Double]] ) = {

    	printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
    	printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }

    def printPiPij(outPiPij: PrintStream, stateOrderString: String, 
		pia: Array[Array[Double]],
		pijab: Array[Array[Array[Double]]] ) = {

    	printhOrPi(outPiPij, stateOrderString, pia, "Pi" )
    	printJOrPij(outPiPij, stateOrderString, pijab, "Pij")
    }

    def printIndependentMCsamplings(outSamples: PrintStream, stateOrderString: String, 
		independentSamplings: IndependentSamplings ) = {

	val nSamplings = independentSamplings.size
      //val nSamples = independentSamplings(0).size
	val nUnits = independentSamplings(0)(0).configuration.size

    //State(val configuration: Array[Byte], val energy: Double, val kT: Double, val annealingRate: Double, val step: Int, finalkT: Double = 1.0 )
	if( stateOrderString == null || stateOrderString.size == 0) {
	  //outSamples.print("# State_Order_String: %s\n".format(stateOrderString) )
		cfor(0)(i => i < nSamplings, i => i + 1)( i => {
			//outSamples.print("# Independent_Sampling: %d / %d\n".format(i, nSamplings))
			val nSamples = independentSamplings(i).size
			cfor(0)(j => j < nSamples, j => j + 1)( j => {
				val state = independentSamplings(i)(j)
	  		      //outSamples.print("> %d/%d_%d/%d  sampling=%d/%d sample=%d/%d step=%d E=%g kT=%g annealing_rate=%g\n".format(
			      //	i, nSamplings, j, nSamples, 
			      //        i, nSamplings, j, nSamples, state.step, state.energy, state.kT, state.annealingRate) )
				cfor(0)(k => k < nUnits, k => k + 1)( k => {
					outSamples.print("%d ".format( state.configuration(k) ) ) } )
				outSamples.print("\n")
			} )
		} )
	  
	} else {
		cfor(0)(i => i < nSamplings, i => i + 1)( i => {
			val nSamples = independentSamplings(i).size
			cfor(0)(j => j < nSamples, j => j + 1)( j => {
				val state = independentSamplings(i)(j)
	  			outSamples.print("> %d/%d_%d/%d  sampling=%d/%d sample=%d/%d step=%d E=%g kT=%g annealing_rate=%g\n".format(
					i, nSamplings, j, nSamples,
					i, nSamplings, j, nSamples, state.step, state.energy, state.kT, state.annealingRate) )
				outSamples.print("%s\n".format( state.configuration.map(s 
					=> stateOrderString(s).toString).reduceLeft(_ + _) ) )
			} )
		} )
	}
	
    }

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

    def proposeState(currentState: Int, proposedDistributionAtUnit: DenseVector[Double]): Byte = {
	// x <= Uniform(x, y) < y; see breeze.stats.distributions.{Uniform, RandBasis.uniform}
	// but it is better to code 	u = Uniform(x, y) if u != y
        def uniformGTxLEy(x: Double, y: Double): Rand[Double] = 
	  for {
		r <- Uniform(x, y)	if r < y
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

    def totalE(configuration: Array[Byte],
	interactions: Interactions ) = {

	val nUnits = configuration.size
	val hia = interactions.hia
	val Jijab = interactions.Jijab

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
	- ne.sum
	*/

	val ne = Array.range(0, nUnits).map ( i => {
	  	    val ia = configuration(i)
		    val nei = Array.range(0, nUnits).map ( j => {
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
	- ne.sum
    }

    def deltaE(configuration: Array[Byte], 
	interactions: Interactions,
	position: Int, proposedState: Byte ) = {

	val i = position
	val currentState = configuration(i)
	val nUnits = configuration.size

	val hia = interactions.hia
	val Jijab = interactions.Jijab
	val dNE = // (new Array[Double](nUnits)).zipWithIndex.map ( dNEj => {
		  //   val j = dNEj._2 
	    Array.range(0, nUnits).map ( j => {
	  	if ( i != j ) { 
	    		val jb = configuration(j)
	    		val (ij, a, b) = pairIndex(i, j, currentState, jb)
	    		val (ji, c, d) = pairIndex(i, j, proposedState, jb)
	    		Jijab(ij)(c, d) - Jijab(ij)(a, b)
	  	} else {
	  		hia(i)(proposedState) - hia(i)(currentState)
	  	}
		
	    } )  
	- dNE.sum 
    }

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

    def deltaEs(configuration: Array[Byte], 
	interactions: Interactions,
	position: Int) = {

	val i = position
	val currentState = configuration(i)
	val nUnits = configuration.size
	val nStates = interactions.hia(0).size

	val hia = interactions.hia
	val Jijab = interactions.Jijab
	val dNEs = Array.range(0, nUnits).map ( j => {
		  	if ( i != j ) { 
		    		val jb = configuration(j)
				val (ij, ia, ja) = pairIndex(i, j, 0, 1)
				Array.range(0, nStates).map ( s => {
					    if ( ia == 0 ) {
		    				Jijab(ij)(s, jb) - Jijab(ij)(currentState, jb)
					    } else {
		    				Jijab(ij)(jb, s) - Jijab(ij)(jb, currentState)
					    }
					  } )
		 	} else {
				Array.range(0, nStates).map ( s => {
	  				hia(i)(s) - hia(i)(currentState)
				  } )
  			}
		
	    	} )  

	Array.range(0, nStates).map ( s => {
		val dNE = Array.range(0, nUnits).map ( j => { 
				dNEs(j)(s)
			} )
		- dNE.sum
		} )
    }

    // x <= Uniform(x, y) < y
    def kernelMH1Step(
	currentState: State,
	proposedDistributions: ProposedDistributions,
	//proposedDistributionsAtAllUnits: Array[DenseVector[Double]],
	//logProposedDistributionsAtAllUnits: Array[DenseVector[Double]] ,
	//log1_ProposedDistributionsAtAllUnits: Array[DenseVector[Double]] ,
	interactions: Interactions,
	siteLocation: Int = -1	// In the case of -1, the siteLocation will be randomly determined.
				 ): Rand[State] = {

      val energy = currentState.energy
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val proposedDistributionsAtAllUnits = proposedDistributions.proposedDistributionsAtAllUnits
      val nUnits = proposedDistributionsAtAllUnits.size
      val newkT = if ( kT <= finalkT ) finalkT else {
		     val t = currentState.annealingRate * (kT/finalkT - 1.0)
		     if ( t < 0.00001 ) finalkT else (t + 1.0) * finalkT
		}
      for {
	u <- Uniform(0.0, 1.0) if (u < 1.0)

        configuration = currentState.configuration	//currentState.configuration.clone
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
      //   		proposeState(maxP - uniformS, currentS, proposedDistributionsAtAllUnits(position))
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
		State(configuration, energy + dE, newkT, 
			currentState.annealingRate, currentState.step + 1, currentState.finalkT)
		//println("position: %d , state:  %d => %d".format(position, currentS, configuration(position) ));
		//proposedState
	    } else {
		State(configuration, energy, newkT, 
			currentState.annealingRate, currentState.step + 1, currentState.finalkT)
	    }
    }

    def kernelMH(
	currentState: State,
	proposedDistributions: ProposedDistributions,
	//proposedDistributionsAtAllUnits: Array[DenseVector[Double]],
	//logProposedDistributionsAtAllUnits: Array[DenseVector[Double]] ,
	//log1_ProposedDistributionsAtAllUnits: Array[DenseVector[Double]] ,
	interactions: Interactions ): Rand[State] = {

      val nUnits = proposedDistributions.proposedDistributionsAtAllUnits.size

      val kernel1Step = kernelMH1Step(_: State, proposedDistributions, interactions) 
    //val dupCurrentState = currentState.copy(configuration = currentState.configuration.clone)
      /* The following must be within the for block;  the thread stops before <-.
      val iter = markovChain(currentState, kernel1Step)
      val state = if( nUnits > 1 )
      			iter.drop( nUnits - 1 ).next
		  else
			iter.next
      */
      for{ r <- Uniform(0.0, 1.0)	// dummy statement
    	   dupCurrentState = currentState.copy(configuration = currentState.configuration.clone)
	   mc = markovChain(dupCurrentState, kernel1Step)
	   state = mc.drop( nUnits - 1 ).next
		} yield { state }
    }

    def kernelMultiBlockMH(
	currentState: State,
	proposedDistributions: ProposedDistributions,
	interactions: Interactions ): Rand[State] = {

      @annotation.tailrec
      def recCallKernel(site: Int, state: State): State = {
	if (site >= 0) {
    	  val kernel = kernelMH1Step(_: State, proposedDistributions, interactions, siteLocation = site)
	  val newState = kernel(state).draw
	  recCallKernel(site - 1, newState)
	} else {
	  state
	}
      }

    //val configuration = currentState.configuration.clone
      val nUnits = currentState.configuration.size
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val newkT = if ( kT <= finalkT ) finalkT else {
		     val t = math.pow(currentState.annealingRate, nUnits.toDouble) * (kT/finalkT - 1.0)
		     if ( t < 0.00001 ) finalkT else (t + 1.0) * finalkT
		}
      for {
	xx <- Uniform(0.0, 1.0)	// dummy statement

        dupCurrentState = currentState.copy(configuration = currentState.configuration.clone, kT = newkT, finalkT = newkT)
	state = recCallKernel(nUnits - 1, dupCurrentState)
      } yield { 
	/**/
	state.copy(finalkT = finalkT)
	/**/
      }
    }

    def kernelGibbs(
	currentState: State,
	interactions: Interactions ): Rand[State] = {
    //val configuration = currentState.configuration.clone
      val energy = currentState.energy
      val nUnits = currentState.configuration.size
      val kT = currentState.kT
      val finalkT = currentState.finalkT
      val newkT = if ( kT <= finalkT ) finalkT else {
		     val t = math.pow(currentState.annealingRate, nUnits.toDouble) * (kT/finalkT - 1.0)
		     if ( t < 0.00001 ) finalkT else (t + 1.0) * finalkT
		}
      for {
	xx <- Uniform(0.0, 1.0)	// dummy statement

        configuration = currentState.configuration.clone
        dE = Array.range(0, nUnits).map( i => {
		val dEs = deltaEs(configuration, interactions, i )
		val mindEs = dEs.min
		val unnormP = dEs.map( e => math.exp( - (e - mindEs ) / newkT ) )
		val z = unnormP.sum
		val p = unnormP.map( x => x / z )
		val pDVec = DenseVector(p)
		val stateChosen: Byte =  proposeState(-10, pDVec)
		configuration(i) = stateChosen		// Be caution
		dEs(stateChosen.toInt)
	      } )
        newE = energy + dE.sum
      } yield { 
	/**/
	State(configuration, newE, newkT, 
		currentState.annealingRate, currentState.step + nUnits, currentState.finalkT)
	/**/
      }
    }

    def markovChain(initial: State, kernel: State => Rand[State] ): Iterator[State] = {
 	//val init = State(initial.configuration.clone, initial.energy,
	//		 initial.kT, initial.annealingRate, initial.step, initial.finalkT)
	val init = initial.copy(configuration = initial.configuration.clone )
	//
  	// import breeze.stats.distributions.MarkovChain
	// MarkovChain(init)(kernel).steps			// for breeze 0.13.1
	// or
	val kern = kernel(_: State).draw()			// return T rather than Rand[T]
	val iter = scala.collection.Iterator.iterate(init)(kern)
	val initValue = iter.next		// remove the initial value.
	iter
	// or
	// val kern = kernel(_: State).draw()			// return T rather than Rand[T]
	// scala.collection.Iterator.iterate(init)(kern).drop(1)	// remove the initial value.
    }

   /*
    def sampling(markovProcess: Process[State], 
		nInitialIterations: Int, 
		everyNIterations: Int,
		nSamples: Int 
 		): Array[State] = {

	val iterator = markovProcess.steps.drop(nInitialIterations)
	val samples = (new Array[State](nSamples)).map(sample => {
		val s = iterator.drop(everyNIterations - 1).next	// Error ?
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
 		): Array[State] = {

	//val iterator = markovChain.drop(nInitialIterations)
	val iterator = if ( nInitialIterations > 0 ) {
		val i = (nInitialIterations - 1 ) / 2
		val j = nInitialIterations - 1 - i
		val state = markovChain.drop(i).next
		require( state.kT == state.finalkT )
		markovChain.drop(j)				// Error; see the manual for Iterator
		//
	} else {
		markovChain
	}
	val samples = (new Array[State](nSamples)).map(sample => {
		val s = iterator.drop(everyNIterations - 1).next		// Error ?
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
 		): Tuple2[Array[State], Iterator[State]] = {

	//val iterator = markovChain.drop(nInitialIterationsPerUnit)
	val iterator = if ( nInitialIterationsPerUnit > 0 ) {
		val i = (nInitialIterationsPerUnit - 1 ) / 2
		val j = nInitialIterationsPerUnit - 1 - i
		Range(0, i).foreach( x => { val drop = markovChain.next } )
		val state = markovChain.next
		require( state.kT == state.finalkT )
		Range(0, j).foreach( x => { val drop = markovChain.next } )
		markovChain
	      //Or
	      //val markovChain2 = markovChain.drop(i)
	      //val state = markovChain2.next
	      //require( state.kT == state.finalkT )
	      //val markovChain3 = markovChain2.drop(j)
	      //markovChain3
	} else {
		markovChain
	}
	val samples = (new Array[State](nSamples)).map(sample => {
		Range(0, everyNIterationsPerUnit - 1).foreach( x => { val drop = iterator.next } )
		val s = iterator.next
		s.copy(configuration = s.configuration.clone)
	      //State(s.configuration.clone, s.energy, s.kT, s.annealingRate, s.step, s.finalkT )
	} )
	// Or recursive call must be used.
	(samples, iterator)
    }

    def nNonEquilibrium(independentSamplings: IndependentSamplings ) = {

	val nRuns = independentSamplings.size
	val nSamples = independentSamplings(0).size
	val firstSamples = independentSamplings.map{ s => s(0).energy }.toArray
	val lastSamples = independentSamplings.map{ s => s(nSamples - 1).energy }.toArray
	val statFirstSamples = meanAndVariance(firstSamples)
	val statLastSamples = meanAndVariance(lastSamples)
	val avFirstSamples = statFirstSamples.mean
	val avLastSamples = statLastSamples.mean

	//val stdDifAv = math.sqrt( (statFirstSamples.variance + statLastSamples.variance) / nRuns )
	//val nDifAv = if ( (avFirstSamples - avLastSamples).abs > stdDifAv * 2.0) 1 else 0

	val stdDifAv = math.sqrt( (statFirstSamples.variance + statLastSamples.variance) )
	val difAv = firstSamples.zip(lastSamples).map( m => if ( (m._1 - m._2).abs > stdDifAv * 2.0) 1 else 0 )

	val nDifAv = difAv.take(5).sum

      //nDifAv

	// t-test for autocorrelation

	val intervalAC = nSamples - 10
	//val intervalAC = math.min(900, nSamples - 10)

      if ( intervalAC < 0 ) {
	nDifAv

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
				{	val (e1, e2) = x
					val corrMat = breeze.stats.corrcoeff ( DenseMatrix(e1, e2).t )
					corrMat(0, 1)
				} )
	val autoCorrelated  = independentSampledAutoCorr.map ( c => if (isCorrelated(c, nSamples - intervalAC)) 1 else 0 )
	val nCorr = autoCorrelated.take(5).sum
	/* */

	/*
	val nRunsArray = Array.range(0, nRuns).par
	val nACArray = Array.range(0, nSamples - intervalAC).par
	val independentSampledEnergies = independentSamplings.map( s => s.map( si => si.energy ) )
	val statIndependentSampledEnergies = independentSampledEnergies.map(s => meanAndVariance(s) )
	val independentSampledDEnergies = nRunsArray.map( s => independentSampledEnergies(s).map
						(ei => ei - statIndependentSampledEnergies(s).mean) )
	val independentSampledAutoCorr = nRunsArray.map( s => {
			val corr = nACArray.map( i => {
			  independentSampledDEnergies(s)(i) * independentSampledDEnergies(s)(i + intervalAC) 
			  } ).sum / (nSamples - intervalAC)
			val corrCoeff = corr / statIndependentSampledEnergies(s).variance
			corrCoeff * corrCoeff
		} )

	val autoCorrelated  = independentSampledAutoCorr.map ( c => if (isCorrelated(c, nSamples - intervalAC)) 1 else 0 )
	val nCorr = autoCorrelated.take(10).sum
	*/

	scala.sys.process.stderr.println(nDifAv, nCorr )

	nDifAv + nCorr
      }
    }

    def randomConfiguration(nUnits: Int, nStatesofUnit: Int,
	proposedDistributionsAtAllUnits: Array[DenseVector[Double]] ) = {

      // Probably: x <= Uniform(x, y) < y
      def uniformGExLTy(x: Double, y: Double): Rand[Double] = 
	for {
		r <- Uniform(x, y)	if r < y
        } yield { r }

      val conf = new Array[Byte](nUnits)
      if ( proposedDistributionsAtAllUnits == null ) {
	conf.map( x => { 
		( ( ( uniformGExLTy(0.0, 1.0).sample() * nStatesofUnit ).floor.toInt ) % nStatesofUnit ).toByte
		} )
      } else {
	//conf.zipWithIndex.map( x => {
	//	val position = x._2
	Array.range(0, nUnits).map( position => {
		if ( proposedDistributionsAtAllUnits(position).size <= 0 ) {
			( ( ( uniformGExLTy(0.0, 1.0).sample() * nStatesofUnit ).floor.toInt ) % nStatesofUnit ).toByte
		} else {
                	proposeState(currentState = -1 , proposedDistributionsAtAllUnits(position))
		}
	} )
      }
    }

    def  frequenciesAtUnitInSamples( independentSamplings: IndependentSamplings,		
		nStatesOfUnit: Int
		): Array[DenseVector[Double]] = {

	  val nUnits = independentSamplings(0)(0).configuration.size
	  //val samples = independentSamplings.flatten
	  val samples = independentSamplings.toArray.flatten

	  val pia = (new Array[DenseVector[Double]](nUnits)).map( _ => (new DenseVector[Double](nStatesOfUnit)) )

	  val weight = 1.0 / samples.size

	/* samples cannot be ParArray. */
	  samples.foreach( s => {
		s.configuration.zipWithIndex.foreach( c => pia(c._2)(c._1) += weight )
	  } )
	  pia
    }

    def pairwiseFrequenciesInSamples( independentSamplings: IndependentSamplings,		
		nStatesOfUnit: Int
		): Array[DenseMatrix[Double]] = {

	  val nUnits = independentSamplings(0)(0).configuration.size
	  val nPairs = ( nUnits * ( nUnits - 1)) / 2
	  //val samples = independentSamplings.flatten
	  val samples = independentSamplings.toArray.flatten

	  val pijab = (new Array[DenseMatrix[Double]](nPairs)).map(_ => new DenseMatrix[Double](nStatesOfUnit, nStatesOfUnit) )

	  val weight = 1.0 / samples.size

	/* samples cannot be ParArray. */
	  samples.foreach( s => {
		val configuration = s.configuration
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

  }

  class MCMC( val outputFiles: HashMap[String, File], 
	      val stateOrderString: String,
	      val interactions: MCMC.Interactions,
	      val proposedPia: Array[DenseVector[Double]] ) {

    import miyazawa.potts.MCMC._

    //assert( nUnits == proposedDistributionAtUnit.size )
    //assert( nStatesOfUnit == proposedDistributionAtUnit(0).size )
    val nUnits = interactions.hia.size
    val nStatesOfUnit = interactions.hia(0).size

    //val (proposedDistributionsAtAllUnits, 
    //	logProposedDistributionsAtAllUnits, log1_ProposedDistributionsAtAllUnits ) = 
    val proposedDistributions = this.proposeDistributions(proposedPia)

    def this(outputFiles: HashMap[String, File],
		stateOrderString: String,
		hia: Array[Array[Double]],
		Jijab: Array[Array[Array[Double]]],
		proposedPia: Array[Array[Double]] ) = {
      this(outputFiles, stateOrderString,
		MCMC.Interactions( hia.map( hi => DenseVector(hi) ), 
		//Jijab.map( Jij => DenseMatrix(Jij:_*) )  ), 						// for scala 2.11
		  Jijab.map( Jij => new DenseMatrix(Jij(0).size, Jij.size, Jij.flatten ).t )  ), 	// for scala 2.13
		proposedPia.map(pi => DenseVector(pi)) )
    }

    def this(outputFiles: HashMap[String, File], 
		stateOrderString: String,
		interactions: MCMC.Interactions ) = {
	this( outputFiles, stateOrderString, interactions,
		{
		  val hia = interactions.hia
		  val p = 0.6
		  val uniform = (1.0/ hia(0).size) * (1.0 - p)

          	  hia.map( hi => {
	     		val maxhi = breeze.linalg.max(hi)
	       		val exphi = breeze.numerics.exp(hi - maxhi)	//hi.map( math.exp(_ - maxhi) )
			val z = breeze.linalg.sum(exphi)
			exphi.map{ x => (x / z) * p + uniform }
          	  } )
		}
	)
    }

    def this(outputFiles: HashMap[String, File], 
		stateOrderString: String,
		hia: Array[Array[Double]],
		Jijab: Array[Array[Array[Double]]] ) = {
	this(outputFiles, stateOrderString, hia, Jijab,
		{
		  val p = 0.6
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

    def proposeDistributions(proposedPia: Array[DenseVector[Double]]) = {

	val p = 0.9 
	val eps =  (1.0 / proposedPia(0).size) * (1.0 - p)
	val proposedP = proposedPia.map ( pi => { 
				//val newpi = pi.map{ pia => pia * p + eps }
				val newpi = pi * p + eps
				val sum = breeze.linalg.sum(newpi)
				newpi.map{ pia => pia / sum }
				} )
	/*
	val proposedP = proposedPia.map ( pi => { 
			pi.map( _ => 1.0/pi.size )
		} )
	*/
	val logProposedP = proposedP.map(pi => pi.map(pia => math.log(pia)) )
	val log1_proposedP = proposedP.map(Pi => {
		Pi.map( Pia => math.log(1.0 - Pia) )
	  } )
	ProposedDistributions(proposedP, logProposedP, log1_proposedP)
    }

    def kernelMH(currentState: State ): Rand[State] =
	  MCMC.kernelMH(currentState, 
		this.proposedDistributions,
		//proposedDistributionsAtAllUnits, 
		//logProposedDistributionsAtAllUnits,
		//log1_ProposedDistributionsAtAllUnits,
		this.interactions )

    def kernelMultiBlockMH(currentState: State ): Rand[State] =
	  MCMC.kernelMultiBlockMH(currentState,
	 	this.proposedDistributions, this.interactions )

    def kernelGibbs(currentState: State ): Rand[State] =
	  MCMC.kernelGibbs(currentState, 
		this.interactions )

    def markovChain(initial: State, kernel: State => Rand[State] = kernelMultiBlockMH ): Iterator[State] =
	  MCMC.markovChain(initial, kernel)
	
    def runMC(  initialConfigurations: Array[Array[Byte]],
		nInitialIterationsPerUnit: Int = 100, 
		everyNIterationsPerUnit: Int = 5,
		nSamples: Int = 1000,
	      //nIndependentMC: Int = 10,
		initialT: Double = 1.2,
		annealingRate: Double = 0.99,
		finalT: Double = 1.0,
		maxExtendedIterations: Int = 10,
		kernel: State => Rand[State] = kernelMultiBlockMH
 		) = {

	val nInitialIterations = nInitialIterationsPerUnit * this.nUnits
	val everyNIterations = everyNIterationsPerUnit * this.nUnits

	val nIndependentMC = initialConfigurations.size

	val (kT, rate) = {
		val t = if ( initialT < finalT ) finalT else initialT 
		val maxRate = math.exp( (math.log(0.00001) - math.log(t / finalT - 1.0 + 0.00001)) / 
					(nInitialIterations / 5.0) )	/* (nInitialIterations / 10.0) ) */
		val r = if ( annealingRate > maxRate ) maxRate else  annealingRate 
		(t, r)
		}

	//val independentMarkovProcesses = (new Array[Process[State]](nIndependentMC) ).map( x => {
	//	val configuration = MCMC.randomConfiguration(this.nUnits, this.nStatesOfUnit)

	val initialStates =
		  initialConfigurations.map( x => {
			val configuration = if( x == null || x.size <= 0 ) {
				MCMC.randomConfiguration(this.nUnits, this.nStatesOfUnit, 
					this.proposedDistributions.proposedDistributionsAtAllUnits )
			  } else {
				x.clone
			  }
			val energy = MCMC.totalE( configuration, this.interactions)
			val initial = State(configuration, energy, kT, rate, 0, finalT)
			initial
		   } )
    
	//val independentMarkovProcesses = 
	//	  initialStates.map( initial => markovChain(initial) )

	@annotation.tailrec
	def equilSampling(mthIteration: Int, 
		//independentMC: Array[Process[State]],
		initialStates: Array[State],
		nInitialIterationsPerUnit: Int ): (Array[MCMC.State], Int, IndependentSamplings ) = {

	 val independentMC = 
		    initialStates.map( initial => MCMC.markovChain(initial, kernel) ).par	// markovChain: Iterator[State]
		//(new IndependentMC(initialStates.size)).zipWithIndex.map( x => { MCMC.markovChain(initialStates(x._2), kernel) } )
		  //initialStates.map( initial => markovChain(initial) ).par

	 @annotation.tailrec
	 def equilSamplingWithIterator(nthIteration: Int, 
		independentMC: IndependentMC,
		nInitialIterationsPerUnit: Int, 
		newNInitialIterationsPerUnit: Int ): (Int, Int, IndependentSamplings, IndependentMC) = {
	
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
			nthIteration + 1, independentMCcont, newNInitialIterationsPerUnit, newNInitialIterationsPerUnit ) 
	  }

	 }

	  val (nNonEquil, nthIteration, independentSamplings, independentMCcont) = equilSamplingWithIterator(
			mthIteration, independentMC, nInitialIterationsPerUnit, newNInitialIterationsPerUnit = 0) 

	  val newNInitialIterationsPerUnit = 0
	  val newNInitialIterations = newNInitialIterationsPerUnit * this.nUnits

          //val nNonEquil = nNonEquilibrium(independentSamplings)

	  if ( nNonEquil <= 0 || nthIteration >= maxExtendedIterations ) {
		(initialStates, nthIteration, independentSamplings)
	  } else {
	  // slightly heat 
	    val (newkT, newRate) = {
		  if ( newNInitialIterations >= 1 ) {
		    val t = 1.2 * finalT
		    val maxRate = math.exp( (math.log(0.00001) - math.log(t / finalT - 1.0 + 0.00001)) / 
					(newNInitialIterations / 5.0) )	/* (newNInitialIterations / 10.0) ) */
		    (t, maxRate)
		  } else {
		    (finalT, 0.0)
		  }
		}

	    val lastStates = independentSamplings.map( s => { 
			val last =  s( s.size -1)
			last.copy( configuration = last.configuration.clone )
			//MCMC.State(last.configuration.clone, last.energy, newkT, newRate, last.step, last.finalkT)
  			} ).toArray

	    //val independentMarkovProcesses = lastStates.map ( s => markovChain(s))

	    //equilSampling(nthIteration + 1, independentMarkovProcesses, nInitialIterations ) 
	    //equilSampling(nthIteration + 1, independentMarkovProcesses, nInitialIterations = 0 ) 

	    equilSampling(nthIteration + 1, lastStates, nInitialIterationsPerUnit = newNInitialIterationsPerUnit ) 
	  }

	}

	//equilSampling(1, independentMarkovProcesses, nInitialIterations ) 
	equilSampling(1, initialStates, nInitialIterationsPerUnit ) 
    }

    def  frequenciesAtUnitInSamples( independentSamplings: IndependentSamplings ): Array[DenseVector[Double]] =
    	MCMC.frequenciesAtUnitInSamples( independentSamplings, this.nStatesOfUnit)

    def pairwiseFrequenciesInSamples( independentSamplings: IndependentSamplings ): Array[DenseMatrix[Double]] =
    	MCMC.pairwiseFrequenciesInSamples( independentSamplings, this.nStatesOfUnit)

  }


  object BM {

	case class BMInteractions(val phia: Array[DenseVector[Double]], val phijab: Array[DenseMatrix[Double]], gauge: String = "ungauged")
	/*
  	case class BMState (bmInteractions: BMInteractions, 
			v: BMInteractions, m: BMInteractions, grad: BMInteractions, 
			learningRate: Double, betaV: Double, betaM: Double, eps: Double, step: Int )
	//betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-8
	*/
  	case class BMState (bmInteractions: BMInteractions, 
			grad: BMInteractions,
			v: BMInteractions, m: BMInteractions,
			betaV: Double, betaM: Double,
			learningRate: Double,
			minLearningRate: Double, maxLearningRate: Double,
			rateDecrease: Double, rateIncrease: Double,
			learningRates: BMInteractions,
			eps: Double, step: Int )

	// 
	def createBMState (bmInteractions: BMInteractions, 
			grad: BMInteractions, 
			v: BMInteractions, m: BMInteractions, 
			betaV: Double, betaM: Double,
			learningRate: Double, 
			eps: Double, step: Int ) = {

		val phi: Array[DenseVector[Double]] = Array()
		val phij: Array[DenseMatrix[Double]] = Array()
		val learningRates = BMInteractions(phi, phij)

		BMState(bmInteractions = bmInteractions, grad = grad, 
			v = v, m = m, 
			betaV = betaV, betaM = betaM, 
			learningRate = learningRate, 
			minLearningRate = learningRate, maxLearningRate = learningRate,
			rateDecrease = 1.0, rateIncrease = 1.0,
			learningRates = learningRates,
			eps = eps, step = step)
	}

	// For ModAdam
	def createBMState (bmInteractions: BMInteractions, 
			grad: BMInteractions, 
			v: BMInteractions, m: BMInteractions, 
			betaV: Double, betaM: Double,
			learningRate: Double, 
			learningRates: BMInteractions,
			eps: Double, step: Int ) = {

		BMState(bmInteractions = bmInteractions, grad = grad, 
			v = v, m = m, 
			betaV = betaV, betaM = betaM, 
			learningRate = learningRate, 
			minLearningRate = learningRate, maxLearningRate = learningRate,
			rateDecrease = 1.0, rateIncrease = 1.0,
			learningRates = learningRates,
			eps = eps, step = step)
	}

	// For RPROP-LR
	def createBMState (bmInteractions: BMInteractions, 
			grad: BMInteractions, 
			learningRate: Double, 
			minLearningRate: Double, maxLearningRate: Double,
			rateDecrease: Double, rateIncrease: Double,
			learningRates: BMInteractions,
			eps: Double, step: Int ) = {
		/*
	  v = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
		 (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) )
			) 
		*/

		val phi: Array[DenseVector[Double]] = Array()
		val phij: Array[DenseMatrix[Double]] = Array()
		val v = BMInteractions(phi, phij)

		BMState(bmInteractions = bmInteractions, grad = grad, 
			v = v, 
			m = learningRates, 
			betaV = rateDecrease, betaM = rateIncrease, 
			learningRate = learningRate, 
			minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
			rateDecrease = rateDecrease, rateIncrease = rateIncrease,
			learningRates = learningRates,
			eps = eps, step = step)
	}

	def hJToPhi (hia: Array[DenseVector[Double]],  Jijab: Array[DenseMatrix[Double]],
		fia: Array[DenseVector[Double]] ) = {
		
		val nUnits = hia.size
  		val nPairs = Jijab.size
		val phijab = Jijab
		val phia =  hia.zipWithIndex.map( x => {
				val hi = x._1
				val i = x._2
				//val xj = new Array[DenseVector[Double]](nUnits)
				//xj.zipWithIndex.map( xx => {
				//  	val j = xx._2
				Array.range(0, nUnits).map( j => {
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

	def phiTohJ (phia: Array[DenseVector[Double]],  phijab: Array[DenseMatrix[Double]],
		fia: Array[DenseVector[Double]] ) = {
		
		val nUnits = phia.size
  		val nPairs = phijab.size
		val Jijab = phijab
		val hia =  phia.zipWithIndex.map( x => {
				val phi = x._1
				val i = x._2
				val xj = new Array[DenseVector[Double]](nUnits)
				xj.zipWithIndex.map( xx => {
				    	val j = xx._2
				    	if ( i == j ) {
						phi
				    	} else {
						val (ij, a, b) = pairIndex(i, j, 0, 1)
						  if ( a == 0 )
							- phijab(ij) * fia(j)
						  else
							- phijab(ij).t * fia(j)
				    	}
				    } ).reduce(_ + _)
				
			} )
		(hia, Jijab)
	}

	def toIsingGauge( hia: Array[DenseVector[Double]], Jijab: Array[DenseMatrix[Double]] ) = {
	  //pairIndex: (Int, Int, Int, Int) => Tuple3[Int, Int, Int] 

	  import breeze.linalg._

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
	  val nPairs = Jijab.size
		
	  val indexIJ = Array.range(0, nPairs).map( ij => inversePairIndex(ij) )   

	  val fullJija_ = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
				{
					if ( i == j ) {
						DenseVector.fill(nStatesOfUnit)(0.0)
					} else {
						val (ij, a, b) = pairIndex(i, j, 0, 1)
						val fullJijab = if( a == 0 ) Jijab(ij) else Jijab(ij).t
						// fullJijab * fia(j)
						breeze.linalg.sum(fullJijab(* , ::) ) / nStatesOfUnit.toDouble
					}
				} )	
			} )

	// fullJij_b(i)(j)_(b) =  (fia(i).t * fullJijab(i)(j))(b) 
	//			= ( fia(i).t * fullJijab(j)(i).t )(b) = (fullJijab(j)(i) * fia(i)).t(b) 
	//			= fullJija_(j)(i).t(b)

	  val fullJij_b = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
				{
					if ( i == j ) {
						DenseVector.fill(nStatesOfUnit)(0.0).t
					} else {
						//fia(i).t * fullJijab(i)(j)
						fullJija_(j)(i).t
					}
				} )	
			} )
	  val fullJij__ = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
				{
					if ( i == j ) {
						0.0
					} else {
						//fia(i).t * fullJija_(i)(j)
						breeze.linalg.sum(fullJija_(i)(j)) / nStatesOfUnit
					}
				} )	
			} )

	  val Jijab_Ising = Array.range(0, nPairs).map( ij =>
			{
				val (i, j) = indexIJ(ij)
				val x = Jijab(ij)(::, *) - fullJija_(i)(j)
				val y = x.t(::, *) - fullJij_b(i)(j).t
				val yt = y.t
				yt + fullJij__(i)(j)
				//(Jijab(ij)(::,*) - fullJija_(i)(j))(*,::) - fullJij_b(i)(j) + fullJij__(i)(j)
			} )

	  val hia_Ising = Array.range(0, nUnits).map( i =>
			{
				//hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - fia(i).t * hia(i)
				hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - breeze.linalg.sum(hia(i)) / nStatesOfUnit
			} )

	  (hia_Ising, Jijab_Ising)
	}

	def toCompositionGauge( hia: Array[DenseVector[Double]], Jijab: Array[DenseMatrix[Double]],
		fia: Array[DenseVector[Double]] ) = {
	  //pairIndex: (Int, Int, Int, Int) => Tuple3[Int, Int, Int] 

	  import breeze.linalg._

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
	  val nPairs = Jijab.size
		
	  val indexIJ = Array.range(0, nPairs).map( ij => inversePairIndex(ij) )   

	  val fullJija_ = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
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
	//			= ( fia(i).t * fullJijab(j)(i).t )(b) = (fullJijab(j)(i) * fia(i)).t(b) 
	//			= fullJija_(j)(i).t(b)

	  val fullJij_b = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
				{
					if ( i == j ) {
						DenseVector.fill(nStatesOfUnit)(0.0).t
					} else {
						//fia(i).t * fullJijab(i)(j)
						fullJija_(j)(i).t
					}
				} )	
			} )
	  val fullJij__ = Array.range(0, nUnits).map( i =>
			{
				Array.range(0, nUnits).map( j =>
				{
					if ( i == j ) {
						0.0
					} else {
						fia(i).t * fullJija_(i)(j)
					}
				} )	
			} )

	  val Jijab_comp = Array.range(0, nPairs).map( ij =>
			{
				val (i, j) = indexIJ(ij)
				val x = Jijab(ij)(::, *) - fullJija_(i)(j)
				val y = x.t(::, *) - fullJij_b(i)(j).t
				val yt = y.t
				yt + fullJij__(i)(j)
				//(Jijab(ij)(::,*) - fullJija_(i)(j))(*,::) - fullJij_b(i)(j) + fullJij__(i)(j)
			} )

	  val hia_comp = Array.range(0, nUnits).map( i =>
			{
				hia(i) + fullJija_(i).reduce(_ + _) - fullJij__(i).sum - fia(i).t * hia(i)
			} )

	  (hia_comp, Jijab_comp)
	}

	def toZeroSumGauge( hia: Array[DenseVector[Double]], Jijab: Array[DenseMatrix[Double]] ) = {

	  import breeze.linalg._

          val nUnits = hia.size
          val nStatesOfUnit = hia(0).size
	  val nPairs = Jijab.size
	
	  val nStatesOfUnitSquare = nStatesOfUnit * nStatesOfUnit
	  val Jijab_ZeroSum = Array.range(0, nPairs).map( ij => {
				val Jij__ = breeze.linalg.sum(Jijab(ij)) / nStatesOfUnitSquare.toDouble
				Jijab(ij) - Jij__
			      } )
	  val hia_ZeroSum = Array.range(0, nUnits).map( i => {
				val hi_ = breeze.linalg.sum(hia(i)) / nStatesOfUnit.toDouble
				hia(i) - hi_
			    } )

	  (hia_ZeroSum, Jijab_ZeroSum)
	}

	def toGauge( hia: Array[DenseVector[Double]], Jijab: Array[DenseMatrix[Double]], gauge: String, 
			fia: Array[DenseVector[Double]] = Array(DenseVector[Double]()) ) = {

	    if ( gauge == "ungauged" || gauge == "unused" ) {
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

	def newPhiThroughhJ(bmInteractions: BMInteractions, fia: Array[DenseVector[Double]], gauge: String):
		(MCMC.Interactions, BMInteractions) = {

	    val phi_ = gauge.slice(0,4)

	    if ( gauge == "ungauged" || gauge == "unused" ) {
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
	def newPhiThroughhJ(bmInteractions: BMInteractions, fia: Array[DenseVector[Double]]):
		(MCMC.Interactions, BMInteractions) = {
		newPhiThroughhJ(bmInteractions, fia, gauge = "ungauged")
	}

	def newPhiThroughhJ(bmState: BMState, fia: Array[DenseVector[Double]], gauge: String):
		(MCMC.Interactions, BMState) = {

	    val bmInteractions = bmState.bmInteractions
	    val (newInteractions, newBMInteractions) = newPhiThroughhJ(bmInteractions, fia, gauge)

	    val newBMState = BMState( bmInteractions = newBMInteractions, 
			grad = bmState.grad,
			v = bmState.v, m = bmState.m,
			betaV = bmState.betaV, betaM = bmState.betaM, 
			learningRate = bmState.learningRate,
			minLearningRate = bmState.minLearningRate, maxLearningRate = bmState.maxLearningRate,
			rateDecrease = bmState.rateDecrease, rateIncrease = bmState.rateIncrease,
			learningRates = bmState.learningRates,
			eps = bmState.eps, step = bmState.step)

	    (newInteractions, newBMState)
	}
	def newPhiThroughhJ(bmState: BMState, fia: Array[DenseVector[Double]]):
		(MCMC.Interactions, BMState) = {
		newPhiThroughhJ(bmState, fia, gauge = "ungauged")
	}

	def softThresholdingForGL1(oldBMInteractions: BMInteractions,
		optMethod: String, regTerm: String, 
		propL1h: Double, propL1J: Double, 
		lambdaPhi: Double, lambdaPhij: Double,
		newBMStateWithoutSoft: BMState, learningRates: BMInteractions ) = {

	// The increments for interactions must be proportional to learningRates that depend on each parameter.

	// Not completed yet; Only for optMethod == "NAG", L1 corrections are included to new v.

	  if(regTerm == "L2" || ( (regTerm == "L1L2" || regTerm == "GL1L2") && ( propL1h == 0.0 && propL1J == 0.0) ) ) {
		newBMStateWithoutSoft
	  } else {
	    val pL1h = if ( regTerm == "GL1" ) 1.0 
			else if ( regTerm == "GL1L2" ) {
			  assert(propL1h >= 0.0 && propL1h <= 1.0)
			  propL1h
			} else {
			  sys.error(s"Not supported: ${regTerm}")
			  propL1h
			}
	    val pL1J = if ( regTerm == "GL1" ) 1.0 
			else if ( regTerm == "GL1L2" ) {
			  assert(propL1J >= 0.0 && propL1J <= 1.0)
			  propL1J
			} else {
			  sys.error(s"Not supported: ${regTerm}")
			  propL1J
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
			val gammai = learningRates.phia(i) * lambdaPhi * pL1h

			val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
			/*
			val signi = if (l2NormPhia(i) == 0.0 ) 
					DenseVector.zeros[Double](nStatesOfUnit)
				     else oldPhi / l2NormPhia(i)
			*/

			val phiWithoutL1 = phi + gammai *:* signi

			/* V2: incorrect
			val phiInSmallCase = phiWithoutL1.mapPairs( (s, p) => { 
						if(phiWithoutL1(s).abs <= gammai(s)) 0.0 else 1.0
					} )
			val suffLarge = breeze.linalg.sum( phijInSmallCase )

			if ( suffLarge == 0 ) {
			  phiInSmallCase
			} else {
			*/
			    phiWithoutL1.mapPairs( (s, p) => {
					if ( p > gammai(s) ) p - gammai(s)
					else if ( p < - gammai(s) ) p + gammai(s)
					else 0.0
					/*
					val dp =
					    if ( l2NormPhia(i) == 0.0 ) { 
						if( p > gammai(s) )
							p - gammai(s)
						else if ( p < - gammai(s) )
							p + gammai(s)
						else
							0.0
					    } else {
						p - gammai(s) * (oldPhi(s) / l2NormPhia(i) )
					    }

					if ( oldPhi(s) > 0.0 ) {
						if ( dp > 0.0 )
							dp
						else
							0.0
					} else if ( oldPhi(s) < 0.0 ) {
						if ( dp < 0.0 )
							dp
						else
							0.0
					} else {
						dp
					}

					*/
				} )
		     // }
		    } )
		}

	    val l2NormPhijab = oldBMInteractions.phijab.map( phij => math.sqrt( breeze.linalg.sum(phij.map(x => x * x) ) ) ) 
	    val phijab = if (pL1J <= 0.0 ) {
		    newBMInteractionsWithoutSoft.phijab
		} else {
		    newBMInteractionsWithoutSoft.phijab.zipWithIndex.map( xx => {
			val (phij, ij) = xx
			val oldPhij = oldBMInteractions.phijab(ij)
			val gammaij = learningRates.phijab(ij) * lambdaPhij * pL1J

			//val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
			val signij = if (l2NormPhijab(ij) == 0.0 ) 
					DenseMatrix.zeros[Double](nStatesOfUnit,nStatesOfUnit)
				     else oldPhij / l2NormPhijab(ij)

		  	val phijWithoutL1 = phij + gammaij *:* signij

			/* V2: incorrect
			val phijInSmallCase = phijWithoutL1.mapPairs( (s, p) => { 
						val (a, b) = s
						if(phijWithoutL1(a,b).abs <= gammaij(a,b)) 0.0 else 1.0
					} )
			val suffLarge = breeze.linalg.sum( phijInSmallCase )

			if ( suffLarge == 0 ) {
			  phijInSmallCase
			} else {
			*/
				
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
			  phijWithoutL1.mapPairs( (s, p) => {
					val (a, b) = s
					/*
					if ( p > gammaij(a, b) ) p - gammaij(a,b)
					else if ( p < - gammaij(a,b) ) p + gammaij(a,b)
					else 0.0
					*/

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
				} )
			//}

			} )
		}
					
	    val newBMInteractionsWithL1 = BMInteractions(phia, phijab)

	  val newV = 
		newBMStateWithoutSoft.v
	  /*
	    if ( optMethod == "NAG" ) {

		/* not completed; probably the following correction is not needed. */
	        val learnRate = newBMStateWithoutSoft.learningRate  

	        val vPhia = if (pL1h <= 0.0) {
		    newBMStateWithoutSoft.v.phia
		} else {
		    newBMStateWithoutSoft.v.phia.zipWithIndex.map( xx => {
			val (v, i) = xx
			val oldPhi = oldBMInteractions.phia(i)

			val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

			val vWithoutL1 = v +  learnRate * lambdaPhi * pL1h * signi

			val gammai = learningRates.phia(i) * lambdaPhi * pL1h
			val phi = newBMInteractionsWithoutSoft.phia(i)
			val phiWithoutL1 = phi + gammai *:* signi

			val vL1 = phiWithoutL1.mapPairs( (s, p) => {
					if ( p > gammai(s) ) - gammai(s)
					else if ( p < - gammai(s) ) gammai(s)
					else -p
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

			val vWithoutL1 = v +  learnRate * lambdaPhij * pL1J * signij

			val gammaij = learningRates.phijab(ij) * lambdaPhij * pL1J
			val phij = newBMInteractionsWithoutSoft.phijab(ij)
			val phijWithoutL1 = phij + gammaij *:* signij

			val vL1 = phijWithoutL1.mapPairs( (s, p) => {
					val (a, b) = s
					if ( p > gammaij(a, b) ) - gammaij(a,b)
					else if ( p < - gammaij(a,b) ) gammaij(a,b)
					else  - p
				} )
			vWithoutL1 + vL1
			} )
		}

	    	val vWithSoft = BMInteractions(vPhia, vPhijab) 
	    	vWithSoft
	    } else { 
		/* The correction for v is ignored for other methods */
		newBMStateWithoutSoft.v
	    }
	  */

	    BMState( bmInteractions = newBMInteractionsWithL1, 
		grad = newBMStateWithoutSoft.grad,
		v = newV,
		m = newBMStateWithoutSoft.m, 
		betaV = newBMStateWithoutSoft.betaV, betaM = newBMStateWithoutSoft.betaM,
		learningRate = newBMStateWithoutSoft.learningRate,
		minLearningRate = newBMStateWithoutSoft.minLearningRate, maxLearningRate = newBMStateWithoutSoft.maxLearningRate,
		rateDecrease = newBMStateWithoutSoft.rateDecrease, rateIncrease = newBMStateWithoutSoft.rateIncrease,
		learningRates = newBMStateWithoutSoft.learningRates,
		eps = newBMStateWithoutSoft.eps, step = newBMStateWithoutSoft.step)

	  }

	}

	def softThresholdingForL1(oldBMInteractions: BMInteractions,
		optMethod: String, regTerm: String,
		propL1h: Double, propL1J: Double,
		lambdaPhi: Double, lambdaPhij: Double,
		newBMStateWithoutSoft: BMState, learningRates: BMInteractions ) = {

	// The increments for interactions must be proportional to learningRates that depend on each parameter.

	// Not completed yet; Only for optMethod == "NAG", L1 corrections are included to new v.

	  if(regTerm == "L2" || ( (regTerm == "L1L2" || regTerm == "GL1L2") && (propL1h == 0.0 && propL1J == 0.0) ) ) {
		newBMStateWithoutSoft
	  } else {
	    val pL1h = if ( regTerm == "L1" ) 1.0 
			else if ( regTerm == "L1L2" ) {
			  assert(propL1h >= 0.0 && propL1h <= 1.0)
			  propL1h
			} else {
			  sys.error(s"Not supported: ${regTerm}")
			  propL1h
			}
	    val pL1J = if ( regTerm == "L1" ) 1.0 
			else if ( regTerm == "L1L2" ) {
			  assert(propL1J >= 0.0 && propL1J <= 1.0)
			  propL1J
			} else {
			  sys.error(s"Not supported: ${regTerm}")
			  propL1J
			}

	    val newBMInteractionsWithoutSoft = newBMStateWithoutSoft.bmInteractions

	    val nStatesOfUnit = learningRates.phia(0).size

	    val phia = if (pL1h <= 0.0 ) {
		    newBMInteractionsWithoutSoft.phia
		} else {
		    newBMInteractionsWithoutSoft.phia.zipWithIndex.map( xx => {
			val (phi, i) = xx
			val oldPhi = oldBMInteractions.phia(i)
			val gammai = learningRates.phia(i) * lambdaPhi * pL1h

			val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

			val phiWithoutL1 = phi + gammai *:* signi
			phiWithoutL1.mapPairs( (s, p) => {
					if ( p > gammai(s) ) p - gammai(s)
					else if ( p < - gammai(s) ) p + gammai(s)
					else 0.0
				} )
			} )
		}

	    val phijab = if (pL1J <= 0.0 ) {
		    newBMInteractionsWithoutSoft.phijab
		} else {
		    newBMInteractionsWithoutSoft.phijab.zipWithIndex.map( xx => {
			val (phij, ij) = xx
			val oldPhij = oldBMInteractions.phijab(ij)
			val gammaij = learningRates.phijab(ij) * lambdaPhij * pL1J

			val signij = oldPhij.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 
			val phijWithoutL1 = phij + gammaij *:* signij
			phijWithoutL1.mapPairs( (s, p) => {
					val (a, b) = s
					if ( p > gammaij(a, b) ) p - gammaij(a,b)
					else if ( p < - gammaij(a,b) ) p + gammaij(a,b)
					else 0.0
				} )
			} )
		}
					
	    val newBMInteractionsWithL1 = BMInteractions(phia, phijab)

	  val newV = 
	    if ( optMethod == "NAG" ) {

	        val learnRate = newBMStateWithoutSoft.learningRate  

	        val vPhia = if (pL1h <= 0.0 ) {
		    newBMStateWithoutSoft.v.phia
		} else {
		    newBMStateWithoutSoft.v.phia.zipWithIndex.map( xx => {
			val (v, i) = xx
			val oldPhi = oldBMInteractions.phia(i)

			val signi = oldPhi.map( x => {if ( x > 0.0 ) 1.0 else if (x < 0.0 ) -1.0 else 0.0 } ) 

			val vWithoutL1 = v +  learnRate * lambdaPhi * pL1h * signi

			val gammai = learningRates.phia(i) * lambdaPhi * pL1h
			val phi = newBMInteractionsWithoutSoft.phia(i)
			val phiWithoutL1 = phi + gammai *:* signi

			val vL1 = phiWithoutL1.mapPairs( (s, p) => {
					if ( p > gammai(s) ) - gammai(s)
					else if ( p < - gammai(s) ) gammai(s)
					else -p
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

			val vWithoutL1 = v +  learnRate * lambdaPhij * pL1J * signij

			val gammaij = learningRates.phijab(ij) * lambdaPhij * pL1J
			val phij = newBMInteractionsWithoutSoft.phijab(ij)
			val phijWithoutL1 = phij + gammaij *:* signij

			val vL1 = phijWithoutL1.mapPairs( (s, p) => {
					val (a, b) = s
					if ( p > gammaij(a, b) ) - gammaij(a,b)
					else if ( p < - gammaij(a,b) ) gammaij(a,b)
					else  - p
				} )
			vWithoutL1 + vL1
			} )
		}

	    	val vWithSoft = BMInteractions(vPhia, vPhijab) 
	    	vWithSoft
	    } else { 
		/* The correction for v is ignored for other methods */
		newBMStateWithoutSoft.v
	    }

	    BMState( bmInteractions = newBMInteractionsWithL1, 
		grad = newBMStateWithoutSoft.grad,
		v = newV,
		m = newBMStateWithoutSoft.m, 
		betaV = newBMStateWithoutSoft.betaV, betaM = newBMStateWithoutSoft.betaM,
		learningRate = newBMStateWithoutSoft.learningRate,
		minLearningRate = newBMStateWithoutSoft.minLearningRate, maxLearningRate = newBMStateWithoutSoft.maxLearningRate,
		rateDecrease = newBMStateWithoutSoft.rateDecrease, rateIncrease = newBMStateWithoutSoft.rateIncrease,
		learningRates = newBMStateWithoutSoft.learningRates,
		eps = newBMStateWithoutSoft.eps, step = newBMStateWithoutSoft.step)
	  }

	}

	// GL1 means group L1 for J; R = sum_i sum_a |hia|  +    sum_ij sqrt sum_ab Jijab^2 
	def dLdphiWithGL1L2(propL1h: Double, propL1J: Double,
		fia: Array[DenseVector[Double]], fijab: Array[DenseMatrix[Double]],
		pia: Array[DenseVector[Double]], pijab: Array[DenseMatrix[Double]],
		phia: Array[DenseVector[Double]],  phijab: Array[DenseMatrix[Double]], 
		//hia: Array[DenseVector[Double]],  Jijab: Array[DenseMatrix[Double]], 
		lambdaPhi: Double, lambdaPhij: Double  ) = {

		//val (phia, phijab) = hJToPhi (hia,  Jijab, fia)

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
			          fi - pi - lambdaPhi * ( propL2h * phi )
				} else {
				  /* for L1 */
				  val signPhi = phi.map( x => {
					if (x > 0.0) 1.0 
					else if( x < 0.0) -1.0 
					else { 0.0 }
					} )
				  /*
				  val signPhi = phi.map( x => {
					if ( l2NormPhia(i) == 0.0 ) {
						0.0
					} else {
						x / l2NormPhia(i)
					}
				  } )
				  */
						
				  if (propL2h == 0.0 )
			            fi - pi - lambdaPhi * ( propL1h * signPhi )
				  else
			            fi - pi - lambdaPhi * ( propL1h * signPhi + propL2h * phi )
				}
			} )

		val l2NormPhijab = phijab.map( phij => math.sqrt( breeze.linalg.sum(phij.map(x => x * x) ) ) ) 

		val dldphijab = new Array[DenseMatrix[Double]](nPairs)
		cfor(0)( i => i < nUnits, i => i + 1)( i => {
		  cfor(0)( j => j < i, j => j + 1)( j => {
			val (ij, a, b) = pairIndex(i, j, 0, 1)
			dldphijab(ij) =
			  if ( propL1J == 0.0 ) {
			  	if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * ( propL2J * phijab(ij) )
				} else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * ( propL2J * phijab(ij) )
				}
			  } else {
			  	/* for L1
			  	val signPhij = phijab(ij).map( x => {
					if (x > 0.0) 1.0 
					else if( x < 0.0) -1.0 
					else { 0.0 }
					} )
			  	*/
			  	val signPhij = phijab(ij).map( x => {
					if ( l2NormPhijab(ij) == 0.0 ) {
						0.0
					} else {
						x / l2NormPhijab(ij)
					}

				} )

				if ( propL2J == 0.0 ) 
				  if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij )
				  } else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij )
				  }
				else
				  if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
				  } else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
				  }
			  }
		  } )
		} )
		(dldphia, dldphijab)
	}

	def dLdphiWithL1L2(propL1h: Double, propL1J: Double,
		fia: Array[DenseVector[Double]], fijab: Array[DenseMatrix[Double]],
		pia: Array[DenseVector[Double]], pijab: Array[DenseMatrix[Double]],
		phia: Array[DenseVector[Double]],  phijab: Array[DenseMatrix[Double]], 
		//hia: Array[DenseVector[Double]],  Jijab: Array[DenseMatrix[Double]], 
		lambdaPhi: Double, lambdaPhij: Double  ) = {

		//val (phia, phijab) = hJToPhi (hia,  Jijab, fia)

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
			          fi - pi - lambdaPhi * ( propL2h * phi )
				} else {
				  val signPhi = phi.map( x => {
					if (x > 0.0) 1.0 
					else if( x < 0.0) -1.0 
					else { 0.0 }
					} )
				  if ( propL2h == 0.0 )
			            fi - pi - lambdaPhi * ( propL1h * signPhi )
				  else
			            fi - pi - lambdaPhi * ( propL1h * signPhi + propL2h * phi )
				}
			} )

		val dldphijab = new Array[DenseMatrix[Double]](nPairs)
		cfor(0)( i => i < nUnits, i => i + 1)( i => {
		  cfor(0)( j => j < i, j => j + 1)( j => {
			val (ij, a, b) = pairIndex(i, j, 0, 1)
			dldphijab(ij) = 
			  if ( propL1J == 0.0 ) {
			  	if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * ( propL2J * phijab(ij) )
				} else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * ( propL2J * phijab(ij) )
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
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij )
				  } else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij )
				  }
				else
				  if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
				  } else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * (propL1J * signPhij + propL2J * phijab(ij) )
				  }
			}
		  } )
		} )
		(dldphia, dldphijab)
	}

	def dLdphiWithL1(fia: Array[DenseVector[Double]], fijab: Array[DenseMatrix[Double]],
		pia: Array[DenseVector[Double]], pijab: Array[DenseMatrix[Double]],
		phia: Array[DenseVector[Double]],  phijab: Array[DenseMatrix[Double]], 
		//hia: Array[DenseVector[Double]],  Jijab: Array[DenseMatrix[Double]], 
		lambdaPhi: Double, lambdaPhij: Double ) = {

		//val (phia, phijab) = hJToPhi (hia,  Jijab, fia)
		
		// Please notice that softThresholding will be taken care of later. 

		val nUnits = fia.size
  		val nPairs = fijab.size

		val dldphia = phia.zip(pia).zip(fia).map( x => {
					val phi = x._1._1
					val pi = x._1._2
					val fi = x._2
				        val signPhi = phi.map( x => if (x > 0.0) 1.0 else if( x < 0.0) -1.0 else 0.0)
				        fi - pi - lambdaPhi * signPhi 
				} )

		val dldphijab = new Array[DenseMatrix[Double]](nPairs)
		cfor(0)( i => i < nUnits, i => i + 1)( i => {
		  cfor(0)( j => j < i, j => j + 1)( j => {
			val (ij, a, b) = pairIndex(i, j, 0, 1)
			val signPhij = phijab(ij).map( x => if (x > 0.0) 1.0 else if( x < 0.0) -1.0 else 0.0)
			dldphijab(ij) = if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
						fia(i) * (pia(j) - fia(j)).t - lambdaPhij * signPhij
				} else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t +
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * signPhij
				}
		  } )
		} )
		(dldphia, dldphijab)
	}

	def dLdphiWithL2(fia: Array[DenseVector[Double]], fijab: Array[DenseMatrix[Double]],
		pia: Array[DenseVector[Double]], pijab: Array[DenseMatrix[Double]],
		phia: Array[DenseVector[Double]],  phijab: Array[DenseMatrix[Double]], 
		//hia: Array[DenseVector[Double]],  Jijab: Array[DenseMatrix[Double]], 
		lambdaPhi: Double, lambdaPhij: Double  ) = {

		//val (phia, phijab) = hJToPhi (hia,  Jijab, fia)
		
		val nUnits = fia.size
  		val nPairs = fijab.size

		val dldphia = phia.zip(pia).zip(fia).map( x => {
					val phi = x._1._1
					val pi = x._1._2
					val fi = x._2
					fi - pi - lambdaPhi * phi  
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

		val dldphijab = new Array[DenseMatrix[Double]](nPairs)
		cfor(0)( i => i < nUnits, i => i + 1)( i => {
		  cfor(0)( j => j < i, j => j + 1)( j => {
			val (ij, a, b) = pairIndex(i, j, 0, 1)
			dldphijab(ij) = if ( a == 0 ) {
					fijab(ij) - pijab(ij) + (pia(i) - fia(i)) * fia(j).t + 
					fia(i) * (pia(j) - fia(j)).t - lambdaPhij * phijab(ij)
				} else {
					fijab(ij) - pijab(ij) + fia(j) * (pia(i) - fia(i)).t + 
						(pia(j) - fia(j)) * fia(i).t - lambdaPhij * phijab(ij)
				}
		  } )
		} )
		(dldphia, dldphijab)
	}

	//val (dldphia, dldphijab) = dLdphi(fia, fijab, pia, pijab, hia, Jijab, lambda )

	def learningByNAG(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ): 
							Tuple2[BMState, Option[BMInteractions]] = {

	  assert(minOrMax.abs == 1)
	  assert(0.0 <= bmState.betaV && bmState.betaV < 1.0)
	// Be careful of the definition of "leaningRate". 
	  val learningRate = bmState.learningRate * ( 1.0 - bmState.betaV )

	  val step = bmState.step + 1
	/* The following code was used until 191016
	  val betaV = if ( step <= 10 ) { 		// betaV == mu
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

	//println( bmState.v.phia.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
	  val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + learningRate * g
		} ) 
	//println( vphia.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
	//println( bmState.v.phijab.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )
	  val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + learningRate * g
		} ) 
	//println( vphijab.map(v => v.map(x => x.abs)).map(v => breeze.linalg.max(v)).max )

	  val phia = Array.range(0, vphia.size).map( i => {
			val s = bmState.bmInteractions.phia(i)
			s + minOrMax.toDouble *
				( betaV * vphia(i) * normalizeBetaVStep + learningRate * grad.phia(i) )
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
	  val phijab = Array.range(0, vphijab.size).map( ij => {
			val s = bmState.bmInteractions.phijab(ij)
			s + minOrMax.toDouble *
				( betaV * vphijab(ij) * normalizeBetaVStep + learningRate * grad.phijab(ij) )
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

	  ( createBMState( bmInteractions = BMInteractions(phia, phijab),
			grad = grad,
			v = BMInteractions(vphia, vphijab), bmState.m,
			betaV = bmState.betaV, betaM = bmState.betaM,
			learningRate = bmState.learningRate,
			learningRates = learningRates,
			eps = bmState.eps, step = step ),
		Option(learningRates) )

	}

	def learningByAdadelta(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ):
							Tuple2[BMState, Option[BMInteractions]] = {

	  assert(minOrMax.abs == 1)
	  val eps = bmState.eps
	  val beta = bmState.betaV
	  assert( bmState.betaM == bmState.betaV )
	  val betaV = beta
	  val betaM = beta

	  val step = bmState.step + 1

	/* modified from the original */

	  val nUnits = bmState.bmInteractions.phia.size
	  val nPairs = (nUnits * ( nUnits - 1)) / 2 
	  val nStatesOfUnit = bmState.bmInteractions.phia(0).size
	  val learningRatePhia = 
		if( step <= 1) {
			(new Array[DenseVector[Double]](nUnits)).map( x => {
				DenseVector.fill(nStatesOfUnit)(bmState.learningRate)
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
			(new Array[DenseMatrix[Double]](nPairs)).map( x => {
				DenseMatrix.fill(nStatesOfUnit, nStatesOfUnit)(bmState.learningRate)
			} )
		} else {
			bmState.m.phijab.zip(bmState.v.phijab).map( x => { 
				val m = x._1
				val v = x._2
				breeze.numerics.sqrt(m + eps) /:/
				  breeze.numerics.sqrt(v + eps)
				} ) 
		}

	  val g2phia = bmState.v.phia.zip(grad.phia).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
		} ) 
	  val g2phijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
		} ) 

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
			betaM * m + (1.0 - betaM) * (d *:* d)
		} ) 
	  val dphijab2 = bmState.m.phijab.zip(dphijab).map( dx => {
			val m = dx._1
			val d = dx._2
			betaM * m + (1.0 - betaM) * (d *:* d)
		} ) 

	  val newBMState = createBMState( bmInteractions = BMInteractions(phia, phijab),
				grad = grad,
				v = BMInteractions(g2phia, g2phijab), m = BMInteractions(dphia2, dphijab2),
				betaV = betaV, betaM = betaM,
				learningRate = bmState.learningRate, 
				learningRates = learningRates,
				eps = eps, step = step )

	  val optionLearningRates: Option[BMInteractions] = None
	  (newBMState, optionLearningRates)	// Adadelta canot be used for the soft thresholding function for L1.
	}

	def learningByAdam(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ):
							Tuple2[BMState, Option[BMInteractions]] = {

	  assert(minOrMax.abs == 1)
	  val learningRate = bmState.learningRate
	  val betaM = bmState.betaM
	  val betaV = bmState.betaV
	  val eps = bmState.eps
	  val step = bmState.step + 1
	  assert(0.0 <= betaM && betaM < 1.0)
	  assert(0.0 <= betaV && betaV < 1.0)

	  val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
			val m = dx._1
			val g = dx._2
			betaM * m + (1.0 - betaM) * g
		} ) 
	  val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
			val m = dx._1
			val g = dx._2
			betaM * m + (1.0 - betaM) * g
		} ) 
	  val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
		} ) 
	  val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
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
			s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
		} ) 
	  val phijab = mphijab.zipWithIndex.map( x => {
			val (m, ij) = x
			val invVt = invSqrtVtPhijab(ij)
			val s = bmState.bmInteractions.phijab(ij)
			s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
		} ) 

	  val learningRatePhia = invSqrtVtPhia.map( invVt => {
				(learningRate * (1.0 - betaM) / betaMt) * invVt
			} )
	  val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
				(learningRate * (1.0 - betaM) / betaMt) * invVt
			} )

	  val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
	  val optionLearningRates: Option[BMInteractions] = None

	  ( createBMState( bmInteractions = BMInteractions(phia, phijab),
		grad = grad,
		v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
		betaV = betaV, betaM = betaM,
		learningRate = learningRate,
		learningRates = learningRates,
		eps = eps, step = step),
	        optionLearningRates	)	// Adam cannot be used for the soft thresholding function for L1.
	}

	def learningByModAdam(bmState: BMState, grad: BMInteractions, minOrMax: Int = 1 ):
							Tuple2[BMState, Option[BMInteractions]] = {

	  assert(minOrMax.abs == 1)
	  val learningRate = bmState.learningRate
	  val betaM = bmState.betaM
	  val betaV = bmState.betaV
	  val eps = bmState.eps
	  val step = bmState.step + 1
	  assert(0.0 <= betaM && betaM < 1.0)
	  assert(0.0 <= betaV && betaV < 1.0)

	  val mphia = bmState.m.phia.zip(grad.phia).map( dx => {
			val m = dx._1
			val g = dx._2
			betaM * m + (1.0 - betaM) * g
		} ) 
	  val mphijab = bmState.m.phijab.zip(grad.phijab).map( dx => {
			val m = dx._1
			val g = dx._2
			betaM * m + (1.0 - betaM) * g
		} ) 
	  val vphia = bmState.v.phia.zip(grad.phia).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
		} ) 
	  val vphijab = bmState.v.phijab.zip(grad.phijab).map( dx => {
			val v = dx._1
			val g = dx._2
			betaV * v + (1.0 - betaV) * (g *:* g)
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
			s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
		} ) 
	  val phijab = mphijab.zipWithIndex.map( x => {
			val (m, ij) = x
			val invVt = invSqrtVtPhijab(ij)
			val s = bmState.bmInteractions.phijab(ij)
			s + (minOrMax * learningRate / betaMt) * ( m *:* invVt )
		} ) 

	  val learningRatePhia = invSqrtVtPhia.map( invVt => {
				(learningRate * (1.0 - betaM) / betaMt) * invVt
			} )
	  val learningRatePhijab = invSqrtVtPhijab.map( invVt => {
				(learningRate * (1.0 - betaM) / betaMt) * invVt
			} )

	  val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
	*/
	//
	  val vtphia = vphia.map( v => breeze.linalg.max(v) ).max 
	  val vtphijab = vphijab.map( v =>  breeze.linalg.max(v) ).max
	  val vtphi = math.max(vtphia, vtphijab) / betaVt

	  val invSqrtVtPhi = 1.0 / (math.sqrt(vtphi) + eps)

	  val phia = mphia.zipWithIndex.map( x => {
			val (m, i) = x
			val s = bmState.bmInteractions.phia(i)
			s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m 
		} ) 
	  val phijab = mphijab.zipWithIndex.map( x => {
			val (m, ij) = x
			val s = bmState.bmInteractions.phijab(ij)
			s + (minOrMax * learningRate / betaMt * invSqrtVtPhi) * m
		} ) 

	  val actualLearnRate = learningRate * ((1.0 - betaM) / betaMt) * invSqrtVtPhi

	  val learningRatePhia = mphia.map( x => x.map( y =>  actualLearnRate ) )
	  val learningRatePhijab = mphijab.map( x => x.map( y =>  actualLearnRate ) )

	  val learningRates = BMInteractions(learningRatePhia, learningRatePhijab)
	//

	  ( createBMState( bmInteractions = BMInteractions(phia, phijab),
			grad = grad,
			v = BMInteractions(vphia, vphijab), m = BMInteractions(mphia, mphijab),
			betaV = betaV, betaM = betaM,
			learningRate = learningRate,
			learningRates = learningRates,
			eps = eps, step = step),
		Option(learningRates) )

	}

	// J. P. Barton1, E. De Leonardis, A. Coucke, and S. Cocco
	// 16_Bioinformatics_32_3089-3097 and supplementary_data
	// By Matteo Figliuzzi, Pierre Barrat-Charlaix, and Martin Weigt
	// Mol. Biol. Evol. 35, 1018-1027, 2018
	def learningByRPROPLR(bmState: BMState, grad: BMInteractions, 
			regTerm: String, propL1h: Double, propL1J: Double, minOrMax: Int = 1 ):
						Tuple2[BMState, Option[BMInteractions]] = {

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
	  val minLearnRatePhia = bmState.minLearningRate
	  val minLearnRatePhijab = bmState.minLearningRate
	  val maxLearnRatePhia = bmState.maxLearningRate
	  val maxLearnRatePhijab = bmState.maxLearningRate

	  val rateRange = (minR: Double, maxR: Double, r: Double) =>  
					if ( r > maxR )
						maxR
					else if ( r < minR )
						minR 
					else
						r
	  val rateRangePhia = rateRange(minLearnRatePhia, maxLearnRatePhia, _: Double)
	  val rateRangePhijab = rateRange(minLearnRatePhijab, maxLearnRatePhijab, _: Double)

	  assert(minOrMax.abs == 1)
	  val rateDecrease = bmState.rateDecrease	// bmState.betaV
	  val rateIncrease = bmState.rateIncrease	// bmState.betaM
	  //assert( rateIncrease > rateDecrease )
	  assert( rateIncrease > 1.0 && rateDecrease < 1.0 )
	  val prevLearningRates = bmState.learningRates	// bmState.m
	  val prevGrad = bmState.grad

	  val step = bmState.step + 1

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
				/*	does not converge.
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
				/*	does not converge.
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

	  val newBMState = createBMState( bmInteractions = BMInteractions(phia, phijab),
			grad = grad,
			learningRate = bmState.learningRate,
			minLearningRate = bmState.minLearningRate, maxLearningRate = bmState.maxLearningRate,
			rateDecrease = bmState.rateDecrease, rateIncrease = bmState.rateIncrease,
			learningRates = learningRates,
			eps = bmState.eps, step = step)

	  (newBMState, Option(learningRates))	// RPROP-LR cannot be used for the soft thresholding function for L1.
	}

    def initializeBMStateForNAG( bmInteractions: BMInteractions, learningRate: Double,
	betaV: Double = 0.9, betaM: Double = 0.0, eps: Double = 1.0e-8, step: Int = 0 ) = {

        initializeBMStateForAdam( bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForAdam( bmInteractions: BMInteractions, learningRate: Double,
	betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-8, step: Int = 0 ) = {

	val nUnits = bmInteractions.phia.size 
	val nPairs = (nUnits * (nUnits - 1)) / 2
	val nStatesofUnit = bmInteractions.phia(0).size

	val grad = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
		 (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )
	val v = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
		 (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) ) 
	val m = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
		 (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )

	val learningRates = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.fill(nStatesofUnit)(learningRate) ), 
	    (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.fill(nStatesofUnit, nStatesofUnit)(learningRate) ) )

	val bmState = BM.createBMState( bmInteractions = bmInteractions, 
		grad = grad, v = v, m = m, betaV = betaV, betaM = betaM, 
		learningRate = learningRate, learningRates = learningRates, eps = eps, step = step)

	bmState
    }

    def initializeBMStateForModAdam( bmInteractions: BMInteractions, learningRate: Double,
	betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-8, step: Int = 0 ) = {

        initializeBMStateForAdam( bmInteractions, learningRate, betaV, betaM, eps, step )
    }

    def initializeBMStateForAdadelta( bmInteractions: BMInteractions, learningRate: Double,
	betaV: Double = 0.9, betaM: Double = 0.9, eps: Double = 1.0e-8, step: Int = 0 ) = {

	assert( betaM == betaV )
        initializeBMStateForAdam( bmInteractions, learningRate,
		betaV, betaM, eps, step )
    }

    def initializeBMStateForRPROPLR( bmInteractions: BMInteractions, learningRate: Double,
	minLearningRate: Double = 0.000001, maxLearningRate: Double = 10.0,
	rateDecrease: Double = 0.5, rateIncrease: Double = 1.2, eps: Double = 1.0e-8, step: Int = 0 ) = {
	// betaV: rateDecrease , betaM: rateInccrease
        // bmState.m: prevLearningRates

	val nUnits = bmInteractions.phia.size 
	val nPairs = (nUnits * (nUnits - 1)) / 2
	val nStatesofUnit = bmInteractions.phia(0).size

	val grad = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.zeros[Double](nStatesofUnit) ), 
	    (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.zeros[Double](nStatesofUnit, nStatesofUnit) ) )

	val learningRates = BM.BMInteractions( (new Array[DenseVector[Double]](nUnits)).map ( ve => DenseVector.fill(nStatesofUnit)(learningRate) ), 
	    (new Array[DenseMatrix[Double]](nPairs)).map ( ma => DenseMatrix.fill(nStatesofUnit, nStatesofUnit)(learningRate) ) )

	val bmState = BM.createBMState( bmInteractions = bmInteractions, 
		grad = grad, 
		learningRate = learningRate,
		minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
		rateDecrease = rateDecrease, rateIncrease = rateIncrease,
		learningRates = learningRates, eps = eps, step = step)

	bmState
    }

    def initializehJ ( fia: Array[DenseVector[Double] ], small_eps: Double = -1.0) = {

	val nStatesOfUnit = fia(0).size
        val eps = if (small_eps <= 0.0) (1.0 / nStatesOfUnit / 100.0 ) else small_eps
	val nUnits = fia.size 
	val nPairs = (nUnits * (nUnits - 1)) / 2
	val Jijab = (new Array[DenseMatrix[Double]](nPairs)).map( 
			x => DenseMatrix.zeros[Double](nStatesOfUnit, nStatesOfUnit) )
	val hia = fia.map ( fi => { 
			val hi = breeze.numerics.log(fi + eps)
			hi - breeze.linalg.sum(hi) / hi.size	// Ising gaige
		} )
		 
	(hia, Jijab)
    }

    def initialize(fia: Array[DenseVector[Double]],
	initialInteractions: MCMC.Interactions,
	optMethod: String,
	learningRate: Double,
	minLearningRate: Double, maxLearningRate: Double,
	rateDecrease: Double, rateIncrease: Double,
	betaV: Double, betaM: Double, eps: Double = 1.0e-8, step: Int = 0 ) = {
	//betaV: Double = 0.999, betaM: Double = 0.9, eps: Double = 1.0e-8, step: Int = 0 ) = {

	val nUnits = fia.size 
	val nPairs = (nUnits * (nUnits - 1)) / 2
	val nStatesofUnit = fia(0).size
	
	val interactions = initialInteractions

	val (phia, phijab) = BM.hJToPhi (interactions.hia, interactions.Jijab, fia )
	val bmInteractions = BM.BMInteractions( phia, phijab, interactions.gauge )

    	val bmState = 
	    if ( optMethod == "NAG" )
		initializeBMStateForNAG(bmInteractions, learningRate, betaV, betaM, eps, step )
	    else if ( optMethod == "ModAdam" )
		initializeBMStateForModAdam(bmInteractions, learningRate, betaV, betaM, eps, step )
	    else if ( optMethod == "Adam" )
		initializeBMStateForAdam(bmInteractions, learningRate, betaV, betaM, eps, step )
	    else if ( optMethod == "Adadelta" )
		initializeBMStateForAdadelta(bmInteractions, learningRate, betaV, betaM, eps, step )
	    else if  ( optMethod == "RPROP-LR" || optMethod == "MF" )
		initializeBMStateForRPROPLR(bmInteractions = bmInteractions, learningRate = learningRate, 
		minLearningRate = minLearningRate, maxLearningRate = maxLearningRate,
		rateDecrease = rateDecrease, rateIncrease = rateIncrease, 
		eps = eps, step = step)
	    else {
		sys.error("Not supported: " + optMethod) 
		initializeBMStateForAdam(bmInteractions, learningRate, betaV, betaM, eps, step )
	    }

	(interactions, bmState)
    }

    def runBM1Step( outputFiles: HashMap[String, File],
		stateOrderString: String,
		interactions: MCMC.Interactions,
		initialConfigurations: Array[Array[Byte]],
		nInitialIterationsPerUnit: Int = 100, 
		everyNIterationsPerUnit: Int = 5,
		nSamples: Int = 1000,
	      //nIndependentMC: Int = 10,
		initialT: Double = 1.2,
		annealingRate: Double = 0.99,
		maxExtendedIterations: Int = 10,
	      //kernel: MCMC.State => Rand[MCMC.State],
		mcmcKernel: String = "MultiBlockMH", // or "MH" // or "Gibbs"

		fia: Array[DenseVector[Double]],
		fijab: Array[DenseMatrix[Double]],

		regTerm: String,
		propL1h: Double, propL1J: Double,

		lambdaPhi: Double, lambdaPhij: Double,
		bmState: BMState,
		optMethod: String
		) = {

	require( mcmcKernel == "MH" || mcmcKernel == "MultiBlockMH" || mcmcKernel == "Gibbs", 
			"Not supported: mcmcKernel == MH, MultiBlockMH or Gibbs" )

  	val mcmc = new MCMC(outputFiles, stateOrderString, interactions, proposedPia = fia )

	val (initialStates, nExtendedIterations, independentSamplings) = 
	  if ( mcmcKernel == "MH" ) {
	    mcmc.runMC( initialConfigurations,
		nInitialIterationsPerUnit,
		everyNIterationsPerUnit,
		nSamples,
		initialT,
		annealingRate,
		finalT = 1.0,
		maxExtendedIterations = maxExtendedIterations,
		kernel = mcmc.kernelMH)
       	  } else if ( mcmcKernel == "MultiBlockMH" ) {
	    mcmc.runMC( initialConfigurations,
		nInitialIterationsPerUnit,
		everyNIterationsPerUnit,
		nSamples,
		initialT,
		annealingRate,
		finalT = 1.0,
		maxExtendedIterations = maxExtendedIterations,
		kernel = mcmc.kernelMultiBlockMH)
          } else { 	// if ( mcmcKernel == "Gibbs" ) {
	    mcmc.runMC( initialConfigurations,
		nInitialIterationsPerUnit,
		everyNIterationsPerUnit,
		nSamples,
		initialT,
		annealingRate,
		finalT = 1.0,
		maxExtendedIterations = maxExtendedIterations,
		kernel = mcmc.kernelGibbs)
          }


	val pia = mcmc.frequenciesAtUnitInSamples( independentSamplings )
	val pijab = mcmc.pairwiseFrequenciesInSamples( independentSamplings)
	val ensembleAverages = MCMC.EnsembleAverages(pia, pijab)

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

	val oldBMInteractions = bmState.bmInteractions 

	//val (nextBMStateWithoutSoft: BMState, learningRates: BMInteractions) =
	val nextBMStateWithoutSoftPlus =
	    	if ( optMethod == "NAG" )
			BM.learningByNAG(bmState, grad )
		else if ( optMethod == "ModAdam" ) 
			BM.learningByModAdam(bmState, grad )
		else if ( optMethod == "Adam" ) 
			BM.learningByAdam(bmState, grad )
		else if ( optMethod == "Adadelta" ) 
			BM.learningByAdadelta(bmState, grad )
		else if ( optMethod == "RPROP-LR" || optMethod == "MF" ) 
			BM.learningByRPROPLR(bmState, grad, regTerm, propL1h, propL1J )
		else {
			sys.error("Not supported: " + optMethod) 
			BM.learningByAdam(bmState, grad )
		}
	val nextBMStateWithoutSoft = nextBMStateWithoutSoftPlus._1
	val optionLearningRates = nextBMStateWithoutSoftPlus._2

	val nextBMState =
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

	(mcmc, initialStates, nExtendedIterations, independentSamplings, ensembleAverages, nextBMState)
    }

    def calcKL(pia1: Array[DenseVector[Double]], pijab1: Array[DenseMatrix[Double]],
  	pia2: Array[DenseVector[Double]], pijab2: Array[DenseMatrix[Double]], eps: Double = 0.00001 ) = {

	val KL1 = pia1.zip(pia2).
		map ( pi => breeze.linalg.sum( pi._1 *:* breeze.numerics.log((pi._1 + eps) /:/ (pi._2 + eps )) ) ).sum
	
	val KL2 = pijab1.zip(pijab2).
		map ( pij => breeze.linalg.sum( pij._1 *:* breeze.numerics.log((pij._1 + eps) /:/ (pij._2 + eps )) ) ).sum
	
	val nUnits = pia1.size
	val nPairs = pijab1.size

	( KL1 / nUnits, KL2 / nPairs, eps )
    }

    def printKL(outKL: PrintStream, bmStep: Int, nExtendedIterations: Int, aveKL: (Double, Double, Double)) = {
	if ( bmStep < 0 ) 
	  outKL.print("# step\tnExtendedIterations of MC\t<KLpi>\t<KLpij>\teps= %g\n".format(aveKL._3))
	else
	  outKL.print("%d\t%d\t%g %g\n".format(bmStep, nExtendedIterations, aveKL._1, aveKL._2))
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

    def printLog(logFiles: Tuple2[java.io.File, java.io.File], stateOrderString: String,
		interactions: MCMC.Interactions, 
		independentSamplings: MCMC.IndependentSamplings ) = {

	//val logFileInteractions = new File(outputFiles("outputDir"), logFiles._1)
	//val logFileIndpendentMCsamplings = new File(outputFiles("outputDir"), logFiles._2)
	val logFileInteractions = logFiles._1
	val logFileIndpendentMCsamplings = logFiles._2

	val outhJ =  new PrintStream(logFileInteractions)
    	MCMC.printInteractions(outhJ, stateOrderString, interactions )

	outhJ.close()

	val outSamples =  new PrintStream(logFileIndpendentMCsamplings)
	MCMC.printIndependentMCsamplings(outSamples, stateOrderString,
                		independentSamplings )
	outSamples.close()
    }

  }


  class BM ( 
	val outputFiles: HashMap[String, File],
	val stateOrderString: String,
	val fia: Array[DenseVector[Double]], val fijab: Array[DenseMatrix[Double]], 
	val regTerm: String,
	val propL1h: Double, val propL1J: Double,
	val lambdaPhi: Double, val lambdaPhij: Double,
      //val learningrate: Double,
	val betaV: Double, val betaM: Double, val eps: Double ) {

	//val betaV = 0.999, val betaM = 0.9, val learningRate = 0.000001 to 10, val eps = 1.0e-8 ) = {

    import miyazawa.potts.BM._

	require( regTerm == "L2" || regTerm == "GL1L2" || regTerm == "L1L2" )
	if ( regTerm == "L1L2" || regTerm == "GL1L2" ) {
		assert( propL1h >= 0.0 && propL1h <= 1.0)
		assert( propL1J >= 0.0 && propL1J <= 1.0)
	}
	val nUnits = fia.size
	val nStatesOfUnit = fia(0).size
	val nPairs = (nUnits * (nUnits - 1)) / 2

    def this( outputFiles: HashMap[String, File],
	stateOrderString: String,
	fia: Array[Array[Double]], fijab: Array[Array[Array[Double]]], 
	regTerm: String, 
	propL1h: Double, propL1J: Double,
	lambdaPhi: Double, lambdaPhij: Double,
	betaV: Double, betaM: Double, eps: Double ) =
		this( outputFiles, stateOrderString, fia.map(fi => DenseVector(fi)),
			//fijab.map( fij => DenseMatrix(fij:_*) ), 			// for scala 2.11
			  fijab.map( fij => new DenseMatrix(fij(0).size, fij.size, fij.flatten).t ), 	// for scala 2.13
				regTerm, propL1h, propL1J, lambdaPhi, lambdaPhij, betaV, betaM, eps )

    def runBM(
		optionInitialInteractions: Option[MCMC.Interactions],	// If optionInitialInteractions == None/null, fia will be used to generate hia.
		initialConfigurations: Array[Array[Byte]],	// If initialConfigurations(i) == null or _.size == 0, randomConfiguration will be called.

		nInitialIterationsPerUnit: Int = 100, 
		everyNIterationsPerUnit: Int = 10,

		// nSamples or initialNSamples and maxNSamples must be provided.
		// nSamples != 0 means initialNSamples=maxNSamples=nSamples.
		nSamples: Int = 0,
		initialNSamples: Int = 0,
		incrementSamples: Double = 0.0, 
		maxNSamples: Int = 0,

	      ////nIndependentMC: Int = 10,
		initialT: Double = 1.2,
		annealingRate: Double = 0.99,
		maxExtendedIterations: Int = 10,
		mcmcKernel: String = "MultiBlockMH", // or "MH" // or "Gibbs"

	      //fia: Array[DenseVector[Double]],
	      //fijab: Array[DenseMatrix[Double]],

	      //regTerm: String,
	      //propL1h: Double, propL1J: Double,
	      //lambdaPhi: Double, lambdaPhij: Double,
	      ////bmState: BMState, 

		gauge: String = "ungauged",
		optMethod: String = "ModAdam",
		learningRate: Double = 1.0e-2,
		minLearningRate: Double = 1.0e-6,
		maxLearningRate: Double = 10.0,
		rateDecrease: Double = 0.5,
		rateIncrease: Double = 1.2,
	
		epsForKL: Double = 1.0e-5,
		nBestKLs: Int = 1,
		minLearningStepsForBestKL: Int = 1000,
		minLearningSteps: Int = 1200,
		maxNoLearnings: Int = 100,
		logInterval: Int = 10
		) = {

	require( mcmcKernel == "MH" || mcmcKernel == "MultiBlockMH" || mcmcKernel == "Gibbs", 
			"Not supported: mcmcKernel == MH, MultiBlockMH or Gibbs" )

	// If the Ising gauge is used, convergence will be not good.
	//require( gauge == "ungauged" )

	if ( regTerm == "L1L2" ) {
		require( gauge == "ungauged" || gauge == "unused" )
	} else {
		require( gauge == "ungauged" || gauge == "unused" || 
			gauge == "phi_zeroSum" || gauge == "phi_ZeroSum" )	// gauge for phi and phij
	}

	val nIndependentMC = initialConfigurations.size

        val fileKL = new File(outputFiles("outputDir"), 
		(s"KL_of_each_step_${this.regTerm}_${this.propL1h}_${this.propL1J}_%.0e_%.0e_${gauge}_${optMethod}_${learningRate}_${betaV}_${betaM}.txt").format(this.lambdaPhi, this.lambdaPhij))
	val outKL =  new PrintStream(fileKL)
      //val epsForKL = 1.0e-5
	outKL.print("# %s  propL1h= %g propL1J= %g  lambdaPhi= %g  lambdaPhij= %g  %s  %s  learningRate= %g  betaV= %g  betaM= %g  eps= %g\n".format(
		this.regTerm, this.propL1h, this.propL1J, this.lambdaPhi, this.lambdaPhij, gauge, optMethod, 
		learningRate, this.betaV, this.betaM, this.eps) )
	if ( optMethod == "RPROP-LR" ) {
		outKL.print("# For RPROP-LR: minLearningRate= %g  maxLearningRate= %g  rateDecrease= %g  rateIncrease= %g\n".format(
			minLearningRate, maxLearningRate, rateDecrease, rateIncrease) )
	}

	outKL.print(s"# MCMC Kernel: ${mcmcKernel}\n")
	outKL.print(s"# nInitialIterationsPerUnit=${nInitialIterationsPerUnit}  everyNIterationsPerUnit=${everyNIterationsPerUnit}\n")

	outKL.print("# initialT = %g  annealingRate = %g\n".format(initialT, annealingRate))
	outKL.print("# maxExtendedIterations = %d\n".format(maxExtendedIterations))
	outKL.print("# nIndependentMC = %d".format(nIndependentMC))
	if ( nSamples != 0 )
			outKL.print("  nSamples for each MC run = %d\n".format(nSamples))
		else
			outKL.print("  initialNSamples = %d  incrementSamples = %g  maxNSamples = %d  for each MC run\n".format(
					initialNSamples, incrementSamples, maxNSamples))
	outKL.print("# minLearningStepsForBestKL = %d  minLearningSteps = %d  maxNoLearnings = %d\n".format(minLearningStepsForBestKL, minLearningSteps, maxNoLearnings) )
	printKL(outKL, bmStep = -1, nExtendedIterations = 0, (0.0, 0.0, epsForKL) ) 	// print a title

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

	val initialInteractions = { 
		val hJ =
		  if (optionInitialInteractions == null || optionInitialInteractions == None ) {
			initializehJ (fia )
		  } else {
			val interactions = optionInitialInteractions.get
			( interactions.hia, interactions.Jijab )
		  }
		val phiPhij = hJToPhi (hJ._1,  hJ._2, fia)
		val bmInteractions = BMInteractions(phiPhij._1, phiPhij._2)
		val (newInteractions, newBMInteractions) = newPhiThroughhJ(bmInteractions, fia, gauge)
		newInteractions
	}

	val (interactions, bmState) = BM.initialize(this.fia,
				initialInteractions,
				optMethod, 
				learningRate, 
				minLearningRate, maxLearningRate, 
				rateDecrease, rateIncrease, 
				this.betaV, this.betaM, this.eps, step = 0 )

	val runBM1Step = BM.runBM1Step(
		this.outputFiles,
		this.stateOrderString,
		_: MCMC.Interactions,
		_: Array[Array[Byte]],
	 	nInitialIterationsPerUnit,
		everyNIterationsPerUnit,
		_: Int,
	      //nIndependentMC,
		initialT,
		annealingRate,
		maxExtendedIterations,
		mcmcKernel,

		this.fia, this.fijab,
		this.regTerm,
		this.propL1h, this.propL1J,
		this.lambdaPhi, this.lambdaPhij,
		_: BM.BMState,
		optMethod
		)

	@annotation.tailrec
	def  iterateBMSteps(
		runBM1Step: (MCMC.Interactions, Array[Array[Byte]], Int, BM.BMState ) =>
			(MCMC, Array[MCMC.State], Int, MCMC.IndependentSamplings, MCMC.EnsembleAverages, BM.BMState) ,
		interactions: MCMC.Interactions, 
		initialConfigurations: Array[Array[Byte]],

		initialNSamples: Int,
		incrementSamples: Double, 
		maxNSamples: Int,

		bmState: BM.BMState,
		bestKLs: Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ],
		bestKL:  Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ] ) : 
		(MCMC.Interactions, Array[Array[Byte]]) = {
		//(MCMC.Interactions, BM.BMState) = {

	  val nSamples = { val n = initialNSamples + (incrementSamples * bmState.step).floor.toInt
			   if ( n > maxNSamples ) maxNSamples else n
			 }

	  val (mcmc, initialStates, nExtendedIterations, independentSamplings, ensembleAverages, nextBMState) = 
		runBM1Step( interactions, initialConfigurations, nSamples, bmState )

	  if(bmState.step == 0) {
		val pia = mcmc.proposedDistributions.proposedDistributionsAtAllUnits
		val nUnits = pia.size
		val pijab = new Array[DenseMatrix[Double]]( (nUnits * (nUnits - 1)) / 2)
		cfor(0)( i => i < pia.size, i => i + 1)( i => {
		  cfor(0)( j => j < i, j => j + 1)( j => {
			val (ij, a, b) = pairIndex(i, j, 0, 1)
			if ( a == 0 ) 
				pijab(ij) = pia(i) * pia(j).t
			else
				pijab(ij) = pia(j) * pia(i).t
		  } )
		} )
	  	val aveKL = calcKL(this.fia, this.fijab, pia, pijab, epsForKL)
	  	printKL(outKL, 0, 0, aveKL)
	  }

	  val aveKL = calcKL(this.fia, this.fijab, ensembleAverages.pia, ensembleAverages.pijab, epsForKL)
	  printKL(outKL, bmState.step, nExtendedIterations, aveKL)

	  val bestKLsRevised =
	    if ( bmState.step % logInterval == 0 ) {
		val logFileInteractions = new File(outputFiles("outputDir"), s"log_Interactions_hJ_${bmState.step}.txt")
		val logFileIndpendentMCsamplings = new File(outputFiles("outputDir"), s"log_MC_samples_${bmState.step}.fasta")
		BM.printLog( (logFileInteractions, logFileIndpendentMCsamplings),
			stateOrderString, interactions, independentSamplings)

		//if ( bestKLs.size < (logInterval / 2.0) ) deleteBestKLs(bestKLs) else bestKLs
		if ( bestKLs.size < (logInterval / 2.0) )
			fillNone(nBestKLs, Queue[Option[Tuple5[Double, Double, Int, Option[File], Option[File] ] ] ]() )
		else
			bestKLs
	    } else {
		bestKLs
	    }

	  val fileInteractions = new File(outputFiles("outputDir"), s"Interactions_hJ_${bmState.step}.txt")
	  val fileIndpendentMCsamplings = new File(outputFiles("outputDir"), s"MC_samples_${bmState.step}.fasta")

	  val (newBestKLs, oneOfBestKLs) = keepBestKLs(bestKLsRevised, Tuple5(aveKL._1, aveKL._2, bmState.step, Option(fileInteractions), Option(fileIndpendentMCsamplings) ) )
	  if( oneOfBestKLs ) {
 	    BM.printLog( (fileInteractions, fileIndpendentMCsamplings),
		stateOrderString, interactions, independentSamplings)
	  }

	  val (newBestKL, oneOfBestKL) = if ( nSamples >= maxNSamples && bmState.step >= minLearningStepsForBestKL ) {
						val none: Option[File] = None
						keepBestKLs(bestKL, Tuple5(aveKL._1, aveKL._2, bmState.step, none, none)) 
					 } else {
						(bestKL, false)
					 }

	  val lastConfigurations = independentSamplings.map( x => x(nSamples - 1).configuration ).toArray
	  if(bmState.step >= math.max(minLearningStepsForBestKL, minLearningSteps) && nSamples >= maxNSamples && bmState.step - newBestKL(0).get._3 >= maxNoLearnings ) {

	  	(interactions, lastConfigurations)

	  } else {
	    val (newInteractions, newBMState) = BM.newPhiThroughhJ(nextBMState, this.fia, gauge) 
	  
	    iterateBMSteps(runBM1Step, newInteractions, lastConfigurations, 
		initialNSamples, incrementSamples, maxNSamples,
		newBMState, newBestKLs, newBestKL)
	  }

	}
	
	val (newInteractions, newBMState) = if(gauge == interactions.gauge ) {
					(interactions, bmState)
				} else {
					BM.newPhiThroughhJ(bmState, this.fia, gauge) 
				}
	val results  = if ( nSamples == 0 ) {
			  iterateBMSteps(runBM1Step, newInteractions, initialConfigurations,
			  initialNSamples, incrementSamples, maxNSamples,
			  newBMState, bestKLs, bestKL)
		       } else {
			  iterateBMSteps(runBM1Step, newInteractions, initialConfigurations,
			  initialNSamples = nSamples, incrementSamples = 0.0, maxNSamples = nSamples,
			  newBMState, bestKLs, bestKL)
		       }

	outKL.close()

	results
    }


  }

}
