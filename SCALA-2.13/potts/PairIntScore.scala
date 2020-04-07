
package org.sanzo.potts

import java.io.PrintStream

import breeze.linalg.DenseMatrix
import breeze.linalg.max
import breeze.numerics.exp


object PairIntScore {

  // unpertubed state: Jij , pertubed state: Jij = 0

  def logAvExpNJij(jij: DenseMatrix[Double], pij: DenseMatrix[Double]) = {
	val nJij = - jij
	val maxNJij = breeze.linalg.max(nJij)
	val expNJij = breeze.numerics.exp(nJij - maxNJij )
	val avFrac = breeze.linalg.sum(expNJij *:* pij)
	val logAv = maxNJij + math.log(avFrac)
	logAv
  }

  def deltaLogLikelihood(jij: DenseMatrix[Double], pij: DenseMatrix[Double], nSamples: Double) = {
	val logAv = logAvExpNJij(jij, pij)
	val deltaL = - nSamples * breeze.linalg.sum(jij *:* pij) - logAv

	deltaL
  }

  def chi2ValForLR(deltaMaxLogLikelihood : Double) = (- deltaMaxLogLikelihood * 2.0) 

  def chi2CriticalValue(p: Double, df: Double) = {
	val chi2 = new breeze.stats.distributions.ChiSquared(df)
	val chi2CV = chi2.inverseCdf(1.0 - p)

	chi2CV
  }

  def calcChi2ValForLR(jijab: Array[DenseMatrix[Double]], pijab: Array[DenseMatrix[Double]], nSamples: Double) = {
	val nPairs = jijab.size
	val dLogLikelihoodij = Array.range(0, nPairs).map( ij => {
				  deltaLogLikelihood(jijab(ij), pijab(ij), nSamples)	
				} )
	val chi2ValForLRij = Array.range(0, nPairs).map( ij => {
			      chi2ValForLR( dLogLikelihoodij(ij) )
			} )
	val df = (jijab(0).rows - 1.0)*(jijab(0).rows - 1.0)
	(dLogLikelihoodij, chi2ValForLRij, df)
  }

  def sortPairIndexWithKey(dLogLikelihoodij: Array[Double], chi2ValForLRij: Array[Double] ) = {

	val nPairs = dLogLikelihoodij.size
	val pairIndexWithKey = Array.range(0, nPairs).map ( ij => {
			val pairIJ = MCMC.inversePairIndex(ij)
			(ij, pairIJ._1, pairIJ._2, dLogLikelihoodij(ij), chi2ValForLRij(ij)) 
		} )
	val pairIndexWithKeySorted = pairIndexWithKey.sortWith( (x, y) => { if ( x._5 < y._5 ) true else false } )

	pairIndexWithKeySorted
  }

  def printPairIndexWithKeySorted(out: PrintStream, pairIndexWithKeySorted: Array[Tuple5[Int, Int, Int, Double, Double]] ) = {
	
	val nPairs = pairIndexWithKeySorted.size
	out.print("#index\ti\tj\tddLogLikelihood\tchi2ValForLR\n")
	pairIndexWithKeySorted.foreach ( tuple => {
		val (ij, i, j, dLogLikelihoodij, chi2ValForLRij) = tuple
		out.print("%d\t%d %d\t%g\t%g\n".format( ij, i, j, dLogLikelihoodij, chi2ValForLRij) )
	} )
  }

}

