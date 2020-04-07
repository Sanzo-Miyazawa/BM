
package org.sanzo.potts

import java.io.PrintStream

import breeze.linalg.DenseMatrix
import breeze.linalg.max
import breeze.numerics.exp


object PairIntScore {

  // unpertubed state: Jij , pertubed state: Jij = 0

  // Free energy perturbation: < exp( -Jij ) >_{Z(h, J)} = Z(Jij=0) / Z(h, J)
  def logEnsAvExpNJij(jij: DenseMatrix[Double], ensemblePij: DenseMatrix[Double]) = {
	val nJij = - jij
	val maxNJij = breeze.linalg.max(nJij)
	val expNJij = breeze.numerics.exp(nJij - maxNJij )
	val avFrac = breeze.linalg.sum(expNJij *:* ensemblePij)
	val logAv = maxNJij + math.log(avFrac)
	logAv
  }

  def deltaLogLikelihoodJij(jij: DenseMatrix[Double], samplePij: DenseMatrix[Double], nSamples: Double, ensemblePij: DenseMatrix[Double] ) = {
	val logAv = logEnsAvExpNJij(jij, ensemblePij)
	val deltaL = (- breeze.linalg.sum(jij *:* samplePij) - logAv ) * nSamples

	deltaL
  }

  def chi2ValForLR(deltaMaxLogLikelihood : Double) = (- deltaMaxLogLikelihood * 2.0) 

  def chi2CriticalValue(p: Double, df: Double) = {
	val chi2 = new breeze.stats.distributions.ChiSquared(df)
	val chi2CV = chi2.inverseCdf(1.0 - p)

	chi2CV
  }

  def calcChi2ValForLRJij(jijab: Array[DenseMatrix[Double]], samplePijab: Array[DenseMatrix[Double]], nSamples: Double, ensemblePijab: Array[DenseMatrix[Double]] ) = {
	val nPairs = jijab.size
	val dLogLikelihoodJij = Array.range(0, nPairs).map( ij => {
				  deltaLogLikelihoodJij(jijab(ij), samplePijab(ij), nSamples, ensemblePijab(ij) )	
				} )
	val chi2ValForLRJij = Array.range(0, nPairs).map( ij => {
			      chi2ValForLR( dLogLikelihoodJij(ij) )
			} )
	val df = (jijab(0).rows - 1.0)*(jijab(0).cols - 1.0)
	(dLogLikelihoodJij, chi2ValForLRJij, df)
  }

  def sortPairIndexWithKey(dLogLikelihoodJij: Array[Double], chi2ValForLRJij: Array[Double] ) = {

	val nPairs = dLogLikelihoodJij.size
	val pairIndexWithKey = Array.range(0, nPairs).map ( ij => {
			val pairIJ = MCMC.inversePairIndex(ij)
			(ij, pairIJ._1, pairIJ._2, dLogLikelihoodJij(ij), chi2ValForLRJij(ij)) 
		} )
	val pairIndexWithKeySorted = pairIndexWithKey.sortWith( (x, y) => { if ( x._5 < y._5 ) true else false } )

	pairIndexWithKeySorted
  }

  def printPairIndexWithKeySorted(out: PrintStream, pairIndexWithKeySorted: Array[Tuple5[Int, Int, Int, Double, Double]] ) = {
	
	out.print("#serial\tindex\ti\tj\tdLogLikelihood\tchi2ValForLR\n")
	/*
	pairIndexWithKeySorted.foreach ( tuple => {
		val (ij, i, j, dLogLikelihoodJij, chi2ValForLRJij) = tuple
		out.print("%d\t%d %d\t%g\t%g\n".format( ij, i, j, dLogLikelihoodJij, chi2ValForLRJij) )
	} )
	*/
	val nPairs = pairIndexWithKeySorted.size
	Range(0, nPairs).foreach( serial => {
		val (ij, i, j, dLogLikelihoodJij, chi2ValForLRJij) = pairIndexWithKeySorted(serial)
		out.print("%d\t%d\t%d %d\t%g\t%g\n".format(serial, ij, i, j, dLogLikelihoodJij, chi2ValForLRJij) )
	} )
  }

}

