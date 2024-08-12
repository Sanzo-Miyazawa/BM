
package org.sanzo.potts {

  object Util  {

        import scala.annotation.tailrec
        import scala.collection.immutable.ArraySeq
        import spire.syntax.cfor.cfor
        import breeze.linalg
        import breeze.linalg.{DenseVector, DenseMatrix}
        import scala.Array.ofDim

        import org.sanzo.potts.MCMC.pairIndex 

        def toProb(fa: Array[Double]): Array[Double] = {
                val total = fa.sum 
                val p = fa.map( a => a / total )
                p
        }
        def toProb(fia: Array[Array[Double]]): Array[Double] = {
                val f0 = new Array[Double](fia(0).size)

                val fa = fia.fold(f0)( (s, x) => s.zip(x).map( (si, xi) => si + xi ) ) 
                toProb(fa)
        }
        def faToProb(fa: ArraySeq[Double]): ArraySeq[Double] = {
                val total = fa.sum 
                val p = fa.map( a => a / total )
                p
        }
        def fiaToProb(fia: ArraySeq[ArraySeq[Double]]): ArraySeq[Double] = {
                val f0 = ArraySeq.range(0, fia(0).size).map(_ => 0.0)
                val fa = fia.fold(f0)( (s, x) => s.zip(x).map( (si, xi) => si + xi ) ) 
                faToProb(fa)
        }
        def toProb(fa: DenseVector[Double]): DenseVector[Double] = {
                val total = linalg.sum(fa)
                val p = fa.map( a => a / total )
                p
        }
        def toProb(fia: Array[DenseVector[Double]]): DenseVector[Double] = {
              //val fa = fia.fold(new DenseVector[Double](fia(0).size))(_ + _) 
                val fa = fia.fold(DenseVector.fill(fia(0).size)(0.0))(_ + _) 
                toProb(fa)
        }
        def toProb(fia: ArraySeq[DenseVector[Double]]): DenseVector[Double] = {
              //val fa = fia.fold(new DenseVector[Double](fia(0).size))(_ + _) 
                val fa = fia.fold(DenseVector.fill(fia(0).size)(0.0))(_ + _) 
                toProb(fa)
        }

        def energyOfRandomSeq(pa: DenseVector[Double],                               //nContactsPerSite: Double,
                              hia: ArraySeq[DenseVector[Double]], Jijab: ArraySeq[DenseMatrix[Double]] ):
                              Tuple3[Double, Double, Double] = {

            def isContact(i: Int, j: Int, thresholdContactRank: Double ) = { true }


          //set_nContactsPerSite(nContactsPerSite)  /**/

          //val f = totalAADelCounts
          //val p = countsToProb(f)

            val p = pa
            val nAADel = hia(0).size
            require( p.size == nAADel)
            val nSites = hia.size

            val hi_ = hia.map( hi => hi.t * p )        //meanhi(p, h)
            val dhia = hia.zip(hi_).map ( (ia, i_) => (ia - i_) )


            val Jija_ = Jijab.map( Jij => Jij * p )
            val Jij_a = Jijab.map( Jij => p.t * Jij )

            /* Jij = sum_{ai,aj}Jij(ai,aj) p(ai)p(aj) = sum_{ai,aj}Jji(aj,ai) p(ai)p(aj) = Jji */

            val Jij__ = Jija_.map( a_ => p.t * a_ )

            val dJijab = Jijab.zip(Jij__).map( (ijab, ij__) => ijab - ij__ )

            val deltaJija_ = Jija_.zip(Jij__).map ( (a_, ave) => a_ - ave )
            val deltaJij_a = Jij_a.zip(Jij__).map ( (_a, ave) => _a - ave )

            val dJija_ = ArraySeq.range(0, nSites).map( i =>
                           {
                             ArraySeq.range(0, nSites).map( j =>
                               {
                                 if ( i != j ) {
                                   val (ij, a, b) = pairIndex(i, j, 0, 1)
                                   if ( a == 0 ) {
                                     deltaJija_(ij) 
                                   } else {
                                     deltaJij_a(ij).t 
                                   }
                                 } else {
                                   null
                                 }
                               } )
                           } )

            val thresholdContactRank = scala.Double.PositiveInfinity  //val thresholdContactRank = this.thresholdContactRank( nContactsPerSite )

            @tailrec 
            def sum_i (i: Int, 
                                mean: Double, oneBodyVar: DenseVector[Double], twoBodyVar: DenseMatrix[Double] ): 
                                Tuple3[Double, DenseVector[Double], DenseMatrix[Double]] = {
                        
                        @tailrec
                        def sum_j (j: Int, 
                                meanJi: Double, oneBodyVar: DenseVector[Double], twoBodyVar: DenseMatrix[Double] ): 
                                Tuple3[Double, DenseVector[Double], DenseMatrix[Double]] = {

                                @tailrec
                                def sum_k (k: Int, 
                                        oneBodyVar: DenseVector[Double], twoBodyVar: DenseMatrix[Double] ): 
                                        Tuple2[DenseVector[Double], DenseMatrix[Double]] = {

                                        if( k < 0 ) {
                                                (oneBodyVar, twoBodyVar)
                                        } else if ( k == i || k == j ) {
                                                sum_k(k - 1, oneBodyVar, twoBodyVar)
                                        } else if ( isContact(i, j, thresholdContactRank ) && 
                                                        isContact(i, k, thresholdContactRank ) ) {
                                                val oneBodyV = ( dJija_(i)(j) *:* dJija_(i)(k) )

                                                sum_k(k - 1, oneBodyVar + oneBodyV, twoBodyVar )
                                        } else {
                                                sum_k(k - 1, oneBodyVar, twoBodyVar)
                                        }
                                }

                                if( j < 0 ) {
                                        (meanJi, oneBodyVar, twoBodyVar)        
                                } else if ( j == i ) {
                                        sum_j(j - 1, meanJi, oneBodyVar, twoBodyVar)
                                } else if ( isContact(i, j, thresholdContactRank ) ) {
                                        val (oneBodyVar_ij, twoBodyVar_ij) = sum_k(nSites - 1, oneBodyVar, twoBodyVar)

                                      //val ij = pairIndex(i, j)
                                        val (ij, a0, b1) = pairIndex(i, j, 0, 1)

                                        val oneBodyV = (dhia(i) *:* dJija_(i)(j)) * 2.0
                                        val twoBodyV = 
                                          if ( a0 == 0 ) {
                                            (dJijab(ij) *:* dJijab(ij)) * 0.5        
                                          } else {
                                            (dJijab(ij) *:* dJijab(ij)).t * 0.5  
                                          }

                                        sum_j(j - 1, meanJi + Jij__(ij), 
                                          oneBodyVar_ij + oneBodyV, twoBodyVar_ij + twoBodyV )
                                } else {
                                        sum_j(j - 1, meanJi, oneBodyVar, twoBodyVar)
                                }
                        }
                        if ( i < 0 ) {
                                (mean, oneBodyVar, twoBodyVar )
                        } else {
                                val (mJi, oneBodyVar_i, twoBodyVar_i) = sum_j(nSites -1, 0.0, oneBodyVar, twoBodyVar)

                                val oneBodyV = dhia(i) *:* dhia(i)

                                sum_i(i - 1, mean + hi_(i) + mJi * 0.5, oneBodyVar_i + oneBodyV, twoBodyVar_i )
                        }
            }

            val (mean, oneBodyVar, twoBodyVar) = sum_i(nSites - 1, 0.0,
                oneBodyVar = new DenseVector[Double](nAADel),
                twoBodyVar = new DenseMatrix[Double](nAADel, nAADel) )

            val meanPhi = - mean
            val varPhi = p.t * ( oneBodyVar + twoBodyVar * p )

            val totalPhi = meanPhi - varPhi
            (totalPhi, meanPhi, varPhi)

        /**************************************************
          //set_nContactsPerSite(nContactsPerSite)  /**/

            val f = totalAADelCounts
            val p = countsToProb(f)
            val nAADel = hia(0).size
            require( p.size == nAADel)

            val h = hia
            val J = Jijab
            val nSites = hia.size

            val hi = meanhi(p, h)
            val (j_ija_, j_ij) = meanJija_AndJij(p, J, nSites)

            val Jija_ = j_ija_
            val Jij = j_ij

            val oneBodyVar = new Array[Double](nAADel)
            val twoBodyVar = ofDim[Double](nAADel, nAADel)

            val thresholdContactRank = 0 //val thresholdContactRank = this.thresholdContactRank( nContactsPerSite )

            @tailrec 
            def sum_i (i: Int, 
                                //hi: Array[Double], Jij: Array[Double], Jija: Array[Array[Double]], 
                                mean: Double, oneBodyVar:  Array[Double], twoBodyVar: Array[Array[Double]] ): 
                                Tuple3[Double, Array[Double], Array[Array[Double]]] = {
                        
                        @tailrec
                        def sum_j (j: Int, 
                                        //hi: Array[Double], Jij: Array[Double], Jija: Array[Array[Double]], 
                                                meanJi: Double, oneBodyVar:  Array[Double], twoBodyVar: Array[Array[Double]] ): 
                                Tuple3[Double, Array[Double], Array[Array[Double]]] = {

                                @tailrec def sum_k (k: Int, 
                                        //hi: Array[Double], Jij: Array[Double], Jija: Array[Array[Double]], 
                                        oneBodyVar:  Array[Double], twoBodyVar: Array[Array[Double]] ): 
                                        Tuple2[Array[Double], Array[Array[Double]]] = {
                                        if( k < 0 ) {
                                                (oneBodyVar, twoBodyVar)        
                                        } else if ( k == i || k == j ) {
                                                sum_k(k - 1, oneBodyVar, twoBodyVar)
                                        } else if ( isContact(i, j, thresholdContactRank ) && 
                                                        isContact(i, k, thresholdContactRank ) ) {
                                                val ij = pairIndex(i, j)
                                                val ik = pairIndex(i, k)
                                                cfor(0) ( a => a < nAADel, a => a + 1) ( a => {
                                                        val dJija = Jija_(i)(j)(a) - Jij(ij)
                                                        val dJika = Jija_(i)(k)(a) - Jij(ik)
                                                        oneBodyVar(a) += dJija * dJika
                                                } )
                                                sum_k(k - 1, oneBodyVar, twoBodyVar )
                                        } else {
                                                sum_k(k - 1, oneBodyVar, twoBodyVar)
                                        }
                                }

                                if( j < 0 ) {
                                        (meanJi, oneBodyVar, twoBodyVar)        
                                } else if ( j == i ) {
                                        sum_j(j - 1, meanJi, oneBodyVar, twoBodyVar)
                                } else if ( isContact(i, j, thresholdContactRank ) ) {
                                        sum_k(nSites - 1, oneBodyVar, twoBodyVar)

                                        val ij = pairIndex(i, j)
                                        cfor(0) ( a => a < nAADel, a => a + 1) ( a => {
                                                val dhia = hia(i)(a) - hi(i)
                                                val dJija = Jija_(i)(j)(a) - Jij(ij)
                                                oneBodyVar(a) += dhia * dJija * 2
                                        } )
                                        cfor(0) ( a => a < nAADel, a => a + 1) ( a => {
                                        cfor(0) ( b => b < nAADel, b => b + 1) ( b => {
                                                val dJijab = Jijab(ij)(a)(b) - Jij(ij)
                                                twoBodyVar(a)(b) += dJijab * dJijab * 0.5
                                        } )
                                        } )
                                        sum_j(j - 1, meanJi + Jij(ij), oneBodyVar, twoBodyVar )
                                } else {
                                        sum_j(j - 1, meanJi, oneBodyVar, twoBodyVar)
                                }
                        }
                        if ( i < 0 ) {
                                (mean, oneBodyVar, twoBodyVar )
                        } else {
                                val mJi = sum_j(nSites -1, 0.0, oneBodyVar, twoBodyVar)._1

                                cfor(0) ( a => a < nAADel, a => a + 1) ( a => {
                                  val dhia = hia(i)(a) - hi(i)
                                  oneBodyVar(a) += dhia * dhia
                                } )

                                sum_i(i - 1, mean + hi(i) + mJi * 0.5, oneBodyVar, twoBodyVar )
                        }
            }
            val mean = sum_i(nSites - 1, 0.0, oneBodyVar, twoBodyVar)._1
            val meanPhi = - mean
            cfor(0) ( a => a < nAADel, a => a + 1) ( a => {
                        cfor(0) ( b => b < nAADel, b => b + 1) ( b => {
                                twoBodyVar(a)(b) = twoBodyVar(a)(b) * p(b)
                        } )
                        oneBodyVar(a) += twoBodyVar(a).sum
                        oneBodyVar(a) = oneBodyVar(a) * p(a)
            } )
            val varPhi = oneBodyVar.sum
            val totalPhi = meanPhi - varPhi
            (totalPhi, meanPhi, varPhi)
          ************************************************/

        }

        def energyOfRandomSeq(pa: DenseVector[Double],                               //nContactsPerSite: Double,
                              hia: Array[DenseVector[Double]], Jijab: Array[DenseMatrix[Double]] ):
                              Tuple3[Double, Double, Double] = {

            energyOfRandomSeq( pa,
                               hia = ArraySeq( hia* ),
                               Jijab = ArraySeq (Jijab* )
                             )
        }

      /*
        def energyOfRandomSeq(pa: DenseVector[Double],                               //nContactsPerSite: Double,
                              hia: ArraySeq[ArraySeq[Double]], Jijab: ArraySeq[ArraySeq[ArraySeq[Double]]] ):
                              Tuple3[Double, Double, Double] = {
          val hia_dvector = hia.map( a => DenseVector( a ) )
          val Jijab_dmatrix = Jijab.map( ab => new DenseMatrix(ab(0).size, ab.size, ab.flatten.toArray).t )

          energyOfRandomSeq(pa,                                                     //nContactsPerSite: Double,
                            hia_dvector, Jijab_dmatrix )
        }
      */

  }

}
