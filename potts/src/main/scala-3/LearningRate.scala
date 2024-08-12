

package org.sanzo.potts {

  import java.io.PrintStream

  class LearningRateForRPROPLR (
        val minLearningRate: Double = 1.0e-06, val maxLearningRate: Double = 10.0,
        val rateDecrease: Double = 0.5,        val rateIncrease: Double = 1.2) {

    assert( rateIncrease > 1.0 && rateDecrease < 1.0 )

    def printParameters(out: PrintStream): Unit = {
         out.print("# For RPROP-LR: learningRateForRPROPLR.minLearningRate= %g  maxLearningRate= %g  rateDecrease= %g  rateIncrease= %g\n".format(
                        this.minLearningRate, this.maxLearningRate,
                        this.rateDecrease, this.rateIncrease) )
    }

    def lr(r: Double, minLR: Double = this.minLearningRate , maxLR: Double = this.maxLearningRate) = {
    
        if ( r > maxLR )
            maxLR
        else if ( r < minLR )
            minLR
        else
            r
    }
  }

  class LearningRate(val warmupEpochs: Double = 0.0, val maxLR: Double = 0.0, val maxLREpochs: Double = scala.Int.MaxValue.toDouble,
                     val coolingMethod: String = "",
              // for coolingMethod = "rate" 
                     val coolRate: Double = 0.98, val minLR: Double = 0.0,
              // for coolingMethod = "(1+at)^b"
                     val constPerEpoch: Double = 0.1, val power: Double = -0.5,
              // for coolingMethod = "1/(1+at)"
                     val epochsUntilHalf: Double = 10.0, 
              //
                     val stepsPerEpoch: Double = 1.0, val stepOffset: Double = 0.0
                     ) {

    require ( stepsPerEpoch > 0.0 )
    require ( warmupEpochs >= 0.0 )
    require ( maxLR >= 0.0 )        // ( maxLR > 0.0 )
    require ( maxLREpochs >= 0.0 )

    // learning rate decay may be not good for SGD methods such as Adam, because KL increases although interactions don't change much. 
  //method match {
    coolingMethod match {
        case "rate" => require( coolRate <= 1.0 && ( (coolRate > 0.0 && minLR >= 0.0) || (coolRate == 0.0 && minLR > 0.0) ) )
        case "(1+at)^b" => require( constPerEpoch > 0.0 && power < 0.0 )
        case "1/(1+at)" => require( epochsUntilHalf > 0.0 )                  // Obsolete
        case ""         => require(maxLREpochs >= scala.Int.MaxValue.toDouble, "Must: maxLREpochs >> 0" )
        case _          => sys.error ( s"*** No such a learning rate scheduler: ${coolingMethod}\n" ) 
    }

    def printParameters(out: PrintStream): Unit = {

                out.print("# learningRate.coolingMethod= %s  stepsPerEpoch= %f  warmupEpochs= %f  maxLR= %g  maxLREpochs= %g".format(
                  this.coolingMethod, this.stepsPerEpoch, this.warmupEpochs,
                  this.maxLR, this.maxLREpochs ) )
                if ( this.coolingMethod == "rate" ) {
                  out.print("  coolRate= %f  minLR= %g".format(this.coolRate, this.minLR) )
                } else if ( this.coolingMethod == "(1+at)^b" ) {
                  out.print("  constPerEpoch= %f  power= %f".format(this.constPerEpoch, this.power) )
                } else if ( this.coolingMethod == "1/(1+at)" ) {
                  out.print("  epochsUntilHalf= %f".format(this.epochsUntilHalf) )
                }
                out.print("\n")
    }

    def newLearningRate(warmupEpochs: Double = this.warmupEpochs, 
                        maxLR: Double = this.maxLR, maxLREpochs: Double = this.maxLREpochs,
                      //method: String = this.method
                        coolingMethod: String = this.coolingMethod,
                        coolRate: Double = this.coolRate, minLR: Double = this.minLR,
                        constPerEpoch: Double = this.constPerEpoch, power: Double = this.power,
                        epochsUntilHalf: Double = this.epochsUntilHalf, 
                        stepsPerEpoch: Double = this.stepsPerEpoch,
                        stepOffset: Double = this.stepOffset
                        ) = {
        new LearningRate(warmupEpochs, maxLR, maxLREpochs,
                       //method,
                         coolingMethod,
                         coolRate, minLR,
                         constPerEpoch, power,
                         epochsUntilHalf,
                         stepsPerEpoch, stepOffset)
    }

    def lr(step: Int, stepsPerEpoch: Double = this.stepsPerEpoch,
                   warmupEpochs: Double = this.warmupEpochs,
                   maxLR: Double = this.maxLR, maxLREpochs: Double = this.maxLREpochs,
                 //method: String = this.method
                   coolingMethod: String = this.coolingMethod,
                   coolRate: Double = this.coolRate, minLR: Double = this.minLR,
                   constPerEpoch: Double = this.constPerEpoch, power: Double = this.power,
                   epochsUntilHalf: Double = this.epochsUntilHalf,
                   stepOffset: Double = this.stepOffset ) = {

        val dStep = step - this.stepOffset
        require ( dStep > 0 )
        require ( stepsPerEpoch > 0.0 )
        require ( warmupEpochs >= 0.0 )
        require ( maxLR > 0.0 )
        require ( maxLREpochs >= 0.0 )
        val epoch = dStep / stepsPerEpoch
        if ( epoch <= warmupEpochs ) {
            maxLR * epoch / warmupEpochs
        } else if ( epoch - warmupEpochs <= maxLREpochs ) {
            maxLR
        } else {
            coolingMethod match {
              case "(1+at)^b" => {
                require(constPerEpoch > 0.0 && power < 0.0 )
                val rate = scala.math.pow(1.0 + constPerEpoch * (epoch - warmupEpochs - maxLREpochs) , power)
                maxLR * rate
              }
              case "1/(1+at)" => { 
                require( epochsUntilHalf > 0.0 )
                val rate = 1.0 / ( 1.0 + (epoch - warmupEpochs - maxLREpochs ) / epochsUntilHalf )
                maxLR * rate
              }
              case "rate" => {
                require( coolRate <= 1.0 && ( (coolRate > 0.0 && minLR >= 0.0) || (coolRate == 0.0 && minLR > 0.0) ) )
                require( maxLR >= minLR )
                if ( coolRate >= 1.0 ) {
                    maxLR
                } else if ( coolRate <= 0.0 ) {
                    minLR
                } else {
                    val rate = scala.math.pow(coolRate, epoch - warmupEpochs - maxLREpochs)
                    math.max(maxLR * rate, minLR)
                  //maxLR * rate + minLR * (1.0 - rate)
                }
              }
              case _ => sys.error ( "*** No such a learning rate scheduler: ${coolingMethod}\n" ) 
            }
        }
    }

  }

}
