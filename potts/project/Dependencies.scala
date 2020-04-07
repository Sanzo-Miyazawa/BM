import sbt._

object Dependencies {
/* for scala 2.11
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.0.8"

  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.0"
  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.0.8"

  lazy val spire         = "org.typelevel" %% "spire" % "0.15.0"
  lazy val spireMacros   = "org.typelevel" %% "spire-macros" % "0.15.0"
  lazy val breeze        = "org.scalanlp" %% "breeze" % "0.13.2"
  lazy val breezeNatives = "org.scalanlp" %% "breeze-natives" % "0.13.2"
  lazy val breezeViz     = "org.scalanlp" %% "breeze-viz" % "0.13.2"

  lazy val jfreechart	 = "org.jfree" % "jfreechart" % "1.0.19"

  lazy val biojava		= "org.biojava" % "biojava" % "5.3.0"
  lazy val biojavaStructure	= "org.biojava" % "biojava-structure" % "5.3.0"
  lazy val biojavaAlignment	= "org.biojava" % "biojava-alignment" % "5.3.0"
  lazy val biojavaCore		= "org.biojava" % "biojava-core" % "5.3.0"

  lazy val javacpp		= "org.bytedeco" % "javacpp" % "1.3.2"
//lazy val javacppPresets	= "org.bytedeco" % "javacpp-presets" % "1.3"
//lazy val javacppPresetsPlatform	= "org.bytedeco" % "javacpp-presets-platform" % "1.3"
  lazy val gsl			= "org.bytedeco.javacpp-presets" % "gsl" % "2.2.1-1.3"
  lazy val gslPlatform		= "org.bytedeco.javacpp-presets" % "gsl-platform" % "2.2.1-1.3"
*/

/* for scala 2.13 */
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.0.8"

  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.0"
  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.0.8"

  lazy val spire		= "org.typelevel" %% "spire" % "0.17.0-M1"
  lazy val spireMacros		= "org.typelevel" %% "spire-macros" % "0.17.0-M1"
  lazy val breeze		= "org.scalanlp" %% "breeze" % "1.0"
  lazy val breezeNatives	= "org.scalanlp" %% "breeze-natives" % "1.0"
  lazy val breezeViz		= "org.scalanlp" %% "breeze-viz" % "1.0"

  lazy val jfreechart	 	= "org.jfree" % "jfreechart" % "1.0.19"

  lazy val biojava		= "org.biojava" % "biojava" % "5.3.0"
  lazy val biojavaStructure	= "org.biojava" % "biojava-structure" % "5.3.0"
  lazy val biojavaAlignment	= "org.biojava" % "biojava-alignment" % "5.3.0"
  lazy val biojavaCore		= "org.biojava" % "biojava-core" % "5.3.0"

  lazy val javacpp		= "org.bytedeco" % "javacpp" % "1.3.2"
//lazy val javacppPresets	= "org.bytedeco" % "javacpp-presets" % "1.3"
//lazy val javacppPresetsPlatform	= "org.bytedeco" % "javacpp-presets-platform" % "1.3"
  lazy val gsl			= "org.bytedeco.javacpp-presets" % "gsl" % "2.2.1-1.3"
  lazy val gslPlatform		= "org.bytedeco.javacpp-presets" % "gsl-platform" % "2.2.1-1.3"
/**/

}
