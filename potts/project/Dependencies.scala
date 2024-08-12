import sbt._

object Dependencies {
/* for scala 2.11
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.0.8"
  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.0.8"

  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.0"
  lazy val scalaCollectionCompat = "org.scala-lang.modules" %% "scala-collection-compat" % "2.5.0"

  lazy val junit = "com.novocode" % "junit-interface" % "0.11" % "test"

  lazy val spire         = "org.typelevel" %% "spire" % "0.15.0"
  lazy val spireMacros   = "org.typelevel" %% "spire-macros" % "0.15.0"
  lazy val breeze        = "org.scalanlp" %% "breeze" % "1.0"
  lazy val breezeNatives = "org.scalanlp" %% "breeze-natives" % "1.0"
  lazy val breezeViz     = "org.scalanlp" %% "breeze-viz" % "1.0"

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

/* for scala 2.13
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.0.8"
  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.0.8"

  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.0"
  lazy val scalaCollectionCompat = "org.scala-lang.modules" %% "scala-collection-compat" % "2.5.0"

  lazy val spire		= "org.typelevel" %% "spire" % "0.17.0"
  lazy val spireMacros		= "org.typelevel" %% "spire-macros" % "0.17.0"
  lazy val breeze		= "org.scalanlp" %% "breeze" % "1.2"
  lazy val breezeNatives	= "org.scalanlp" %% "breeze-natives" % "1.2"
  lazy val breezeViz		= "org.scalanlp" %% "breeze-viz" % "1.2"

  lazy val netlibJava		= "com.googlecode.netlib-java" % "netlib-java" % "1.1"
  lazy val netlibNativeRef	= "com.github.fommil.netlib" % "netlib-native_ref-linux-x86_64" % "1.1"
  lazy val netlibNativeSystem	= "com.github.fommil.netlib" % "netlib-native_system-linux-x86_64" % "1.1"

  lazy val jfreechart	 	= "org.jfree" % "jfreechart" % "1.0.19"

  lazy val biojava		= "org.biojava" % "biojava" % "5.3.0"
  lazy val biojavaStructure	= "org.biojava" % "biojava-structure" % "5.3.0"
  lazy val biojavaAlignment	= "org.biojava" % "biojava-alignment" % "5.3.0"
  lazy val biojavaCore		= "org.biojava" % "biojava-core" % "5.3.0"

  lazy val evilplot		= "com.cibo" %% "evilplot" % "latest.release"	//"0.2.0"
  lazy val scalaView		= "com.github.darrenjw" %% "scala-view" % "latest.release"	//"0.6-SNAPSHOT"

  lazy val javacpp		= "org.bytedeco" % "javacpp" % "1.3.2"
//lazy val javacppPresets	= "org.bytedeco" % "javacpp-presets" % "1.3"
//lazy val javacppPresetsPlatform	= "org.bytedeco" % "javacpp-presets-platform" % "1.3"
  lazy val gsl			= "org.bytedeco.javacpp-presets" % "gsl" % "2.2.1-1.3"
  lazy val gslPlatform		= "org.bytedeco.javacpp-presets" % "gsl-platform" % "2.2.1-1.3"
*/

/* for scala3 3.0.2 */
/* for scala3 3.1 */
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.2.10"

  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.2.10"
  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.1"
  lazy val scalaCollectionCompat = "org.scala-lang.modules" %% "scala-collection-compat" % "2.5.0"

  lazy val spire		= "org.typelevel" %% "spire" % "0.18.0"
  lazy val spireMacros		= "org.typelevel" %% "spire-macros" % "0.18.0"
  lazy val breeze		= "org.scalanlp" %% "breeze" % "2.1.0"
  lazy val breezeNatives	= "org.scalanlp" %% "breeze-natives" % "2.1.0"
  lazy val breezeViz		= "org.scalanlp" %% "breeze-viz" % "2.1.0"

  lazy val netlibJava		= "com.googlecode.netlib-java" % "netlib-java" % "1.1"
  lazy val netlibNativeRef	= "com.github.fommil.netlib" % "netlib-native_ref-linux-x86_64" % "1.1"
  lazy val netlibNativeSystem	= "com.github.fommil.netlib" % "netlib-native_system-linux-x86_64" % "1.1"

  lazy val jfreechart	 	= "org.jfree" % "jfreechart" % "1.0.19"

  lazy val biojava		= "org.biojava" % "biojava" % "6.0.1"	//"5.4.0"
  lazy val biojavaStructure	= "org.biojava" % "biojava-structure" % "6.0.1"	//"5.4.0"
  lazy val biojavaAlignment	= "org.biojava" % "biojava-alignment" % "6.0.1"	//"5.4.0"
  lazy val biojavaCore		= "org.biojava" % "biojava-core" % "6.0.1"	//"5.4.0"

  lazy val evilplot		= "com.cibo" %% "evilplot" % "latest.release"	//"0.2.0"
  lazy val scalaView		= "com.github.darrenjw" %% "scala-view" % "latest.release"	//"0.6-SNAPSHOT"

  lazy val javacpp		= "org.bytedeco" % "javacpp" % "1.3.2"
//lazy val javacppPresets	= "org.bytedeco" % "javacpp-presets" % "1.3"
//lazy val javacppPresetsPlatform	= "org.bytedeco" % "javacpp-presets-platform" % "1.3"
  lazy val gsl			= "org.bytedeco.javacpp-presets" % "gsl" % "2.2.1-1.3"
  lazy val gslPlatform		= "org.bytedeco.javacpp-presets" % "gsl-platform" % "2.2.1-1.3"

/**/

}
