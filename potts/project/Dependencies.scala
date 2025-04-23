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

  lazy val biojava              = "org.biojava" % "biojava" % "5.3.0"
  lazy val biojavaStructure     = "org.biojava" % "biojava-structure" % "5.3.0"
  lazy val biojavaAlignment     = "org.biojava" % "biojava-alignment" % "5.3.0"
  lazy val biojavaCore          = "org.biojava" % "biojava-core" % "5.3.0"

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

  lazy val toolkit = "org.scala-lang" %% "toolkit" % "0.1.7"
  lazy val toolkitTest = "org.scala-lang" %% "toolkit-test" % "0.1.7" % Test

  lazy val spire		= "org.typelevel" %% "spire" % "0.18.0"
  lazy val spireMacros		= "org.typelevel" %% "spire-macros" % "0.18.0"
  lazy val breeze		= "org.scalanlp" %% "breeze" % "1.2"
  lazy val breezeNatives	= "org.scalanlp" %% "breeze-natives" % "1.2"
  lazy val breezeViz		= "org.scalanlp" %% "breeze-viz" % "1.2"

  lazy val netlibJava		= "com.googlecode.netlib-java" % "netlib-java" % "1.1"
  lazy val netlibNativeRef	= "com.github.fommil.netlib" % "netlib-native_ref-linux-x86_64" % "1.1"
  lazy val netlibNativeSystem	= "com.github.fommil.netlib" % "netlib-native_system-linux-x86_64" % "1.1"

  lazy val jfreechart	 	= "org.jfree" % "jfreechart" % "1.0.19"

  lazy val smileCore	 	= "com.github.haifengl"  % "smile-core" % "3.0.2"
  lazy val smileScala	 	= "com.github.haifengl" %% "smile-scala" % "3.0.2"
  lazy val javacpp		= "org.bytedeco" % "javacpp"   % "1.5.9"  classifier "linux-x86_64"
  lazy val openblas		= "org.bytedeco" % "openblas"  % "0.3.23-1.5.9" classifier "linux-x86_64"
  lazy val arpackNg		= "org.bytedeco" % "arpack-ng" % "3.9.0-1.5.9"  classifier "linux-x86_64"
  lazy val smileMKL		= "com.github.haifengl" % "smile-mkl" % "3.0.2"

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
/* for scala3 3.2 */
/* for scala3 3.3 */
/* for scala3 3.4 */
/* for scala3 3.6 */
  lazy val scalaTest = "org.scalatest" %% "scalatest" % "3.2.10"

  lazy val scalaCtic  = "org.scalactic" %% "scalactic" % "3.2.10"
  lazy val scalaCheck = "org.scalacheck" %% "scalacheck" % "1.14.1"
  lazy val scalaCollectionCompat = "org.scala-lang.modules" %% "scala-collection-compat" % "2.5.0"

  lazy val toolkit = "org.scala-lang" %% "toolkit" % "0.1.7"
  lazy val toolkitTest = "org.scala-lang" %% "toolkit-test" % "0.1.7" % Test

  lazy val spire		= "org.typelevel" %% "spire" % "0.18.0"
  lazy val spireMacros		= "org.typelevel" %% "spire-macros" % "0.18.0"
  lazy val breeze		= "org.scalanlp" %% "breeze" % "2.1.0"
  lazy val breezeNatives	= "org.scalanlp" %% "breeze-natives" % "2.1.0"
  lazy val breezeViz		= "org.scalanlp" %% "breeze-viz" % "2.1.0"

  lazy val netlibJava		= "com.googlecode.netlib-java" % "netlib-java" % "1.1"
  lazy val netlibNativeRef	= "com.github.fommil.netlib" % "netlib-native_ref-linux-x86_64" % "1.1"
  lazy val netlibNativeSystem	= "com.github.fommil.netlib" % "netlib-native_system-linux-x86_64" % "1.1"

  lazy val netlibBlas		= "dev.ludovic.netlib" % "blas" % "3.0.3"
  lazy val netlibLapack		= "dev.ludovic.netlib" % "lapack" % "3.0.3"
  lazy val netlibArpack		= "dev.ludovic.netlib" % "arpack" % "3.0.3"
  lazy val netlibArpackCombined	= "net.sourceforge.f2j" % "arpack_combined_all" % "0.1"

  lazy val jfreechart	 	= "org.jfree" % "jfreechart" % "1.0.19"

  lazy val smileCore	 	= "com.github.haifengl"  % "smile-core" % "3.0.2"
  lazy val smileScala	 	= "com.github.haifengl" %% "smile-scala" % "3.0.2"
  lazy val javacpp		= "org.bytedeco" % "javacpp"   % "1.5.11"  classifier "linux-x86_64"
  lazy val openblas		= "org.bytedeco" % "openblas"  % "0.3.28-1.5.11" classifier "linux-x86_64"
  lazy val arpackNg		= "org.bytedeco" % "arpack-ng" % "3.9.1-1.5.11"  classifier "linux-x86_64"
  lazy val smileMKL		= "com.github.haifengl" % "smile-mkl" % "3.0.2"

  lazy val biojavaStructure	= "org.biojava" % "biojava-structure" % "7.2.2"
  lazy val biojavaAlignment	= "org.biojava" % "biojava-alignment" % "7.2.2"
  lazy val biojavaCore		= "org.biojava" % "biojava-core" % "7.2.2"
  lazy val biojavaGenome	= "org.biojava" % "biojava-genome" % "7.2.2"
  lazy val biojavaModfinder	= "org.biojava" % "biojava-modfinder" % "7.2.2"
  lazy val biojavaOntology	= "org.biojava" % "biojava-ontology" % "7.2.2"

  lazy val forester		= "org.biojava.thirdparty" % "forester" % "1.005"

  lazy val openchart		= "openchart" % "openchart" % "1.4.2"


  lazy val evilplot		= "com.cibo" %% "evilplot" % "latest.release"	//"0.2.0"
  lazy val scalaView		= "com.github.darrenjw" %% "scala-view" % "latest.release"	//"0.6-SNAPSHOT"

//lazy val javacpp		= "org.bytedeco" % "javacpp" % "1.5.11"
//lazy val javacppPlatform	= "org.bytedeco" % "javacpp-platform" % "1.5.11"
  lazy val gsl			= "org.bytedeco" % "gsl" % "2.8-1.5.11"
  lazy val gslPlatform		= "org.bytedeco" % "gsl-platform" % "2.8-1.5.11"


/**/

}
