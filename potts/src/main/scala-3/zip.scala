
package org.sanzo.potts {

  package zip {

    import java.io.File
    import java.io.{FileInputStream, FileOutputStream}
    import java.io.{InputStream, OutputStream}
    import java.util.zip.{GZIPInputStream, GZIPOutputStream}

    import scala.sys.process.stderr

    def gunzipInputStream( file: File, force: Boolean = false ): InputStream = {
    //val absoluteFile = file.getAbsoluteFile()
      val name = file.getName()
      val inputStream = new FileInputStream(file)

      val zipName = ("""(^.*[\.-]gz$|^.*[\.-_]z$|^.*\.t[ag]z$|.*\.tar\.Z$)""").r

    //if ( name.endsWith(".gz") || force ) {
      if ( zipName.findFirstMatchIn(name).nonEmpty || force ) { 
        GZIPInputStream ( inputStream )
      } else {
        inputStream 
      }
    }
    def gunzipInputStream( inputStream: InputStream ): InputStream = {
      val force = true
      if ( force ) {
        GZIPInputStream ( inputStream )
      } else {
        inputStream
      }
    }

    def gzipInputStream( file: File, force: Boolean = false ): InputStream = {
      gunzipInputStream( file, force )
    }
    def gzipInputStream( inputStream: InputStream ): InputStream = {
      gunzipInputStream( inputStream )
    }

    def gzipOutputStream( outputStream: OutputStream ): GZIPOutputStream = {
      GZIPOutputStream( outputStream )
    }
    def gzipOutputStream( file: File ): GZIPOutputStream = {
      GZIPOutputStream( FileOutputStream(file) )
    }

  }
}
