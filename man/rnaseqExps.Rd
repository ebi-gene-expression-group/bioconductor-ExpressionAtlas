\name{ rnaseqExps }
\alias{ rnaseqExps }
\title{ A SimpleList containing a dummy RNA-seq Expression Atlas dataset }
\description{
    This is a SimpleList object containing dummy data from an RNA-seq Expression
    Atlas experiment, to demonstrate a typical example of the reesults of using
    the \code{getAtlasData()} function for a subset of results obtained using
    \code{searchAtlasExperiments}. This object contains the first 10 rows of
    the original data, to save disk space. To get the full dataset, download it
    using \code{getAtlasData()}.
}
\usage{
    data( "rnaseqExps" )
}
\value{
    A SimpleList with one element per Expression Atlas dataset.
}
