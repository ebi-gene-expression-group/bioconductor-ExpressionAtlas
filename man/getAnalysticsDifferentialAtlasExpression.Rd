\name{ getAnalysticsDifferentialAtlasExpression }
\alias{ getAnalysticsDifferentialAtlasExpression }
\title{ Download differential expression analytics data for an Expression Atlas experiment }
\description{
    This function downloads and returns differential expression analytics data for an 
    Expression Atlas experiment, based on the ArrayExpress accession of
    the experiment.
}
\usage{
    getAnalysticsDifferentialAtlasExpression( experimentAccession )
}
\arguments{
    \item{ experimentAccession }{
        ArrayExpress experiment accession e.g. "E-MTAB-10104"
    }
}
\value{
    A dataframe with \code{p.value} and \code{log2foldchange} values of the genes for each assay group contrast.
}
\examples{
    
    # Download the differential expression analytics data for experiment E-MTAB-10104
    mtab10104_dea <- getAnalysticsDifferentialAtlasExpression( "E-MTAB-10104" )

    # See the entries available
    dim(mtab10104_dea)
    # Prints out the following:
    # [1] 56748    10

}
