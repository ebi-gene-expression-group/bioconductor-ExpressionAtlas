\name{ getNormalisedAtlasExpression }
\alias{ getNormalisedAtlasExpression }
\title{ Download normalised expression data for an Expression Atlas experiment }
\description{
    This function downloads and returns normalised expression data for a
    single Expression Atlas experiment, based on the ArrayExpress accession of
    the experiment.
}
\usage{
    getNormalisedAtlasExpression( experimentAccession, normalisation )
}
\arguments{
    \item{ experimentAccession }{
    ArrayExpress experiment accession e.g. "E-MTAB-4045"
    }
    \item{ normalisation }{
    normalisation method type e.g. "tpm", "fpkm"
    }
}
\value{
    A dataframe with expression values of the genes for each assay group. Entries are
    mean values of replicates.
}
\examples{
    
    # Download the Expression values for experiment E-MTAB-4045
    mtab4045_tpm <- getNormalisedAtlasExpression( "E-MTAB-4045", "tpm" )

    # See the entries available
    dim(results)
    # Prints out the following:
    # [1] 57147    45

}
