\name{ searchAtlasExperiments }
\alias{ searchAtlasExperiments }
\title{ Search for Expression Atlas experiments }
\description{
    This function accepts a query, and optionally a
    secondary filter, and then searches for matching Expression Atlas experiments
    in the EBI RESTful Web Services API
}
\usage{
    searchAtlasExperiments( query, secondaryFilter = NULL, detailed = FALSE  )
}
\arguments{
    \item{ query }{
        Character of sample properties to search Atlas for. These can be
        biological characteristics, experimental treatments, species, etc.
    }
    \item{ secondaryFilter }{
        Optional, a second filter. 
    }
    \item{ detailed }{
        If \code{TRUE}, it will perform a detailed search through Atlas experiment 
        metadata information, and it might take 10-20 seconds. Default: \code{FALSE}.
    }
}
\value{
    A DataFrame containing the Expression Atlas accessions, the species, experiment
    types, and titles of Expression Atlas experiments matching the query.
}
\examples{
    # Search for endoderm experiments (4 results)
    atlasRes_1 <- searchAtlasExperiments( query="endoderm")

    # Search for endoderm experiments in human (2 results)
    atlasRes_2 <- searchAtlasExperiments( query="endoderm", secondaryFilter = "human" )
    # it produces same results as:
    # atlasRes_2 <- searchAtlasExperiments( query="human", secondaryFilter = "endoderm" )

    # Download data for first experiment found.
    if ( nrow( atlasRes_2 ) == 1 ) {
        atlasData <- getAtlasData( atlasRes_2$Accession )
    } else {
        atlasData <- getAtlasData( atlasRes_2$Accession[1] )
    }

}
