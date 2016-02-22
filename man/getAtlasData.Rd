\name{ getAtlasData }
\alias{ getAtlasData }
\title{ Download data from multiple Expression Atlas experiments }
\description{
    This function downloads Expression Atlas experiment summary SimpleList
    objects based on a vector of ArrayExpress experiment accessions, and
    returns a list containing these objects.
}
\usage{
    getAtlasData( experimentAccessions )
}
\arguments{
    \item{ experimentAccessions }{
    Vector of ArrayExpress accessions for experiments to be downloaded.
    }
}
\value{
    A list with one entry per experiment summary SimpleList object. Entries are
    named using the ArrayExpres accession of the respective experiment.
}
\examples{

    # Download some Expression Atlas data into a list.
    myExperimentSummaries <- getAtlasData( 
        c( 
            "E-GEOD-11175", 
            "E-MTAB-3007", 
            "E-GEOD-21070" 
        ) 
    )
}
