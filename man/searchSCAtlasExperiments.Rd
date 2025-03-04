\name{ searchSCAtlasExperiments }
\alias{ searchSCAtlasExperiments }
\title{ Search for Sinle-Cell Expression Atlas experiments }
\description{
    This function accepts a query, and optionally a
    secondary filter, and then searches for matching Single-Cell 
    Expression Atlas experiments in the EBI RESTful Web Services API. 
    The function also performs ontology-based query search against the 
    Single Cell Expression Atlas cell-type wheel.
}
\usage{
    searchSCAtlasExperiments( query, secondaryFilter = NULL  )
}
\arguments{
    \item{ query }{
        Character of sample properties to search Atlas for. These can be
        biological characteristics, experimental treatments, species, etc.
    }
    \item{ secondaryFilter }{
        Optional, a second filter. For instance, if the query is "lung", the
        secondary filter can be "human" to search for human lung experiments.
    }
}
\value{
    A DataFrame containing the Expression Atlas accessions, the species, experiment
    types, and titles of Expression Atlas experiments matching the query.
}
\examples{
    # Search for pancreas experiments
    atlasRes_pancreas <- searchSCAtlasExperiments( query="pancreas")

    # Search for human lung experiments
    atlasRes_human_lung <- searchSCAtlasExperiments( query="lung", secondaryFilter = "human")
}
