# getAtlasExperiment
#   - Download and return the SimpleList object representing a single
#   Expression Atlas experiment.
getAtlasExperiment <- function( experimentAccession ) {
    
    # Make sure the experiment accession is in the correct format.
    if( ! .isValidExperimentAccession( experimentAccession ) ) {

        stop( "Experiment accession not valid. Cannot continue." )
    }
    
    # URL to download Atlas data from.
    urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
    
    # Create filename for R data file.
    atlasExperimentSummaryFile <- paste( 
        experimentAccession,
        "-atlasExperimentSummary.Rdata", 
        sep = "" 
    )

    # Create full URL to download R data from.
    fullUrl <- paste( 
        urlBase, 
        experimentAccession, 
        atlasExperimentSummaryFile, 
        sep = "/" 
    )
    
    message( 
        paste( 
            "Downloading Expression Atlas experiment summary from:\n", 
            fullUrl
        ) 
    )

    # Create connection object for downloading data.
    connection <- url( fullUrl )
    
    # Try download, catching any errors
    loadResult <- try( load( connection ), silent = TRUE )
    
    # Quit if we got an error.
    if( class( loadResult ) == "try-error" ) {

        msg <- geterrmessage()

        warning( 
            paste( 
                paste( 
                    "Error encountered while trying to download experiment summary for",
                    experimentAccession,
                    ":"
                ),
                msg,
                paste( 
                    "There may not currently be an Expression Atlas experiment summary available for ",
                    experimentAccession,
                    ".\nPlease try again later, check the website at ",
                    paste( 
                        "http://www.ebi.ac.uk/gxa/experiments/",
                        experimentAccession,
                        sep = ""
                    ),
                    ", or email us at atlas-feedback@ebi.ac.uk",
                    sep = ""
                ),
                sep = "\n"
            )
        )

        return( )
    }
    
    # Close the connection.
    close( connection )

    # Make sure experiment summary object exists before trying to return it.
    getResult <- try( get( "experiment_summary" ) )

    if( class( getResult ) == "try-error" ) {
        
        stop( 
            "ERROR - Download appeared successful but no experiment summary object was found." 
        )
    }

    # If we're still here, things must have worked ok.
    message( 
        paste( 
            "Successfully downloaded experiment summary object for", 
            experimentAccession
        ) 
    )

    # Return the experiment summary.
    expSum <- get( "experiment_summary" )

    return( expSum )
}


# getAtlasData
#   - Download SimpleList objects for one or more Expression Atlas experiments
#   and return then in a list.
getAtlasData <- function( experimentAccessions ) {

    if( missing( experimentAccessions ) ) {

        stop( "Please provide a vector of experiment accessions to download." )
    }

    # Make sure experimentAccessions is a vector.
    if( ! is.vector( experimentAccessions ) ) {
        
        stop( "Please provide experiment accessions as a vector." )
    }

    # Only use valid accessions to download.
    experimentAccessions <- experimentAccessions[ 
        which( 
            sapply( 
                experimentAccessions, function( accession ) {
                    .isValidExperimentAccession( accession )
                }
            )
        )
    ]

    # The experimentAccessions vector is empty if none of the accessions are
    # valid. Just quit here if so.
    if( length( experimentAccessions ) == 0 ) {
        stop( "None of the accessions passed are valid ArrayExpress accessions. Cannot continue." )
    }
    
    # Go through each one and download it, creating a list.
    # So that the list has the experiment accessions as the names, use them as
    # names for the vector before starting to create the list.
    names( experimentAccessions ) <- experimentAccessions

    experimentSummaryList <- SimpleList( 
        
        lapply( experimentAccessions, function( experimentAccession ) {

            experimentSummary <- getAtlasExperiment( experimentAccession )
        } 
    ) )

    # Remove any null entries, i.e. experiments without R data files available.
    experimentSummaryList <- experimentSummaryList[ ! sapply( experimentSummaryList, is.null ) ]

    return( experimentSummaryList )
}


# .isValidExperimentAccession
#   - Return TRUE if experiment accession matches expected ArrayExpress
#   experiment accession pattern. Return FALSE otherwise.
.isValidExperimentAccession <- function( experimentAccession ) {

    if( missing( experimentAccession ) ) {
        
        warning( "Accession missing. Cannot validate." )

        return( FALSE )
    }

    if( !grepl( "^E-\\w{4}-\\d+$", experimentAccession ) ) {
        
        warning( 
            paste( 
                "\"", 
                experimentAccession, 
                "\" does not look like an ArrayExpress experiment accession. Please check.", 
                sep="" 
            ) 
        )
        
        return( FALSE )

    } else {

        return( TRUE )
    }
}


# searchAtlasExperiments
#   - Search (currently against ArrayExpress API) for datasets in Atlas matching given terms.
searchAtlasExperiments <- function( properties, species = NULL ) {
    
    # Quit if we don't have any search terms
    if( missing( properties ) ) {
        stop( "Please provide at least one search term." )
    }
    
    # If we've got something other than a character vector of properties, also quit.
    if( typeof( properties ) != "character" ) {
        stop( "Please provide search term(s) as a character vector." )
    }
    
    # Make properties URL friendly (e.g. replace spaces with %20
    properties <- sapply( properties, URLencode )

    # If we weren't passed a species, log this.
    if( missing( species ) ) {
    
        message( "No species was provided. Will search for data from all available species." )
    
    } else if( typeof( species ) != "character" ) {
       
        # Quit if the species isn't a character vector.
        stop( "Please provide species as a character vector." )
    
    } else if( length( species ) > 1 ) {
        
        # Only allow one species to be specified.
        stop( "More than one species found. You may only specify one species at a time." )
    }
    
    # ArrayExpress API base URL.
    aeAPIbase <- "http://www.ebi.ac.uk/arrayexpress/xml/v2/experiments?keywords="

    # Construct the query URL
    queryURL <- paste( 
        aeAPIbase,
        paste( properties, collapse = "+OR+" ),
        "&gxa=TRUE",
        sep = ""
    )
    
    # Add the species to the URL if we were passed one.
    if( !missing( species ) ) {
    
        species <- URLencode( species )

        queryURL <- paste( 
            queryURL, 
            "&species=", 
            species, 
            sep = "" 
        )
    }
    
    # Log search is beginning.
    message( "Searching for Expression Atlas experiments matching your query ..." )
    
    # Run the query and download the result.
    response <- GET( queryURL )
    
    # Make sure the HTTP request worked.
    if( status_code( response ) != 200 ) {
        stop( 
            paste( 
                "Error running query. Received HTTP error code", 
                status_code( response ),
                "from server. Please try again later. If you continue to experience problems please email atlas-feedback@ebi.ac.uk"
            )
        )
    } else {
        message( "Query successful." )
    }
    
    # Parse the XML document.
    parsedXML <- xmlParse( content( response ) )

    # Get the root node of the XML.
    allExpsNode <- xmlRoot( parsedXML )
    
    # Get the number of experiments we found.
    numExps <- xmlAttrs( allExpsNode )[ "total" ]

    # If there were no results, quit here.
    if( numExps == 0 ) {
        return( message( "No results found. Cannot continue." ) )
    
    } else {
        message( paste( "Found", numExps, "experiments matching your query." ) )
    }

    # Get a list of all the experiments from the root node.
    allExperiments <- xmlElementsByTagName( allExpsNode, "experiment" )

    # Pull out the title, accession, type and species of each experiment.
    resultsList <- lapply( allExperiments, function( experimentNode ) {

        expAcc <- xmlValue( xmlElementsByTagName( experimentNode, "accession" )$accession )

        expTitle <- xmlValue( xmlElementsByTagName( experimentNode, "name" )$name )

        species <- xmlValue( xmlElementsByTagName( experimentNode, "organism" )$organism )
        
        # Experiment type is e.g. microarray, RNA-seq, ...
        expType <- xmlValue( xmlElementsByTagName( experimentNode, "experimenttype")$experimenttype )
        
        # Return a list with this experiment's collected info.
        list( accession = expAcc, title = expTitle, species = species, expType = expType )

    } )
    
    # Create vectors of all accessions, experiment types, species, and titles.
    allAccessions <- sapply( resultsList, function( x ) { x$accession } )
    allExpTypes <- sapply( resultsList, function( x ) { x$expType } )
    allSpecies <- sapply( resultsList, function( x ) { x$species } )
    allTitles <- sapply( resultsList, function( x ) { x$title } )

    # Remove the names (all names are "experiment") so they don't show up later
    # and confuse things.
    names( allAccessions ) <- NULL
    names( allExpTypes ) <- NULL
    names( allSpecies ) <- NULL
    names( allTitles ) <- NULL
    
    # Create DataFrame containing the above results as columns.
    resultsSummary <- DataFrame( 
        Accession = allAccessions, 
        Species = allSpecies, 
        Type = allExpTypes, 
        Title = allTitles
    )

    # Sort the columns by species, type, then accession.
    resultsSummary <- resultsSummary[ order( resultsSummary$Species, resultsSummary$Type, resultsSummary$Accession ), ]
    
    # Return the DataFrame.
    return( resultsSummary )
}
