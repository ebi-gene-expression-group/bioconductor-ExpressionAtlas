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
                    ", or contact us at https://www.ebi.ac.uk/about/contact/support/gxa",
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
        stop( "None of the accessions passed are valid ArrayExpress/BioStudies accessions. Cannot continue." )
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
#   - Return TRUE if experiment accession matches expected ArrayExpress/BioStudies
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
                "\" does not look like an ArrayExpress/BioStudies/Atlas experiment accession. Please check.",
                sep=""
            )
        )

        return( FALSE )

    } else {

        return( TRUE )
    }
}


# searchAtlasExperiments
#   - Search (currently against BioStudies API) for datasets in Atlas matching given terms.
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

    # BioStudies API base URL.
    bioAPIbase <- "http://www.ebi.ac.uk/biostudies/api/v1/search?query="

    # Construct the query URL
    queryURL <- paste(
        bioAPIbase,
        paste( properties, collapse = "" ),
        "&gxa=TRUE",  #"&link_type=gxa"
        sep = ""
    )

    # Add the species to the URL if we were passed one.
    if( !missing( species ) ) {

        species <- URLencode( species )

        queryURL <- paste(
            queryURL,
            "&organism=",
            species,
            sep = ""
        )
    }

    # Add page size limit to BioStudies API query.
    page_size = 100
    queryURL <- paste(
        queryURL,
        "&pageSize=",
        page_size,
        sep = ""
    )


    # Log search is beginning.
    message( "Searching for Expression Atlas experiments matching your query ..." )

    # Check the query works with httr
    response <- GET( queryURL )

    # Make sure the HTTP request worked.
    if( status_code( response ) != 200 ) {
        stop(
            paste(
                "Error running query. Received HTTP error code",
                status_code( response ),
                "from server. Please try again later. If you continue to experience problems please contact us at https://www.ebi.ac.uk/about/contact/support/gxa"
            )
        )
    } else {
        message( "Query successful." )
    }

    ## Parse the JSON document.
    parsedJSON <- fromJSON( txt = queryURL )

    ## Get the number of experiments we found.
    numExps <- as.numeric( parsedJSON$totalHits )

    # If there were no results, quit here.
    if( numExps == 0 ) {
        return( message( "No results found. Cannot continue." ) )

    } else {
        message( paste( "Found", numExps, "experiments matching your query." ) )
    }

    # check that page size is correct, otherwise quit here.
    if( as.numeric( parsedJSON$pageSize ) != page_size ) {
        return( message( "No results found. Cannot continue." ) )
    }
    # check if total hits is exact number
    if( parsedJSON$isTotalHitsExact != TRUE ) {
        message( "WARNING: Total number of hits reported by BioStudies is not exact" )
    } 

    # traverse BioStudies pages to get the list of experiments
    allExperiments <- c()
    if ( numExps < page_size || numExps%%page_size != 0 ){
        modulus <- 1
    } else {
        modulus <- 0
    }

    number_pages <- floor( numExps/page_size ) + modulus
    for ( page_number in  1:number_pages ){

        pageQueryURL <- paste(
            queryURL,
            "&page=",
            page_number,
            sep = ""
        )
        page_acc <- fromJSON( txt = pageQueryURL )$hits$accession
        allExperiments <- append(allExperiments, page_acc)
    }

    # check the number of experiments is correct, otherwise quit here.
    if( numExps != length(allExperiments) ) {
        return( message( "Number of experiments obtained is not correct. Cannot continue." ) )
    }


    # Create vectors of all accessions, experiment types, species, titles and conect errors.
    allSpecies <- rep( NA, numExps )
    allExpTypes <- rep( NA, numExps )
    allTitles <- rep( NA, numExps )
    allConnErr <- rep( FALSE, numExps )


    # Pull out the title, accession, type and species of each experiment.
    message( paste( "Retrieving information from", length( allExperiments ), "experiments..." ) )
    for ( acc in 1:length( allExperiments ) ) {

        accQueryURL <- paste(
            "http://www.ebi.ac.uk/biostudies/api/v1/studies/",
            allExperiments[ acc ],
            sep = ""
        )
        # Make sure the HTTP request works
        if( status_code( GET( accQueryURL ) ) != 200 ) {
        
            allConnErr[ acc ] <- TRUE
        }
        else{ 
            accQuery <- fromJSON( txt = accQueryURL )

            # get indexes of title, the experiment type and species

            pos_title <-  which( accQuery$section$attributes$name %in% "Title" )
            if ( length( pos_title ) == 1 ){
                allTitles[ acc ] <- accQuery$section$attributes$value[ pos_title ]
            } else if ( length( pos_title ) > 1  ) {
                allTitles[ acc ] <- accQuery$section$attributes$value[ pos_title[1] ]
            } 

            pos_exp_type <-  which( accQuery$section$attributes$name %in% "Study type" )
            if ( length( pos_exp_type ) == 1 ){
                allExpTypes[ acc ] <- accQuery$section$attributes$value[ pos_exp_type ]
            } else if ( length( pos_exp_type) > 1  ) {
                # attept to match to a valid Atlas Experiment Type
                allExpTypes[ acc ] <- getEligibleAtlasExperiment( accQuery$section$attributes$value[ pos_exp_type ] )
            } 

            pos_species <-  which( accQuery$section$attributes$name %in% "Organism" )
            if ( length( pos_species ) == 1 ){
                allSpecies[ acc ] <- accQuery$section$attributes$value[ pos_species ]
            } else if ( length( pos_species ) > 1  ) {
                # return the first species
                allSpecies[ acc ] <- accQuery$section$attributes$value[ pos_species[1] ]  
            } 
        }

    }

    message( "Retrieving information completed."  )

    # Create DataFrame containing the above results as columns.
    resultsSummary <- DataFrame(
        Accession = allExperiments,
        Species = allSpecies,
        Type = allExpTypes,
        Title = allTitles,
        ConnectionError = allConnErr
    )

    # Remove experiments with connection errors.
    if ( any( resultsSummary$ConnectionError ) ) {
        message( "WARNING: One or more studies removed due to Error running query - received HTTP error code" )
        resultsSummary <- resultsSummary[ !resultsSummary$ConnectionError, ]
    }
  
    # Sort the columns by species, type, then accession.
    resultsSummary <- resultsSummary[ order( resultsSummary$Species, resultsSummary$Type, resultsSummary$Accession ), ]

    # Return the DataFrame.
    return( resultsSummary[ c( 'Accession', 'Species', 'Type', 'Title' ) ] )
}


eligibleAtlasExperiments <- c (
    "transcription profiling by array",
    "microRNA profiling by array",
    "antigen profiling",
    "proteomic profiling by mass spectrometer",
    "RNA-seq of coding RNA",
    "RNA-seq of non coding RNA",
    "RNA-seq of total RNA",
    "RNA-seq of coding RNA from single cells",
    "RNA-seq of non coding RNA from single cells"
)


getEligibleAtlasExperiment <- function( experiment_list, valid_experiments = eligibleAtlasExperiments ) {
    for (exptype in experiment_list){
        if ( exptype %in% valid_experiments ){
            # report the first in the list that is valid
            return( exptype )
            break
        } else {
            # otherwise will return the first in the list
        }
    }
    return( experiment_list[1] )
}


# Internal helper to get Experiment Type from accession
.getExperimentType <- function(experimentAccession) {

    # Base URL for downloading data.
    urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
    
    # Define filenames for configuration files.
    configFile <- paste(experimentAccession, "-configuration.xml", sep = "")

    # Define ftpUrl for configuration file.
    ftpUrl <- paste(urlBase, experimentAccession, configFile, sep = "/")

    # Read the XML file directly from the FTP URL
    message("Downloading XML file from FTP...")
    xmlContent <- tryCatch({
        read_xml(ftpUrl)
    }, error = function(e) {
        stop("Failed to download or read the XML file: ", e$message)
    })
    
    # Extract the <configuration> node and the experimentType attribute
    configurationNode <- xml_find_first(xmlContent, "/configuration")
    if (is.na(configurationNode)) {
        stop("No <configuration> node found in the XML file.")
    }
    
    # Get the value of the experimentType attribute
    experimentType <- xml_attr(configurationNode, "experimentType")
    if (is.na(experimentType)) {
        stop("The 'experimentType' attribute is missing in the <configuration> node.")
    }
    
    return(experimentType)
}

getNormalisedAtlasExpression <- function(experimentAccession, normalisation = "tpm") {

    # Ensure the experiment accession is in the correct format.
    if (!.isValidExperimentAccession(experimentAccession)) {
        stop("Experiment accession not valid. Cannot continue.")
    }

    expType <- .getExperimentType(experimentAccession)

    if (expType %in% c("rnaseq_mrna_baseline" )){
        message(paste(experimentAccession," is ", expType , ", will continue downloading data"))
    } else {
        stop("Experiment type not contain normalised expression data.")
    }

    # Ensure the normalisation is valid.
    fileType <- match.arg(normalisation, c("tpm", "fpkm"))

    # Base URL for downloading data.
    urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

    # Define filenames for TPM and FPKM files.
    expressionFile <- paste(experimentAccession, "-", fileType, "s.tsv", sep = "")

    # Create full URLs for TPM and FPKM files.
    expressionUrl <- paste(urlBase, experimentAccession, expressionFile, sep = "/")

    # Initialise a list to hold the results.
    results <- list()

    # Download and read TPM file if requested or in "both" mode.
    message(paste("Downloading expression file from:\n", expressionUrl))
    results <- tryCatch({
        read.table(expressionUrl, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "", comment.char = "#")
    }, error = function(e) {
        warning(paste("Failed to download or read expression file for", experimentAccession, ":", e$message))
        NULL
    })

    # Remove NULL results for files that could not be downloaded.
    results <- Filter(Negate(is.null), results)

    # Check if any files were successfully downloaded.
    if (length(results) == 0) {
        stop("ERROR - Could not download the requested data. Please check the experiment type is RNAseq baseline or try again later.")
    }

    # Message indicating success.
    if ("expression" %in% names(results)) {
        message(paste("Successfully downloaded expression data for", experimentAccession))
    }

    # TPMs and FPKMs can have multiple values in expression comma-separated columns (min, lower quartile, mean, higher quartile and max)
    get_mean_value <- function(col) {
        # Split by commas and return the third value
        sapply(strsplit(col, ","), function(x) x[3])
    }


    # the first two columns are GeneID and Gene.Name, fetching 3rd column as it contains mean values. 
    results[ ,3:ncol(results)] <- lapply(results[ ,3:ncol(results)], get_mean_value)

    return(results)
}


# ----- add visualisation support ----

# at
# getNormalisedAtlasExpression ( experimentAccession , normalisation='TPM' ) { }
# heatmapAtlasExperiment ( experimentAccession  )
#    input is the output of 'getAtlasExperiment' OR 'getNormalisedAtlasExpression'
#    in the first case it will normalise the raw counts to CPM
#    the heatmpa should show Expression Units
# possible packages -> complexHeatmap , or ggpplot2

# getAnalysticsDifferentialAtlasExpression ( experimentAccession ) {}
# volcanoDifferentialAtlasExperiment ( experimentAccession  )
# input is the output of 'getAnalysticsDifferentialAtlasExpression'


# ----- add single-cell support ----

# pm
#searchSCAtlasExperiments <- function( properties, species = NULL ) { }

getAtlasSCExperiment <- function( experimentAccession ) {

    # Make sure the experiment accession is in the correct format.
    if( ! .isValidExperimentAccession( experimentAccession ) ) {

        stop( "Experiment accession not valid. Cannot continue." )
    }

    # URL to download Atlas data from.
    urlBase <- "http://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/sc_experiments"

    # Create filename for anndata file.
    atlasAnnDataFile <- paste(
        experimentAccession,
        ".project.h5ad",
        sep = ""
    )

    # Create full URL to download R data from.
    fullUrl <- paste(
        urlBase,
        experimentAccession,
        atlasAnnDataFile,
        sep = "/"
    )

    message(
        paste(
            "Downloading Single Cell Atlas experiment annData into a temp path from:\n",
            fullUrl
        )
    )

    # Create connection object for downloading data.
    connection <- url( fullUrl )

    # Try download, catching any errors
    #loadResult <- try( load( connection ), silent = TRUE )

    # Download the H5AD file to a temporary location
    tempFile <- tempfile(fileext = ".h5ad")
    download.file(fullUrl, tempFile, mode = "wb")


    # Load the AnnData object (pkgs 'rhdf5' 'HDF5Array' required)
    adata <- readH5AD( tempFile, reader="R" )
    # Reads a H5AD file and returns a SingleCellExperiment object.

    loadResult <- try( readH5AD( tempFile, reader="R" ), silent = TRUE )


    # Quit if we got an error.
    if( class( loadResult ) == "try-error" ) {

        msg <- geterrmessage()

        warning(
            paste(
                paste(
                    "Error encountered while trying to download anndata file for",
                    experimentAccession,
                    ":"
                ),
                msg,
                paste(
                    "There may not currently be a Single-Cell Expression Atlas experiment anndata file available for ",
                    experimentAccession,
                    ".\nPlease try again later, check the website at ",
                    paste(
                        "http://www.ebi.ac.uk/gxa/sc/experiments/",
                        experimentAccession,
                        sep = ""
                    ),
                    ", or contact us at https://www.ebi.ac.uk/about/contact/support/gxasc",
                    sep = ""
                ),
                sep = "\n"
            )
        )

        return( )
    }

    # Close the connection.
    close( connection )


    # Make sure matrix X exists before trying to return it.
    getResult <- try( assays(adata)$X )

    if( class( getResult ) == "try-error" ) {

        stop(
            "ERROR - Download appeared successful but no assay X was found."
        )
    }

    # If we're still here, things must have worked ok.
    message(
        paste(
            "Successfully downloaded assay X object for",
            experimentAccession
        )
    )

    # Return the experiment summary.
    X <- assays(adata)$X 

    return( X )
    # we need to decide if we a returning here a matrix, of the complete SingleCellExperiment object
    # or an object from "anndata" package
}

# iy
# heatmapSCAtlasExperiment ( experimentAccession , genes=NULL )
#  input here would be output of 'getAtlasSCExperiment', aka a SingleCellExperiment object
#  The goal is to plot a Heatmap, from user-defined genes

# dimreadSCAtlasExperiment ( experimentAccession , dim='umap' ) {  }
