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

    max_attempts <- 5
    delay_seconds <- 3

    attempt <- 1
    loadResult <- NULL

    while (attempt <= max_attempts) {
        # Create connection object for downloading data.
        connection <- url(fullUrl)

        # Try download, catching any errors
        loadResult <- try(load(connection), silent = TRUE)

        # Close the connection after each attempt
        close(connection)

        # Check if download was successful
        if (class(loadResult) != "try-error") {
            break
        }

        message(paste("Attempt", attempt, "failed. Retrying in", delay_seconds, "seconds..."))
        Sys.sleep(delay_seconds)

        attempt <- attempt + 1
    }


    # If all attempts failed, display an error message and return
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
    xmlContent <- .downloadXMLFile(ftpUrl)
    
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

.downloadXMLFile <- function(url, attempts = 5, delay = 3) {
    for (i in 1:attempts) {
        result <- tryCatch({
            read_xml(url)
        }, error = function(e) {
            message(paste("Attempt", i, "failed:", e$message))
            NULL
        })
        if (!is.null(result)) return(result)
        Sys.sleep(delay)
    }
    stop("Failed to download XML file after multiple attempts.")
}

.downloadTabularFile <- function(fileUrl, experimentAccession, attempts = 5, delay = 3) {
    message(paste("Downloading expression file from:\n", fileUrl))
    
    for (i in 1:attempts) {
        result <- tryCatch({
            read.table(fileUrl, header = TRUE, sep = "\t", stringsAsFactors = FALSE, na.strings = "", comment.char = "#")
        }, error = function(e) {
            message(paste("Attempt", i, "failed:", e$message))
            NULL
        })
        
        if (!is.null(result)) {
            return(result)  # Return successfully downloaded data
        }
        
        Sys.sleep(delay)  # Wait before retrying
    }
    
    warning(paste("Failed to download or read file for", experimentAccession, "after", attempts, "attempts."))
    return(NULL)  # Return NULL if all attempts fail
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

    # Base URL for downloading data.
    urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

    # Define filenames for configuration files.
    configFile <- paste(experimentAccession, "-configuration.xml", sep = "")
        
    # XML URL for downloading config data.
    xmlUrl <- paste(urlBase, experimentAccession, configFile, sep = "/")

    if (normalisation %in% c("tpm", "fpkm")) {
        # Ensure the normalisation is valid.
        fileType <- match.arg(normalisation, c("tpm", "fpkm"))

        # Define filenames for TPM and FPKM files.
        expressionFile <- paste(experimentAccession, "-", fileType, "s.tsv", sep = "")

        # Create full URLs for TPM and FPKM files.
        expressionUrl <- paste(urlBase, experimentAccession, expressionFile, sep = "/")

        # Initialise a list to hold the results.
        results <- list()

        # Download and read TPM file if requested or in "both" mode.
        results <- .downloadTabularFile(expressionUrl, experimentAccession)

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



        # Read the XML file directly from the FTP URL
        message("Downloading XML file from FTP...")
        xmlContent <- .downloadXMLFile(xmlUrl)

        # Extract the assay_group id and label attributes
        assay_groups <- xml_find_all(xmlContent, ".//assay_group")

        # Create a named vector with id as names and label as values
        assay_vector <- setNames(xml_attr(assay_groups, "label"), xml_attr(assay_groups, "id"))

        # Replace the columns names in results, excluding 'GeneID' and 'Gene.Name' columns
        colnames(results)[-c(1, 2)] <- assay_vector[colnames(results)[-c(1, 2)]]

        return(results)

    } else if (normalisation %in% c("cpm")) {
        # Define filenames for TPM and FPKM files.

        experiment_summary <- getAtlasData( c(experimentAccession))[[experimentAccession]]

        # Get counts.
        counts_matrix <- assay(experiment_summary$rnaseq, "counts")

        # Calculates CPM from counts using EdgeR.
        cpm_matrix <- edgeR::cpm(counts_matrix)

        # Converts to DataFrame to match TPM and FPKM output.
        cpm_df <- as.data.frame(cpm_matrix)

        # Adds two columns GeneID and Gene.Name to DataFrame to match TPM and FPKM output format.
        # cpm_df <- cbind(GeneID = rownames(cpm_df), Gene.Name = NA, cpm_df)

        xmlContent <- .downloadXMLFile(xmlUrl)

        assay_groups <- xml_find_all(xmlContent, ".//assay_group")

        subgroups_list <- list()

        # Loop through each assay group and extract label and assays
        for (group in assay_groups) {
            # Get the label of the assay group
            label <- xml_attr(group, "label")

            # Get all assay values in the group
            assays <- xml_text(xml_find_all(group, "assay"))

            # Add to the list
            subgroups_list[[label]] <- assays
        }

        cpm_mean_df <- data.frame(row.names = rownames(cpm_df))

        # Loop through subgroups and calculate the row-wise mean
        for (group_name in names(subgroups_list)) {
            group_columns <- subgroups_list[[group_name]]

            # Ensure columns exist in 'results' before proceeding
            if(all(group_columns %in% colnames(cpm_df))) {
                # Calculate the row-wise mean for the subgroup columns and assign it to df
                cpm_mean_df[[group_name]] <- rowMeans(cpm_df[, group_columns, drop = FALSE])
            } else {
                print(paste("Error: Columns not found for", group_name))
            }
        }

        cpm_mean_df <- cbind(GeneID = rownames(cpm_mean_df), Gene.Name = NA, cpm_mean_df)

        # assayData(ex[["A-AFFY-44"]])$exprs 
        # this is microarray so not applied here

        return(cpm_mean_df)

    } else {
        stop("Defined Normalisation type ", normalisation ," is not supported, must be a tpm, fpkm or cpm.")
    }

}


heatmapAtlasExperiment <- function(df, 
                                     filename = "heatmap",
                                     save_pdf = FALSE,
                                     show_plot = TRUE,
                                     heatmap_color = "Blues",
                                     top_n = 100,
                                     show_heatmap_title = TRUE ) {  

    if (!is.data.frame(df)) stop("Input must be a dataframe.")
    
    rownames( df ) <- df[[1]]
    
    # Remove the gene id column.
    df[[1]] <- NULL
    geneIDsToGeneNames <- data.frame( id = df[[1]], stringsAsFactors=FALSE )
    
    # Add the Ensembl gene IDs as the row names.
    rownames( geneIDsToGeneNames ) <- rownames(df)

    # Identify rows where the first column (gene names) is NA
    naGeneNameIndices <- which(is.na(geneIDsToGeneNames[, 1]))

    # Replace NA values with the corresponding row names (gene IDs)
    geneIDsToGeneNames[naGeneNameIndices, 1] <- rownames(geneIDsToGeneNames)[naGeneNameIndices]

    # Now remove the Gene.Name column from the data frame.
    df[[1]] <- NULL

    # Converting the data in the data frame to numeric.
    df[] <- lapply(df, as.numeric)


    if(dim(df)[2] > 1) {
        rowVariances <- rowVars( df )
    } else {
        # can't calculate variance of one row, use the FPKM values (alhough the heatmap isn't as valuable then)
        rowVariances <- df[1]
    }

    # Get the indices of the top top_n most variable genes.
    chosenColumnIndices<- order( rowVariances, decreasing = TRUE )[1:top_n]
    topNgeneExpressions <- df[ chosenColumnIndices , ]

    # Get the gene names for the top top_n gene IDs, to use as labels for the
    # heatmape rows.
    topNgeneNames <- geneIDsToGeneNames[ chosenColumnIndices, ]

    # Scale and center the expression levels using Z-score transformation, so they
    # have mean 0 and standard deviation 1. .
    topNgeneExpressions <- t( scale( t( topNgeneExpressions )))

    # Get the assay group labels to use as the labels for the heatmap columns.
    assayGroupLabels <- colnames( topNgeneExpressions )

    # Some nice colours.
    colours <- colorRampPalette( brewer.pal( 9, heatmap_color ) )( top_n )

    imageWidth <- 8
    if( ( length( assayGroupLabels ) / 2 ) > 8 ) {
        imageWidth <- length( assayGroupLabels ) / 2
    }

    # Get the lengths of the longest assay group label
    longestLabel <- max(unlist(lapply( assayGroupLabels, function( x ) nchar( x ) )))

    # Changing image and margin height to get the column (assay
    # group) labels to fit on the page. 
    if( longestLabel / 3 > 8 ) {
        imageHeight <- ( longestLabel / 3 )
        marginHeight <- ( longestLabel / 3 )
    } else {
        imageHeight <- 8
        marginHeight <- 8
    }


    # PDF Output (Only if save_pdf = TRUE)
    if (save_pdf) {
        heatmap_file <- ifelse(grepl("\\.pdf$", filename, ignore.case = TRUE), filename, paste0(filename, ".pdf"))
        pdf(heatmap_file, height=8, width=8)
    }

    title <- paste("Gene Expression for top ", top_n, " Genes.", sep = "")

    # Make the heatmap.
    heatmap.2(
        as.matrix(topNgeneExpressions),
        col = colours,
        labRow = topNgeneNames,
        labCol = assayGroupLabels,
        key = FALSE,
        trace = "none",
        cexRow = 0.4,
        cexCol = 0.7, # hardcoding for now, may need to make this dynamic but requires thinking about.
        cex.main = 0.6,
        margins = c( marginHeight, 6 ),
        main = ifelse(show_heatmap_title, title, "")
    )

    # Close PDF if it was opened
    if (save_pdf) invisible(dev.off())

    # Show plot on screen if requested
    if (show_plot) print("Plotting on screen") 

}



getAnalysticsDifferentialAtlasExpression <- function(experimentAccession) {

    # Ensure the experiment accession is in the correct format.
    if (!.isValidExperimentAccession(experimentAccession)) {
        stop("Experiment accession not valid. Cannot continue.")
    }

    expType <- .getExperimentType(experimentAccession)

    if (expType %in% c("rnaseq_mrna_differential" )){ # add prot differential and microarray later
        message(paste(experimentAccession," is ", expType , ", will continue downloading data"))
    } else {
        stop("Experiment type not contain normalised expression data.")
    }

    # Base URL for downloading data.
    urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"

    # Define filenames for configuration files.
    configFile <- paste(experimentAccession, "-configuration.xml", sep = "")
        
    # XML URL for downloading config data.
    xmlUrl <- paste(urlBase, experimentAccession, configFile, sep = "/")

    fileType <- "analytics.tsv"

    # Define filenames for analytics file.
    analyticsFile <- paste(experimentAccession, "-", fileType, sep = "")

    # Create full URLs for analytics file.
    analyticsUrl <- paste(urlBase, experimentAccession, analyticsFile, sep = "/")

    # Initialise a list to hold the results.
    results <- list()

    # Download and read analytics file.
    results <- .downloadTabularFile(analyticsUrl, experimentAccession)

    # Remove NULL results for files that could not be downloaded.
    results <- Filter(Negate(is.null), results)
    
    xmlContent <- .downloadXMLFile(xmlUrl)

    contrasts <- xml_find_all(xmlContent, ".//contrasts")

    # # Extract contrast IDs and names
    contrast_ids <- xmlContent %>% 
        xml_find_all("//contrast") %>% 
        xml_attr("id")

    contrast_names <- xmlContent %>% 
        xml_find_all("//contrast/name") %>% 
        xml_text()

    # Combine into a data frame
    contrast_df <- data.frame(id = contrast_ids, name = contrast_names, stringsAsFactors = FALSE)

    # Replace the columns names in results, excluding 'GeneID' and 'Gene.Name' columns
    colnames(results) <- sapply(colnames(results), function(col) {
        for (i in seq_len(nrow(contrast_df))) {
            if (grepl(contrast_df$id[i], col)) {
            col <- sub(contrast_df$id[i], contrast_df$name[i], col)
            }
        }
        return(col)
        })

    # Check if any files were successfully downloaded.
    if (length(results) == 0) {
        stop("ERROR - Could not download the requested data. Please check the experiment type is RNAseq differential or try again later.")
    }

    return(results)

}

.sanitize_filename <- function(filename) {
  # Replace invalid characters with an underscore
  sanitized <- gsub("[<>:\"/\\\\|?*]", "_", filename)
  
  # Optionally, trim leading and trailing whitespace
  sanitized <- trimws(sanitized)
  
  # Return sanitized filename
  return(sanitized)
}

volcanoDifferentialAtlasExperiment <- function(df, 
                                     filename_prefix = "volcano-plot",
                                     save_pdf = FALSE,
                                     show_plot = TRUE,
                                     low_fc_colour = "gray",
                                     high_fc_colour = "blue",
                                     cutoff = 1,
                                     show_volcanoplot_title = TRUE ) {  

    # replace NA string with <NA> and remove rows with <NA>, reset rownames
    df[df == "NA"] <- NA
    clean_df <- stats::na.omit(df)
    rownames(clean_df) <- NULL

    for (i in seq(3, ncol(clean_df) - 1, by = 2)) {
        pval_col <- names(clean_df)[i]
        logFC_col <- names(clean_df)[i + 1]

        clean_df[[pval_col]] <- as.numeric(clean_df[[pval_col]])

        # Add a new column as low and high foldchange for colour
        clean_df$color <- ifelse(abs(clean_df[[logFC_col]]) < cutoff, "low", "high")

        title <- paste("Volcano Plot:", logFC_col, "vs", pval_col)
        # Create the plot
        p <- ggplot(clean_df, aes(x = .data[[logFC_col]], y = -log10(.data[[pval_col]]), color = color)) +
            geom_point(alpha = 0.6) +
            geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
            geom_vline(xintercept = c(-cutoff, cutoff), linetype = "dashed", color = high_fc_colour) +
            scale_color_manual(values = c("low" = low_fc_colour, "high" = high_fc_colour), guide = "none") +
            labs(
            title = ifelse(show_volcanoplot_title, title, ""),
            x = logFC_col,
            y = "-log10(p-value)"
            ) +
            theme_minimal()

        # PDF Output (Only if save_pdf = TRUE)
        if (save_pdf) {
            # Save each plot as a separate image
            plot_name <- gsub(".log2foldchange", "", logFC_col)
            filename <- .sanitize_filename(paste0(filename_prefix, "_", plot_name, ".png"))
            message(filename)
            ggsave(filename, plot = p, width = 8, height = 6, bg = "white")
        }

        # Optional: Print the plot in the R console
        if (show_plot) print(p)
    }
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

.ebiAPIsearch <- function (resource, query, secondaryFilter = NULL) {
    # Search the EBI RESTful Web Services API for experiment accessions
    # We will aim to return an output of similar format to the searchAtlasExperiments function:
    #   DataFrame with columns (Accession, Species, Type, Title), all <character>
    # For that we will query Atlas , SCE atlas API

    if (resource == "scxa") {
        domain_source <- "filter=domain_source:(sc-experiments)"
        atlas_api <- "https://www.ebi.ac.uk/gxa/sc/json/experiments"
        atlas_exp_type <- "technologyType"
        # performs additional ontology-based search
        cell_type_wheel <- "https://www.ebi.ac.uk/gxa/sc/json/cell-type-wheel/"
    } else if (resource == "gxa") {
        domain_source <- "filter=domain_source:(atlas-experiments)"
        atlas_api <- "https://www.ebi.ac.uk/gxa/json/experiments"
        atlas_exp_type <- "experimentType"
    } else {
        stop("Resource not supported")
    }

    # Quit if we don't have any query search terms
    if( missing( query ) ) {
        stop( "Please provide at least one query search term." )
    }

    # If we've got something other than a character vector of query, also quit.
    if( typeof( query ) != "character" ) {
        stop( "Please provide query search term(s) as a character vector." )
    }

    # Make 'query' URL friendly (e.g. replace spaces with %20
    query <- sapply( query, URLencode )

    # If we weren't passed a secondaryFilter, log this.
    if( missing( secondaryFilter ) ) {

        message( "No secondaryFilter was provided." )

    } else if( typeof( secondaryFilter ) != "character" ) {

        # Quit if the secondaryFilter isn't a character vector.
        stop( "Please provide secondaryFilter as a character vector." )

    } else if( length( secondaryFilter ) > 1 ) {

        # Only allow one secondaryFilter to be specified.
        stop( "More than one secondaryFilter found. You may only specify one secondaryFilter at a time." )
    }

    # EBI RESTful Web Services API base URL
    ebiAPIbase <- "https://www.ebi.ac.uk/ebisearch/ws/rest/geneExpression?query="

   # Construct the query URL
    queryURL <- paste(
        ebiAPIbase,
        paste( query, collapse = "" ),
        "&size=1000",  
        paste("&", domain_source, sep = ""),
        "&fields=acc",
        sep = ""
    )
    if( ! missing( secondaryFilter ) ) {
        queryURL <- paste(
            queryURL,
            "&filter=",
            paste( secondaryFilter, collapse = "" ),
            sep = ""
        )
    }
    # Log search is beginning.
    message( "Searching for experiments matching your query ..." )

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
    numExps <- as.numeric( parsedJSON$hitCount )



    # search cell type wheel
    if (resource == "scxa") {
        # first search with the 'query' and retrieve accessions
        query_ctw1 <- paste( cell_type_wheel,query, sep = "" )
        response_ctw1 <- GET( query_ctw1 )
         # Make sure the HTTP request worked.
        if( status_code( response_ctw1 ) != 200 ) {
            stop( paste( "Error running query to cell-type wheel, Received HTTP error code", status_code( response_ctw1 )) )
        } else {
            message( "Query to cell-type wheel successful." )
        }
        parsedJSON_ctw1 <- fromJSON( txt = query_ctw1 )
        
        # second, search with the 'secondaryFilter' and retrieve accessions
        if( missing(secondaryFilter) || is.null(secondaryFilter) ) {
            additionalExps <- unique(unlist(parsedJSON_ctw1$experimentAccessions))    

        } else {
            query_ctw2 <- paste( cell_type_wheel,secondaryFilter, sep = "" )
            response_ctw2 <- GET( query_ctw2 )
            if( status_code( response_ctw2 ) != 200 ) { 
                stop( paste( "Error running query to cell-type wheel, Received HTTP error code", status_code( response_ctw2 )) )
            } else { 
                message( "Query to cell-type wheel successful." )
            }
            parsedJSON_ctw2 <- fromJSON( txt = query_ctw2 )
            additionalExps <- intersect( unique(unlist(parsedJSON_ctw1$experimentAccessions))  ,   unique(unlist(parsedJSON_ctw2$experimentAccessions))   )
        }
        additionalNumExps <- length( additionalExps )

    }

    if (resource == "gxa") {
        additionalExps <- 0
        additionalNumExps <- length( additionalExps )
    }



    # If there were no results, quit here.
    if( numExps == 0 && additionalNumExps == 0 ) {
        return( message( "No results found. Cannot continue." ) )

    } else {
        message( paste( "Found", numExps+additionalNumExps, "experiments matching your query." ) )
    }

    ebiResult <- DataFrame(
        Accession = unique( as.character(parsedJSON$entries$id) , additionalExps ),
        Species   = "NA",
        Type      = "NA",
        Title     = "NA"
    )

    # now we need to query Atlas API to get the species, type and title
    atlasAPIresponse <- GET( atlas_api )

    # Make sure the HTTP request worked.
    if( status_code( atlasAPIresponse ) != 200 ) {
        stop(
            paste(
                "Error running query. Received HTTP error code from ATLAS API",
                status_code( atlasAPIresponse ),
                "from server. Please try again later. If you continue to experience problems please contact us at https://www.ebi.ac.uk/about/contact/support/gxa"
            )
        )
    } else {
        message( "Query successful." )
    }

    # Parse the JSON document
    parsedJSON <- fromJSON(content(atlasAPIresponse, "text"))

    experiment_data <- DataFrame(
        Accession = unlist(parsedJSON$experiments[["experimentAccession"]]),
        Species = unlist(parsedJSON$experiments[["species"]]),
        Type = sapply(parsedJSON$experiments[[atlas_exp_type]], paste, collapse = ", "),
        Title = unlist(parsedJSON$experiments[["experimentDescription"]])
    )

    filtered_data <- subset(experiment_data, Accession %in% ebiResult$Accession )

    return( filtered_data )

}


searchSCAtlasExperiments <- function( query, secondaryFilter = NULL ) { 
    # wrapper function to search for Single Cell Atlas experiments
     if( missing( secondaryFilter ) ) {
        scExperiments <- .ebiAPIsearch(resource = "scxa", query)
    } else {
        scExperiments <- .ebiAPIsearch(resource = "scxa", query, secondaryFilter)
    }

    return( scExperiments )
}

searchAtlasExperiments <- function( query, secondaryFilter = NULL ) { 
    # wrapper function to search for Expression Atlas experiments
     if( missing( secondaryFilter ) ) {
        scExperiments <- .ebiAPIsearch(resource = "gxa", query)
    } else {
        scExperiments <- .ebiAPIsearch(resource = "gxa", query, secondaryFilter)
    }

    return( scExperiments )
}




.download_anndata_with_retry <- function(url, destfile, max_attempts = 5) {
  for (attempt in 1:max_attempts) {
    tryCatch(
      {
        download.file(url, destfile, mode = "wb")
        message("anndata file download successful.")
        return(TRUE)
      },
      error = function(e) {
        message(sprintf("Attempt %d/%d to download anndata failed: %s", attempt, max_attempts, e$message))
        if (attempt == max_attempts) stop("Max retries reached. Download of anndata file failed.")
        Sys.sleep(2*attempt)
      }
    )
  }
}




getAtlasSCExperiment <- function( experimentAccession ) {

    # Make sure the experiment accession is in the correct format
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
            "Downloading Single-Cell Atlas experiment annData into a temp path from:\n",
            fullUrl
        )
    )


    # Download the H5AD file to a temporary location
    tempFile <- tempfile(fileext = ".h5ad")

    on.exit( unlink(tempFile) )  # Ensures the temp file is deleted when the function exits

    .download_anndata_with_retry(fullUrl, tempFile)

    # Read a H5AD file and returns a SingleCellExperiment object
    loadResult <- try( readH5AD( file = tempFile, reader="R" , use_hdf5 = TRUE ), silent = TRUE )

    # Quit if we got an error.
    if( class( loadResult ) == "try-error" ) {

        msg <- geterrmessage()

        warning(
            paste(
                paste(
                    "Error encountered while trying to load anndata file for ",
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


    # Make sure 'normalised' expression exists in the SCE object
    getResult <- try( assays(loadResult)$normalised )

    # Here I could check other data needed by plot functions exists,
    # e.g. reducedDim(loadResult, "X_pca") for PCA plot

    if( class( getResult ) == "NULL" || class( getResult ) == "try-error" ) {

        stop(
            "ERROR - Download and load appeared successful but no assay 'normalised' matrix was found."
        )
    }

    # If we're still here, things must have worked ok.
    message(
        paste(
            "Successfully loaded assay 'normalised' from SingleCellExperiment object for ",
            experimentAccession
        )
    )

    # Return SingleCellExperiment object
    return( loadResult )

}


plotDimRedSCAtlasExperiment <- function( sceObject, dimRed, colorby ) {
    
    # Check if the provided dimRed exists in reducedDimNames
    if (!(dimRed %in% reducedDimNames(sceObject))) {
        stop(paste("Error: Dimension reduction method", dimRed, "not found in the object!"))
    }
    
    # Extract dimension reduction coordinates
    dim_coords <- reducedDim(sceObject, dimRed)

    # Check if colorby exists in colData
    if (!(colorby %in% colnames(colData(sceObject)))) {
        stop(paste("Error: Color attribute", colorby, "not found in colData!"))
    }

    # Convert to a data frame
    dimred_df <- data.frame(
        axis1 = dim_coords[, 1], 
        axis2 = dim_coords[, 2], 
        colorBy = as.character( colData(sceObject)[, colorby] )
    )

    # Plot using ggplot2
    ggplot(dimred_df, aes(x = axis1, y = axis2, color = colorBy)) +
        geom_point() +
        labs(
            title = paste("Dimensionality Reduction Plot:", dimRed), 
            x = "Component 1", 
            y = "Component 2",
            color = colorby
        ) +
        theme_minimal()
}



# iy
# heatmapSCAtlasExperiment ( experimentAccession , genes=NULL )
#  input here would be output of 'getAtlasSCExperiment', aka a SingleCellExperiment object
#  The goal is to plot a Heatmap, for user-defined genes
