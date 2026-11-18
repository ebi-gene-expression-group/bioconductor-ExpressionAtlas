# ----- Bulk  ----

# getAtlasExperiment
#   - Download and return the SimpleList object representing a single
#     Expression Atlas experiment.
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
  
  max_attempts <- 3
  delay_seconds <- 3
  
  
  default_time <- getOption("timeout")
  on.exit(options(timeout = default_time), add = TRUE)
  
  attempt <- 1
  loadResult <- NULL
  
  
  for (attempt in seq_len(max_attempts)){
    
    # Try download first time for 1 minute.
    if(attempt == 1){
      message('Downloading....')
    }
    
    # Retry second time for 1 minute. 
    if(attempt == 2){
      message('Trying again...')
    }
    
    # Retry third time for 10 minutes. 
    if(attempt == 3){
      message('The file may be too large.')
      message('Trying last time for 10 minutes....')
      options(timeout = if (attempt == 3) 600 else default_time)
    }
    
    # Create connection object for downloading data.
    connection <- url(fullUrl)
    
    # Try download, catching any errors.
    loadResult <- try(load(connection), silent = TRUE)
    
    # Close the connection after each attempt.
    close(connection)
    
    # Check if download was successful.
    if (!inherits(loadResult, "try-error")) {
      break
    } 
    
    if(attempt < max_attempts){
      message(paste("Attempt", attempt, "failed. Retrying in", delay_seconds, "seconds..."))
      Sys.sleep(delay_seconds)
    }
  }
  
  
  # If all attempts failed, error message displayed. 
  if (inherits(loadResult, "try-error")) {
    
    
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
          ".\nPlease try to find a better connection or try again later, check the website at ",
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
  
  if (inherits(getResult, "try-error")) {
    
    stop(
      "ERROR - Download appeared successful but no experiment summary object was found."
    )
  }
  
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
  
  # Stop if no valid accessions remain after filtering.
  if( length( experimentAccessions ) == 0 ) {
    stop( "None of the accessions passed are valid ArrayExpress/BioStudies accessions. Cannot continue." )
  }
  
  # Named the input vector so the resulting list is keyed by accession.
  names( experimentAccessions ) <- experimentAccessions
  
  # Download the experiment data for each accession and collect it in a SimpleList.
  experimentSummaryList <- SimpleList(
    
    lapply( experimentAccessions, function( experimentAccession ) {
      
      experimentSummary <- getAtlasExperiment( experimentAccession )
    }
    ) )
  
  # Keep only accessions that returned an experiment summary.
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
        "\" does not look like an ArrayExpress/BioStudies experiment accession. Please check.",
        sep=""
      )
    )
    
    return( FALSE )
    
  } else {
    
    return( TRUE )
  }
}

# .getExperimentType 
# - Internal helper to get Experiment Type from accession. 
.getExperimentType <- function(experimentAccession) {
  
  # Base URL for downloading data.
  urlBase <- "ftp://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments"
  
  # Define filenames for configuration files.
  configFile <- paste(experimentAccession, "-configuration.xml", sep = "")
  
  # Define ftpUrl for configuration file.
  ftpUrl <- paste(urlBase, experimentAccession, configFile, sep = "/")
  
  # Read the XML file directly from the FTP URL.
  message("Downloading XML file from FTP...")
  xmlContent <- .downloadXMLFile(ftpUrl)
  
  # Extract the <configuration> node and the experimentType attribute.
  configurationNode <- xml_find_first(xmlContent, "/configuration")
  if (is.na(configurationNode)) {
    stop("No <configuration> node found in the XML file.")
  }
  
  # Get the value of the experimentType attribute.
  experimentType <- xml_attr(configurationNode, "experimentType")
  if (is.na(experimentType)) {
    stop("The 'experimentType' attribute is missing in the <configuration> node.")
  }
  
  return(experimentType)
}

# .downloadXMLFile 
# - Internal helper to download XML file. 

.downloadXMLFile <- function(url, attempts = 7, delay = 3) {
  for (i in 1:attempts) {
    result <- tryCatch({
      read_xml(url)
    }, error = function(e) {
      message(paste("Attempt", i, "failed:", e$message))
      NULL
    })
    if (!is.null(result)) return(result) # Return downloaded data 
    Sys.sleep(delay) # Wait before retrying 
  }
  stop("Failed to download XML file after multiple attempts.")
}

# .downloaTabularFile
# - Internal helper to download TSV file. 
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

# .loadmarkers 
# - Internal helper to download markers TSV file from FTP. 
.loadmarkers <- function(experimentAccession, attempts = 5, delay = 3){
  
  experiment_markers <- paste0(experimentAccession, "-tpms-markers.tsv")
  
  # URL to download the file. 
  url_for_markers <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/experiments/",
    experimentAccession, "/",
    experiment_markers
  )
  
  message("Downloading markers file:\n", url_for_markers)
  
  for (i in seq_len(attempts)) {
    markers <- tryCatch({
      read.delim(
        url_for_markers,
        na.strings = c("NULL", "NA", ""),
        stringsAsFactors = FALSE
      )
    }, error = function(e){
      message("Attempt", i, " failed: ", e$message)
      NULL
    })
    if (!is.null(markers)) {
      return(markers) # Return if file successfully loaded.
    }
    if(i < attempts) {
      Sys.sleep(delay) # Wait before retrying 
    }
  }
  
  warning(
    "Failed to download markers for ",
    experimentAccession, " after ", attempts, " attempts."
  )
  return(NULL)
  
}

# getNormalisedAtlasExpression
# - Download normalized data  for single Expression Atlas experiment 
#   and return data frame with marker genes rank. 
getNormalisedAtlasExpression <- function(experimentAccession, normalisation = "tpm") {
  
  # Ensure the experiment accession is in the correct format.
  if (!.isValidExperimentAccession(experimentAccession)) {
    stop("Experiment accession not valid. Cannot continue.")
  }
  
  # Ensure the experiment type is correct. 
  expType <- .getExperimentType(experimentAccession)
  
  if (expType == "rnaseq_mrna_differential") {
    stop(
      paste(
        "This is a differential RNA-seq experiment.",
        "Use getAnalyticsDifferentialAtlasExpression() instead.",
        sep = "\n"
      )
    )
  }
  
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
  
  # Download markers annotations: 
  markers <- .loadmarkers(experimentAccession)
  
  # Getting ranked gene-assay annotations for merging. 
  ranked_markers <- markers[
    !is.na(markers$marker_gene_rank),
    c("gene_id", "assay", "marker_gene_rank"),
  ]
  
  if (normalisation %in% c("tpm")) {
    # Ensure the normalisation is valid.
    fileType <- match.arg(normalisation, c("tpm"))
    
    # Define filenames for TPM files.
    expressionFile <- paste(experimentAccession, "-", fileType, "s.tsv", sep = "")
    
    # Create full URLs for TPM files.
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
    
    # TPMs can have multiple values in expression comma-separated columns (min, lower quartile, mean, higher quartile and max).
    get_mean_value <- function(col) {
      # Split by commas and return the third value.
      sapply(strsplit(col, ","), function(x) x[3])
    }
    
    # The first two columns are GeneID and Gene.Name, fetching 3rd column as it contains mean values. 
    results[ ,3:ncol(results)] <- lapply(results[ ,3:ncol(results)], get_mean_value)
    
    
    
    # Read the XML file directly from the FTP URL.
    message("Downloading XML file from FTP...")
    xmlContent <- .downloadXMLFile(xmlUrl)
    
    # Extract the assay_group id and label attributes.
    assay_groups <- xml_find_all(xmlContent, ".//assay_group")
    
    # Create a named vector with id as names and label as values.
    assay_vector <- setNames(xml_attr(assay_groups, "label"), xml_attr(assay_groups, "id"))
    
    # Replace the columns names in results, excluding 'GeneID' and 'Gene.Name' columns.
    colnames(results)[-c(1, 2)] <- assay_vector[colnames(results)[-c(1, 2)]]
    
    # Adding markers annotations by matching Gene ID. 
    results <- merge(results,
                     ranked_markers,
                     by.x = "GeneID", 
                     by.y = "gene_id",
                     all.x = TRUE)
    
    # Rename 'assay' to 'tissue' column in data frame.
    names(results)[names(results) == "assay"] <- "tissue"
    
    # Replace NA in gene names with GeneID.
    missing_names <- is.na(results$Gene.Name) |
      results$Gene.Name == ""
    
    results$Gene.Name[missing_names] <-
      results$GeneID[missing_names]
    
    # Set the normalisation label to TPM for the heatmap. 
    attr(results, "normalisation") <- "TPM"
    
    return(results)
    
  } else if (normalisation %in% c("cpm")) {
    
    # Download the experiment summary and extract it by accession.
    experiment_summary <- getAtlasData( c(experimentAccession))[[experimentAccession]]
    
    # Get counts.
    counts_matrix <- assay(experiment_summary$rnaseq, "counts")
    
    # Calculates CPM from counts using EdgeR.
    cpm_matrix <- edgeR::cpm(counts_matrix)
    
    # Converts to DataFrame to match TPM output.
    cpm_df <- as.data.frame(cpm_matrix)
    
    xmlContent <- .downloadXMLFile(xmlUrl)
    
    assay_groups <- xml_find_all(xmlContent, ".//assay_group")
    
    subgroups_list <- list()
    
    # Loop through each assay group and extract label and assays.
    for (group in assay_groups) {
      # Get the label of the assay group.
      label <- xml_attr(group, "label")
      
      # Get all assay values in the group.
      assays <- xml_text(xml_find_all(group, "assay"))
      
      # Add to the list.
      subgroups_list[[label]] <- assays
    }
      
    cpm_mean_df <- data.frame(row.names = rownames(cpm_df))
    
    # Calculate mean CPM for each assay group.
    for (group_name in names(subgroups_list)) {
      group_columns <- subgroups_list[[group_name]]
      
      if (all(group_columns %in% colnames(cpm_df))) {
        cpm_mean_df[[group_name]] <- rowMeans(
          cpm_df[, group_columns, drop = FALSE]
        )
      } else {
        print(paste("Error: Columns not found for", group_name))
      }
    }
    
    # Use GeneID as the default Gene name. 
    cpm_mean_df <- cbind(
      GeneID = rownames(cpm_mean_df),
      Gene.Name = rownames(cpm_mean_df),
      cpm_mean_df
    )
    
    # Download TPM data to look up Gene names. 
    tpm_url <- paste(
      urlBase,
      experimentAccession,
      paste0(experimentAccession, "-tpms.tsv"),
      sep = "/"
    )
    
    tpm_data <- .downloadTabularFile(tpm_url, experimentAccession)
    
    if (!is.null(tpm_data) &&
        all(c("GeneID", "Gene.Name") %in% names(tpm_data))) {
      
      # Keep the first entry for each GeneID.
      gene_lookup <- tpm_data[
        !duplicated(tpm_data$GeneID),
        c("GeneID", "Gene.Name")
      ]
      
      # Match gene names to CPM rows.
      match_index <- match(cpm_mean_df$GeneID, gene_lookup$GeneID)
      gene_names <- as.character(gene_lookup$Gene.Name[match_index])
      
      # Replace GeneID only where a gene name is available.
      has_name <- !is.na(gene_names) & trimws(gene_names) != ""
      cpm_mean_df$Gene.Name[has_name] <- gene_names[has_name]
      
    } else {
      warning("TPM gene names unavailable; using GeneID as Gene.Name.")
    }
    
    # Merging with marker file.  
    cpm_mean_df <- merge(
      cpm_mean_df,
      ranked_markers, 
      by.x = "GeneID",
      by.y = "gene_id",
      all.x = TRUE 
    )
    
    # Renaming assay column to tissue.  
    names(cpm_mean_df)[names(cpm_mean_df) == "assay"] <- "tissue"
    
    # Replacing NA Gene names with GeneID. 
    missing_names_cpm <- is.na(cpm_mean_df$Gene.Name) |
      cpm_mean_df$Gene.Name == ""
    
    cpm_mean_df$Gene.Name[missing_names_cpm] <-
      cpm_mean_df$GeneID[missing_names_cpm]
    
    # Set the normalisation label to CPM for the heatmap.
    attr(cpm_mean_df, "normalisation") <- "CPM"
    
    return(cpm_mean_df)
    
  } else {
    stop(
      "Defined Normalisation type ", normalisation,
      " is not supported, must be a tpm or cpm."
    )
  } 
  
} 

# heatmapAtlasExperiment 
# - function that plots heatmap using normalisation data 
#   from getNormalisedAtlasExpression function. 

heatmapAtlasExperiment <- function(df,
                                   filename = "heatmap",
                                   save_pdf = FALSE,
                                   show_plot = TRUE,
                                   palette = 'viridis',
                                   markers = TRUE,
                                   gene_list = NULL,
                                   top_n = 100,
                                   scaled = FALSE,
                                   show_heatmap_title = TRUE,
                                   geneLabel = c("name", "id")
)
{
  
  # Get one of the supported versions for geneLabel. 
  geneLabel <- match.arg(geneLabel)
  
  # Scaled must be TRUE or FALSE. 
  if (!is.logical(scaled)) {
    stop("'scaled' parameter must be TRUE or FALSE")
  }
  
  # df must be data frame.  
  if (!is.data.frame(df)) stop("Input must be a dataframe.")
  
  # Reading normalisatio label. 
  units <- attr(df, "normalisation", exact = TRUE)
  
  # Matching genes from gene list. 
  user_genes <- df$Gene.Name %in% gene_list |
    df$GeneID %in% gene_list
  
  # Exclude genes that are not in data frame.  
  missing_genes <- setdiff(gene_list,
                           c(df$Gene.Name, df$GeneID))
  
  if(length(missing_genes) > 0) {
    message(
      "These genes are absent from the experiment and cannot be plotted: ",
      paste0(missing_genes, collapse = ", ")
    )
  }
  
  
  # If markers = TRUE, plot rank #1 genes and genes from user list. 
  if(markers) {
    
    top_markers <- !is.na(df$marker_gene_rank) &
      df$marker_gene_rank == 1
    
    df <- df[top_markers | user_genes, ]
    
  } else if(length(gene_list) > 0 ) { 
    # Keep only genes from usrr list. 
    df <- df[user_genes, ]
    
  }
  
  # Stop if no genes remain to plot.
  if(nrow(df) == 0 ) {
    stop("No genes available from requested gene_list")
  }
  
  # Order by tissue. 
  df <- df[
    order(df$tissue, df$marker_gene_rank),
  ]
  
  # Keep one row per GeneID. 
  df <- df[!duplicated(df$GeneID), ]
  
  # Make Gene Name - gene labels. 
  gene_labels <- as.character(df$Gene.Name)
  
  rownames(df) <- df$GeneID
  
  # Remove annotations, only keep expression for heatmap.  
  df <- df[
    ,
    !names(df) %in% c(
      "GeneID", "Gene.Name", "tissue", "marker_gene_rank"
    )
  ]
  
  
  # Make columns numeric.  
  df[] <- lapply(df, function(x) as.numeric(as.character(x)))
  
  # Covert numeric to expression matrix. 
  expression_matrix <- as.matrix(df)
  
  # Select all rows. 
  chosen_rows <- seq_len(nrow(expression_matrix))
  
  # If no markers and no gene list, take top_n genes. 
  if (!markers && length(gene_list) == 0) {
    # Rank genes by variance across columns or by expression for a single column.
    scores <- if (ncol(expression_matrix) > 1) {
      genefilter::rowVars(expression_matrix)
    } else {
      expression_matrix[, 1]
    }
    
    # Exclude missing scores and keep the rows with the highest scores.
    chosen_rows <- utils::head(order(scores, decreasing = TRUE, na.last = NA), top_n)
  }
  
  # If no genes selected, stop. 
  if (length(chosen_rows) == 0) {
    stop("No genes have usable expression scores.")
  }
  
  topNgeneExpressions <- expression_matrix[chosen_rows, ]
  # Sort columns by name. 
  topNgeneExpressions <- topNgeneExpressions[
    , order(colnames(topNgeneExpressions))
  ]
  
  # Get GeneIDs from row names. 
  gene_ids <- rownames(topNgeneExpressions)
  
  # Use GeneID or Gene name as row labels for heatmap.
  if(geneLabel == 'id') {
    topNgeneNames <- gene_ids
  } else {
    topNgeneNames <- gene_labels[chosen_rows] 
  }
  
  missing_ids <- is.na(topNgeneNames) | topNgeneNames == ""
  topNgeneNames[missing_ids] <- gene_ids[missing_ids]
  
  # If several genes have same Gene name, append GeneID to them. 
  shared_name <- duplicated(topNgeneNames) |
    duplicated(topNgeneNames, fromLast = TRUE)
  
  topNgeneNames[shared_name] <- paste(
    topNgeneNames[shared_name],
    gene_ids[shared_name],
    sep = "_"
  )
  
  
  # Scale and center the expression levels using Z-score transformation, so they
  # have mean 0 and standard deviation 1.
  
  if (isTRUE(scaled)) {
    topNgeneExpressions <- t( scale( t( topNgeneExpressions )))
  } else {
    topNgeneExpressions <- log(topNgeneExpressions + 1)
  }
  
  
  # Get the assay group labels to use as the labels for the heatmap columns.
  assayGroupLabels <- colnames( topNgeneExpressions )
  
  min_val <- min(topNgeneExpressions, na.rm = TRUE)
  max_val <- max(topNgeneExpressions, na.rm = TRUE)
  mid_val <- (min_val + max_val) / 2
  
  # Viridis colours.
  colours <- circlize::colorRamp2(
    breaks = c(min_val, mid_val, max_val),
    colors = viridis(3, option = palette)
  )
  
  # Title done according to parameters chosen for heatmap. 
  title <- paste("Gene Expression for top ", top_n, " Genes", sep = "")
  
  if(markers && length(gene_list) > 0) {
    title <- paste("Top ranked genes with gene list")
  } else if (markers && length(gene_list) == 0) {
    title <- ("Top ranked genes ")
  } else if (!is.na(markers) && length(gene_list) > 0) {
    title <- ("Gene list ")
  }
  
  # If scaled = TRUE, z-score displayed on heatmap bar. 
  title_small <- if(isTRUE(scaled)) {
    
    paste("Gene Expression\n", units, "z-score")
    
  } else { # else log + 1
    paste("Gene Expression\nlog(", units, "+ 1)")
  }
  
  
  heatmap_vis <- ComplexHeatmap::Heatmap(
    topNgeneExpressions,
    row_labels = topNgeneNames,
    name = title_small,
    col = colours, 
    use_raster = TRUE,
    raster_quality = 5,
    cluster_rows = !markers && nrow(topNgeneExpressions) >1,  
    cluster_columns = !markers && ncol(topNgeneExpressions) >1,
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = 0.4 * 12),       
    column_names_gp = gpar(fontsize = 0.7 * 12),    
    column_title = ifelse(show_heatmap_title, title, ""),
    column_title_gp = gpar(fontsize = 0.6 * 12),
    show_heatmap_legend = TRUE,
    row_names_max_width = unit(6, "cm")
  )
  
  if (show_plot) {
    ComplexHeatmap::draw(heatmap_vis)
  }
  
  # PDF save (only if user choose save_pdf = TRUE)
  if (save_pdf) {
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
    
    pdf(filename, width = 8, height = 8)
    ComplexHeatmap::draw(heatmap_vis)
    dev.off()
  }
  
  return(heatmap_vis)
  
}

# getAnalyticsDifferentialAtlasExpression 
# - Downloads and returns differential expression analytics data for an 
#   Expression Atlas experiment. 
getAnalyticsDifferentialAtlasExpression <- function(experimentAccession) {
  
  # Ensure the experiment accession is in the correct format.
  if (!.isValidExperimentAccession(experimentAccession)) {
    stop("Experiment accession not valid. Cannot continue.")
  }
  
  # Ensure the experiment type is correct. 
  expType <- .getExperimentType(experimentAccession)
  
  if (expType == "rnaseq_mrna_baseline") {
    stop(
      paste(
        "This is a baseline RNA-seq experiment.",
        "Use getNormalisedAtlasExpression() instead.",
        sep = "\n"
      )
    )
  }
  
  if (expType %in% c("rnaseq_mrna_differential" )){ 
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
  
  # Extract contrast IDs and names.
  contrast_ids <- xmlContent %>% 
    xml_find_all("//contrast") %>% 
    xml_attr("id")
  
  contrast_names <- xmlContent %>% 
    xml_find_all("//contrast/name") %>% 
    xml_text()
  
  # Combine into a data frame.
  contrast_df <- data.frame(id = contrast_ids, name = contrast_names, stringsAsFactors = FALSE)
  
  # Replace the columns names in results, excluding 'GeneID' and 'Gene.Name' columns.
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

# .sanitize-filename
# Internal helper that replace invalid characters with underscores 
# and trims whitespaces. 
.sanitize_filename <- function(filename) {
  # Replace invalid characters with an underscore
  sanitized <- gsub("[<>:\"/\\\\|?*]", "_", filename)
  
  # Optionally, trim leading and trailing whitespace
  sanitized <- trimws(sanitized)
  
  # Return sanitized filename
  return(sanitized)
}

# volcanoDifferentialAtlasExperiment 
# - Function that plots volacno plot from 
#   diffeential expression data from getAnalyticsDifferentialAtlasExpression function. 
volcanoDifferentialAtlasExperiment <- function(df, 
                                               filename_prefix = "volcano-plot",
                                               save_png = FALSE,
                                               show_plot = TRUE,
                                               low_fc_colour = "gray",
                                               high_fc_colour = "blue",
                                               cutoff = 1,
                                               show_volcanoplot_title = TRUE ) {  
  
  # replace NA string with <NA> and remove rows with <NA>, reset rownames.
  df[df == "NA"] <- NA
  clean_df <- stats::na.omit(df)
  rownames(clean_df) <- NULL
  
  for (i in seq(3, ncol(clean_df) - 1, by = 2)) {
    pval_col <- names(clean_df)[i]
    logFC_col <- names(clean_df)[i + 1]
    
    clean_df[[pval_col]] <- as.numeric(clean_df[[pval_col]])
    
    # Add a new column as low and high foldchange for colour.
    clean_df$color <- ifelse(abs(clean_df[[logFC_col]]) < cutoff, "low", "high")
    
    title <- paste(
      strwrap(paste("Volcano Plot:", logFC_col, "vs", pval_col), width = 48),
      collapse = "\n"
    )
    # Create the plot.
    p <- ggplot(clean_df, aes(x = .data[[logFC_col]], y = -log10(.data[[pval_col]]), color = .data$color)) +
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
    
    # PDF Output (Only if save_pdf = TRUE).
    if (save_png) {
      # Save each plot as a separate image.
      plot_name <- gsub(".log2foldchange", "", logFC_col)
      filename <- .sanitize_filename(paste0(filename_prefix, "_", plot_name, ".png"))
      message(filename)
      ggsave(filename, plot = p, width = 8, height = 6, bg = "white")
    }
    
    # Optional: Print the plot in the R console.
    if (show_plot) print(p)
  }
  
  return(p)
}

# ----- Search ----
# .flatten_experiment
# - Internal helper: function to extract all nested elements as a flat string
# used by .search_gxa_experiments_yaml

.flatten_experiment <- function(exp) {
  # extract top-level details
  details_list <- list(
    accession = exp$accession,
    experiment_type = exp$experiment_type,
    organism = exp$organism
  )
  
  # Extract data from assay_groups if present. 
  if (!is.null(exp$assay_groups)) {
    for (group in exp$assay_groups) {
      details_list <- append(details_list, unlist(group, use.names = FALSE))
    }
  }
  
  # Concatenate all extracted values into a single searchable string. 
  return(paste(details_list, collapse = " | "))
}

# .searcg_gxa_experiemnts_yaml 
# - Internal helper to search experiments using YAML. 

.search_gxa_experiments_yaml <- function(query_yaml, attempts = 7, delay = 3) {
  yaml_gxa <- "https://ftp.ebi.ac.uk/pub/databases/microarray/data/atlas/gxa-studies.yaml"
  
  for (i in seq_len(attempts)) {
    yaml_response <- tryCatch(
      GET(yaml_gxa, httr::timeout(60)),
      error = function(e) {
        message(paste("Attempt", i, "failed: ", e$message))
        NULL
      }
    )
    
    if (!is.null(yaml_response)) break
    if (i < attempts) Sys.sleep(delay)
  }
  
  if (is.null(yaml_response)) {
    stop("YAML failed to download after multiple attempts.")
  }
  
  # Check if the request was successful.
  if (status_code(yaml_response) == 200) {
    yaml_text <- content(yaml_response, "text", encoding = "UTF-8")
    
    yaml_data <- tryCatch(
      suppressWarnings( read_yaml(text = yaml_text)),
      error = function(e) {
        message("Error loading YAML: ", e$message)
        NULL
      }
    )
    
    # Ensure all elements remain in their original format.
    yaml_data <- lapply(yaml_data, function(x) {
      if (is.list(x)) {
        return(lapply(x, function(y) {
          if (is.character(y)) {
            return(y)  # Keep as character
          } else if (is.numeric(y)) {
            return(as.numeric(y))  # Keep as numeric
          } else {
            return(y)  # Preserve other types
          }
        }))
      } else {
        return(x)
      }
    })
    
  } else {
    stop("Failed to retrieve YAML file. HTTP status: ", status_code(yaml_response))
  }
  
  
  # Cononvert experiments into a data frame.
  experiments_df <- tidyr::tibble(
    accession = sapply( yaml_data$experiments, function(x) x$accession),
    details = sapply( yaml_data$experiments, .flatten_experiment)  # Flatten all details
  )
  
  # Define search term. 
  search_term <- as.character(query_yaml) 
  
  # search across all details
  matching_accessions <- experiments_df %>%
    dplyr::filter(grepl(search_term, .data$details, ignore.case = TRUE)) %>%
    dplyr::pull(.data$accession)
  
  return(matching_accessions)
  
}

# .ebiAPIsearch
# - Internal helper: Search the EBI RESTful Web Services API for experiment accessions. 
.ebiAPIsearch <- function (resource, query, secondaryFilter = NULL, detailed = FALSE) {

  if (!is.logical(detailed)) {
    stop("'detailed' parameter must be TRUE or FALSE")
  }
  
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
  
  # Quit if we don't have any query search terms.
  if( missing( query ) ) {
    stop( "Please provide at least one query search term." )
  }
  
  # If we've got something other than a character vector of query, also quit.
  if( typeof( query ) != "character" ) {
    stop( "Please provide query search term(s) as a character vector." )
  }
  
  # Make 'query' URL friendly (e.g. replace spaces with %20). 
  query <- sapply( query, URLencode )
  
  # If does not pass a secondaryFilter, log this.
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
  
  # Check the query works with httr.
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
  
  # Parse the JSON document.
  parsedJSON <- fromJSON( txt = queryURL )
  
  # Get the number of experiments found.
  numExps <- as.numeric( parsedJSON$hitCount )
  
  
  
  # search cell type wheel.
  if (resource == "scxa") {
    # first search with the 'query' and retrieve accessions
    query_ctw1 <- paste( cell_type_wheel,as.character(query), sep = "" )
    response_ctw1 <- GET( query_ctw1 )
    # Make sure the HTTP request worked.
    if( status_code( response_ctw1 ) != 200 ) {
      stop( paste( "Error running query to cell-type wheel, Received HTTP error code", status_code( response_ctw1 )) )
    } else {
      message( "Query to cell-type wheel successful." )
    }
    parsedJSON_ctw1 <- fromJSON( txt = query_ctw1 )
    
    # Second, search with the 'secondaryFilter' and retrieve accessions. 
    if( missing(secondaryFilter) || is.null(secondaryFilter) ) {
      additionalExps <- unique(unlist(parsedJSON_ctw1$experimentAccessions))
      
    } else {
      query_ctw2 <- paste( cell_type_wheel,as.character(secondaryFilter), sep = "" )
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
  
  # Search using YAML file. 
  if (resource == "gxa") {
    if ( detailed == FALSE ) {
      additionalExps <- 0
      additionalNumExps <- additionalExps #length( additionalExps )
    } else {
      if( missing(secondaryFilter) || is.null(secondaryFilter) ) {
        additionalExps <- .search_gxa_experiments_yaml( as.character(query) )
      } else {
        additionalExps <- intersect ( .search_gxa_experiments_yaml( as.character(query) ) , .search_gxa_experiments_yaml( as.character(secondaryFilter) ) )
      }
      additionalNumExps <- length( additionalExps )
    }
  }
  
  
  
  # If there were no results, quit here.
  if( numExps == 0 && additionalNumExps == 0 ) {
    return( message( "No results found. Cannot continue." ) )
    
  } else {
    message( paste( "Found", numExps+additionalNumExps, "experiments matching your query." ) )
  }
  
  ebiResult <- DataFrame(
    Accession = unique(c(as.character(parsedJSON$entries$id), as.character(additionalExps) )),
    Species   = "NA",
    Type      = "NA",
    Title     = "NA"
  )
  
  # Query Atlas API to get the species, type and title. 
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
  
  # Parse the JSON document.
  parsedJSON <- fromJSON(content(atlasAPIresponse, "text"))
  
  experiment_data <- DataFrame(
    Accession = unlist(parsedJSON$experiments[["experimentAccession"]]),
    Species = unlist(parsedJSON$experiments[["species"]]),
    Type = sapply(parsedJSON$experiments[[atlas_exp_type]], paste, collapse = ", "),
    Title = unlist(parsedJSON$experiments[["experimentDescription"]])
  )
  
  filtered_data <- experiment_data[experiment_data$Accession %in% ebiResult$Accession, , drop = FALSE]
  
  return( filtered_data )
  
}

# searchSCAtlasExperiement
# - Function to search for Single Cell Atlas experiments. 
searchSCAtlasExperiments <- function( query, secondaryFilter = NULL ) { 
  # Apply secondary filter if supplied. 
  if( missing( secondaryFilter ) ) {
    scExperiments <- .ebiAPIsearch(resource = "scxa", query)
  } else {
    scExperiments <- .ebiAPIsearch(resource = "scxa", query, secondaryFilter)
  }
  
  return( scExperiments )
}

# searchAtlasExperiments 
# - Function to search for Expression Atlas experiments.
searchAtlasExperiments <- function( query, secondaryFilter = NULL, detailed = FALSE ) { 
  # Apply secondary filter, if supplied. 
  if( missing( secondaryFilter ) ) {
    scExperiments <- .ebiAPIsearch(resource = "gxa", query = query, detailed = detailed)
  } else {
    scExperiments <- .ebiAPIsearch(resource = "gxa", query = query, secondaryFilter = secondaryFilter, detailed = detailed )
  }
  
  return( scExperiments )
}

# ----- Single Cell  ----

# .download_anndata_with_entry
# - Internal helper to download Single Cell AnnData file with retries. 
.download_anndata_with_retry <- function(url, destfile, max_attempts = 3) {
  
  default_time <- getOption("timeout")
  on.exit(options(timeout = default_time), add = TRUE)
  
  for (attempt in 1:max_attempts) {
    
    # Try download first time for 1 minute. 
    if(attempt == 1){
      message('Downloading....')
    }
    
    # Retry second time for 1 minute. 
    if(attempt == 2){
      message('First attempt failed.')
      message('Trying again for 5 minutes now...')
      options(timeout = if(attempt == 2) 300 else default_time)
    }
    
    # Retry third time for 10 minutes. 
    if(attempt == 3){
      message('Second attempt failed.')
      message('File may be too large')
      message('Will try last time for 15 minutes....')
      options(timeout = if(attempt == 3) 900 else default_time)
    }
    
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


# getAtlasSCExperiment 
# - Download and return the SingleCellExperiemnt object representing a single
#   Single Cell Expression Atlas experiment.
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
      "Downloading Single-Cell Atlas experiment annData into a temp path from:\n",
      fullUrl
    )
  )
  
  
  # Download the H5AD file to a temporary location.
  tempFile <- tempfile(fileext = ".h5ad")
  
  # Ensure the temp file is deleted when the function exits. 
  on.exit( unlink(tempFile) )  
  
  .download_anndata_with_retry(fullUrl, tempFile)
  
  # Read a H5AD file and returns a SingleCellExperiment object, using python version. 
  loadResult <- try(readH5AD( file = tempFile, reader="python", version ="0.7.6", use_hdf5 = FALSE ), silent = TRUE )
  
  # Check whether gene identifiers survived the import.
  utils::head(rownames(loadResult))
  
  # Quit if we got an error.
  if (inherits(loadResult, "try-error")) {
    
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
  
  
  # Check that 'normalised' expression exists in the SCE object. 
  getResult <- try( assays(loadResult)$normalised )
  
  # Stop if retrieving the normalised assay failed or returned NULL.
  if (is.null(getResult) || inherits(getResult, "try-error")) {
    
    stop(
      "ERROR - Download and load appeared successful but no assay 'normalised' matrix was found."
    )
  }
  
  # Else if retrieved. 
  message(
    paste(
      "Successfully loaded assay 'normalised' from SingleCellExperiment object for ",
      experimentAccession
    )
  )
  
  # Check if the experiement has any cell metadata fields. 
  if (ncol(colData(loadResult)) == 0 ) {
    message("No colData found in the SingleCellExperiment object.")
  }
  
  # Return SingleCellExperiment object. 
  SingleCellExperiment::mainExpName(loadResult) <- experimentAccession
  return( loadResult )
  
}

# plotDimRedSCAtlasExperiment 
# - Plots dimensionality reduction (PCA, t-sne, UMAP) of a SingleCellExperiment object. 
plotDimRedSCAtlasExperiment <- function( sceObject, dimRed, colorby, filename = 'sc-dimred', save_pdf = FALSE, plot_layers = NULL) {
  
  # Check if the provided dimRed exists in reducedDimNames.
  if (!(dimRed %in% SingleCellExperiment::reducedDimNames(sceObject))) {
    stop(paste("Error: Dimension reduction method", dimRed, "not found in the object!"))
  }
  
  # Extract dimension reduction coordinates.
  dim_coords <- SingleCellExperiment::reducedDim(sceObject, dimRed)
  
  # Check if colorby exists in colData.
  if (!(colorby %in% colnames(colData(sceObject)))) {
    stop(paste("Error: Color attribute", colorby, "not found in colData!"))
  }
  
  # Convert to a data frame.
  dimred_df <- data.frame(
    axis1 = dim_coords[, 1], 
    axis2 = dim_coords[, 2], 
    colorBy = as.character( colData(sceObject)[, colorby] )
  )
  
  # Plot using ggplot2.
  plotdot <- ggplot(dimred_df, aes(x = .data$axis1, y = .data$axis2, color = .data$colorBy)) +
    geom_point() +
    labs(
      title = paste("Dimensionality Reduction Plot:", dimRed), 
      x = "Component 1", 
      y = "Component 2",
      color = colorby
    ) +
    theme_minimal()
  
  # If not NULL, plot the theme of users choice. 
  if (!is.null(plot_layers)) {
    plotdot <- plotdot + plot_layers
  }
  
  # PDF is saved (only if save_pdf = TRUE).
   if(save_pdf){
     if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
       filename <- paste0(filename, ".pdf")
     }
  
     ggsave(filename, plot = plotdot, width = 8, height = 6)
  }
   return(plotdot)

}

# heatmapSCAtlasEXperiment 
# - Plots heatmap from from a SingleCellExperiment object. 
heatmapSCAtlasExperiment <- function(
    singleCellExperiment,
    genes = NULL,
    sel.K = NULL,
    scaleNormExp = FALSE,
    show_row_names = FALSE,
    geneLabel = c("name", "id"),
    filename = 'sc-heatmap',
    save_pdf = FALSE 
) {
  
  # Get one of the supported versions for geneLabel.
  geneLabel <- match.arg(geneLabel)
  
  # Use markers if genes is NULL. 
  use_markers <- is.null(genes)
  
  # Checking if SingleCellExperiemnt object was supplied. 
  if (!methods::is(singleCellExperiment, "SingleCellExperiment")) {
    stop("The provided object is not a SingleCellExperiment.")
  }
  
  # Checking if normalized assay available. 
  if (!"normalised" %in%
      SummarizedExperiment::assayNames(singleCellExperiment)) {
    stop(
      "The 'normalised' assay is not available ",
      "in the SingleCellExperiment object."
    )
  }
  
  # Extracting normalized assay. 
  normalised_matrix <- SummarizedExperiment::assay(
    singleCellExperiment,
    "normalised",
    withDimnames = FALSE
  )
  
  # Require the normalised matrix to have at least one row and one column.
  if (nrow(normalised_matrix) == 0 ||
      ncol(normalised_matrix) == 0) {
    stop("The normalised matrix has zero rows or columns.")
  }
  
  # Preserve identifiers and ordering from the imported object.
  gene_ids <- rownames(singleCellExperiment)
  cell_ids <- colnames(singleCellExperiment)
  
  if (is.null(gene_ids) ||
      anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids) > 0) {
    stop("The object must have unique, nonempty gene IDs as row names.")
  }
  
  if (is.null(cell_ids) ||
      anyNA(cell_ids) ||
      any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids) > 0) {
    stop("The object must have unique, nonempty cell IDs as column names.")
  }
  
  # Label the matrix rows with gene IDs and the columns with cell IDs.
  rownames(normalised_matrix) <- gene_ids
  colnames(normalised_matrix) <- cell_ids
  
  # Read symbols, keeping their alignment with matrix rows.
  gene_symbols <- rep(NA_character_, length(gene_ids))
  gene_annotations <- SummarizedExperiment::rowData(
    singleCellExperiment
  )
  
  if ("gene_symbols" %in% colnames(gene_annotations)) {
    gene_symbols <- as.character(gene_annotations$gene_symbols)
  }
  
  valid_symbol <- !is.na(gene_symbols) & nzchar(gene_symbols)
  
  if (use_markers) {
    # Retrieve Atlas marker genes for the selected clustering.
    expName <- SingleCellExperiment::mainExpName(
      singleCellExperiment
    )
    
    # This helper must be available, as in the original function.
    if (length(expName) != 1 ||
        is.na(expName) ||
        !.isValidExperimentAccession(expName)) {
      stop(
        "Invalid experiment accession in ",
        "SingleCellExperiment object."
      )
    }
    
    base_url_path <- paste0(
      "https://ftp.ebi.ac.uk/pub/databases/",
      "microarray/data/atlas/sc_experiments"
    )
    
    exp_url_path <- paste(
      base_url_path,
      expName,
      expName,
      sep = "/"
    )
    
    # Read the cluster table and check that it contains K and sel.K.
    clusters_exp <- utils::read.table(
      paste0(exp_url_path, ".clusters.tsv"),
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    
    if (!all(c("K", "sel.K") %in% colnames(clusters_exp))) {
      stop("The cluster file must contain 'K' and 'sel.K' columns.")
    }
    
    if (is.null(sel.K)) {
      sel.K <- unique(
        clusters_exp$K[which(clusters_exp$sel.K == TRUE)]
      )
    }
    
    if (length(sel.K) != 1 ||
        is.na(sel.K) ||
        !sel.K %in% clusters_exp$K) {
      stop("Select exactly one K present in the cluster file.")
    }
    
    # Read the marker gene table for the selected K.
    markers_exp <- utils::read.table(
      paste0(
        exp_url_path,
        ".marker_genes_",
        sel.K,
        ".tsv"
      ),
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE
    )
    
    if (!all(c("genes", "cluster") %in% colnames(markers_exp))) {
      stop("The marker file must contain 'genes' and 'cluster' columns.")
    }
    
    # Match the object's gene IDs to the Atlas marker table.
    # If a gene appears more than once, use its first entry as before.
    marker_index <- match(
      gene_ids,
      as.character(markers_exp$genes)
    )
    
    keep <- !is.na(marker_index)
    
    if (!any(keep)) {
      stop("None of the Atlas marker gene IDs were found in the object.")
    }
    
    gene_clusters <- as.character(
      markers_exp$cluster[marker_index[keep]]
    )
    
    if (anyNA(gene_clusters) || any(!nzchar(gene_clusters))) {
      stop("Some selected marker genes have missing cluster labels.")
    }
  } else {
    # Accept a mixture of gene IDs and symbols.
    genes <- unique(as.character(genes))
    genes <- genes[!is.na(genes) & nzchar(genes)]
    
    keep <- gene_ids %in% genes |
      (valid_symbol & gene_symbols %in% genes)
    
    if (!any(keep)) {
      stop("None of the provided gene IDs or symbols were found.")
    }
    
    matched <- genes %in% gene_ids |
      genes %in% gene_symbols[valid_symbol]
    
    if (any(!matched)) {
      warning(
        "Genes not found: ",
        paste(genes[!matched], collapse = ", ")
      )
    }
    
    gene_clusters <- rep("selected_genes", sum(keep))
  }
  
  selected_ids <- gene_ids[keep]
  selected_symbols <- gene_symbols[keep]
  
  # Subset before dense conversion to reduce memory use.
  dense_matrix <- as.matrix(
    normalised_matrix[keep, , drop = FALSE]
  )
  
  if (isTRUE(scaleNormExp)) {
    dense_matrix <- t(scale(t(dense_matrix)))
    
    # Constant genes have undefined z-scores.
    if (any(!is.finite(dense_matrix))) {
      stop(
        "Scaling produced non-finite values. Check for constant ",
        "genes or missing/non-finite expression values."
      )
    }
  }
  
  # Use separate display labels; retain gene IDs as matrix row names.
  display_labels <- selected_ids
  
  if (geneLabel == "name") {
    has_symbol <- !is.na(selected_symbols) &
      nzchar(selected_symbols)
    
    display_labels[has_symbol] <- selected_symbols[has_symbol]
    
    repeated <- duplicated(display_labels) |
      duplicated(display_labels, fromLast = TRUE)
    
    display_labels[repeated] <- paste0(
      display_labels[repeated],
      " (", selected_ids[repeated], ")"
    )
    
    display_labels <- make.unique(display_labels)
  }
  
  # Assign a distinct colour to each label in gene_clusters.
  cluster_levels <- unique(gene_clusters)
  cluster_colors <- stats::setNames(
    grDevices::rainbow(length(cluster_levels)),
    cluster_levels
  )
  
  # Annotate each marker gene row with its cluster.
  if (use_markers) {
    row_annotation <- ComplexHeatmap::rowAnnotation(
      cluster = gene_clusters,
      col = list(cluster = cluster_colors)
    )
  } else {
    # Annotate the user-selected gene rows and hide the annotation legend.
    row_annotation <- ComplexHeatmap::rowAnnotation(
      selected_genes = gene_clusters,
      col = list(selected_genes = cluster_colors),
      show_legend = FALSE
    )
  }
  
  
  heat_df <- ComplexHeatmap::Heatmap(
    dense_matrix,
    name = if (isTRUE(scaleNormExp)) "Gene z-score" else "Expression",
    col = viridis::viridis(400),
    use_raster = TRUE,
    raster_quality = 5,
    left_annotation = row_annotation,
    row_labels = display_labels,
    show_row_names = show_row_names,
    show_column_names = FALSE,
    cluster_columns = ncol(dense_matrix) > 1,
    cluster_rows = nrow(dense_matrix) > 1,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_split = gene_clusters
  )
  
  # PDF save (only if user choose save_pdf = TRUE)
  if (save_pdf) {
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
    
    pdf(filename, width = 8, height = 8)
    ComplexHeatmap::draw(heat_df)
    dev.off()
  }
  
  return(heat_df)
}

# dorPlotSCAtlasExperiment 
# - Generates a Dot-Plot from a SingleCellExperiment object.
dotPlotSCAtlasExperiment <- function(
    singleCellExperiment,
    genes,
    sel.K = NULL,
    scaleNormExp = FALSE,
    geneLabel = c("name", "id"),
    filename = 'sc-dotplot',
    save_pdf = FALSE,
    dot_theme = NULL 
) {
  
  # Get one of the supported versions for geneLabel.
  geneLabel <- match.arg(geneLabel)
  
  # Checking if SingleCellExperiemnt object was supplied. 
  if (!methods::is(singleCellExperiment, "SingleCellExperiment")) {
    stop("The provided object is not a SingleCellExperiment.")
  }
  
  if (!"normalised" %in%
      SummarizedExperiment::assayNames(singleCellExperiment)) {
    stop(
      "The 'normalised' assay is not available ",
      "in the SingleCellExperiment object."
    )
  }
  
  normalised_matrix <- SummarizedExperiment::assay(
    singleCellExperiment,
    "normalised",
    withDimnames = FALSE
  )
  
  if (nrow(normalised_matrix) == 0 ||
      ncol(normalised_matrix) == 0) {
    stop("The normalised matrix has zero rows or columns.")
  }
  
  # Use the imported object's identifiers and ordering.
  gene_ids <- rownames(singleCellExperiment)
  cell_ids <- colnames(singleCellExperiment)
  
  if (is.null(gene_ids) ||
      anyNA(gene_ids) ||
      any(!nzchar(gene_ids)) ||
      anyDuplicated(gene_ids) > 0) {
    stop("The object must have unique, nonempty gene IDs as row names.")
  }
  
  if (is.null(cell_ids) ||
      anyNA(cell_ids) ||
      any(!nzchar(cell_ids)) ||
      anyDuplicated(cell_ids) > 0) {
    stop("The object must have unique, nonempty cell IDs as column names.")
  }
  
  gene_symbols <- rep(NA_character_, length(gene_ids))
  gene_annotations <- SummarizedExperiment::rowData(
    singleCellExperiment
  )
  
  if ("gene_symbols" %in% colnames(gene_annotations)) {
    gene_symbols <- as.character(gene_annotations$gene_symbols)
  }
  
  rownames(normalised_matrix) <- gene_ids
  colnames(normalised_matrix) <- cell_ids
  
  # Match requests against both gene IDs and gene symbols.
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  
  valid_symbol <- !is.na(gene_symbols) & nzchar(gene_symbols)
  
  keep <- gene_ids %in% genes |
    (valid_symbol & gene_symbols %in% genes)
  
  if (!any(keep)) {
    stop("None of the provided gene IDs or symbols were found.")
  }
  
  matched <- genes %in% gene_ids |
    genes %in% gene_symbols[valid_symbol]
  
  if (any(!matched)) {
    warning(
      "Genes not found: ",
      paste(genes[!matched], collapse = ", ")
    )
  }
  
  selected_ids <- gene_ids[keep]
  selected_symbols <- gene_symbols[keep]
  
  # Construct labels while retaining IDs for grouping.
  display_labels <- selected_ids
  
  if (geneLabel == "name") {
    has_symbol <- !is.na(selected_symbols) &
      nzchar(selected_symbols)
    
    display_labels[has_symbol] <- selected_symbols[has_symbol]
    
    # Distinguish genes that share the same display label.
    repeated <- duplicated(display_labels) |
      duplicated(display_labels, fromLast = TRUE)
    
    display_labels[repeated] <- paste0(
      display_labels[repeated],
      " (", selected_ids[repeated], ")"
    )
    
    display_labels <- make.unique(display_labels)
  }
  
  # Retrieve Atlas cluster assignments.
  expName <- SingleCellExperiment::mainExpName(
    singleCellExperiment
  )
  
  # This helper must be available, as in the original function.
  if (length(expName) != 1 ||
      is.na(expName) ||
      !.isValidExperimentAccession(expName)) {
    stop("Invalid experiment accession in SingleCellExperiment object.")
  }
  
  base_url_path <- paste0(
    "https://ftp.ebi.ac.uk/pub/databases/",
    "microarray/data/atlas/sc_experiments"
  )
  
  exp_url_path <- paste(
    base_url_path,
    expName,
    expName,
    sep = "/"
  )
  
  clusters_exp <- utils::read.table(
    paste0(exp_url_path, ".clusters.tsv"),
    header = TRUE,
    sep = "\t",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  if (!all(c("K", "sel.K") %in% colnames(clusters_exp))) {
    stop("The cluster file must contain 'K' and 'sel.K' columns.")
  }
  
  if (is.null(sel.K)) {
    sel.K <- unique(
      clusters_exp$K[which(clusters_exp$sel.K == TRUE)]
    )
  }
  
  if (length(sel.K) != 1 ||
      is.na(sel.K) ||
      !sel.K %in% clusters_exp$K) {
    stop("Select exactly one K present in the cluster file.")
  }
  
  selected_clustering <- clusters_exp[
    which(clusters_exp$K == sel.K),
    ,
    drop = FALSE
  ]
  
  if (nrow(selected_clustering) != 1) {
    stop("Expected exactly one row for the selected K.")
  }
  
  # Remaining columns contain cell-specific cluster assignments.
  assignment_columns <- !colnames(selected_clustering) %in%
    c("K", "sel.K")
  
  cluster_data <- data.frame(
    SampleID = colnames(selected_clustering)[assignment_columns],
    Cluster = as.character(
      unlist(
        selected_clustering[
          , assignment_columns, drop = FALSE
        ],
        use.names = FALSE
      )
    ),
    stringsAsFactors = FALSE
  )
  
  if (anyDuplicated(cluster_data$SampleID) > 0) {
    stop("The cluster file contains duplicate cell IDs.")
  }
  
  cluster_index <- match(cell_ids, cluster_data$SampleID)
  
  if (anyNA(cluster_index)) {
    stop(
      sum(is.na(cluster_index)),
      " cell IDs in the object were not found in the Atlas ",
      "cluster file. Check that the cell identifiers match."
    )
  }
  
  cell_clusters <- cluster_data$Cluster[cluster_index]
  
  if (anyNA(cell_clusters) || any(!nzchar(cell_clusters))) {
    stop("Some cells have missing cluster assignments.")
  }
  
  # Subset before converting to a dense matrix.
  dense_matrix <- as.matrix(
    normalised_matrix[keep, , drop = FALSE]
  )
  
  if (isTRUE(scaleNormExp)) {
    dense_matrix <- t(scale(t(dense_matrix)))
  }
  
  df <- reshape2::melt(dense_matrix)
  colnames(df) <- c("Gene", "SampleID", "Expression")
  
  df$Gene <- as.character(df$Gene)
  df$SampleID <- as.character(df$SampleID)
  df$Expression <- as.numeric(df$Expression)
  
  df$Cluster <- cell_clusters[
    match(df$SampleID, cell_ids)
  ]
  
  # Treat clusters as categories, using numeric order where possible.
  cluster_levels <- unique(cell_clusters)
  numeric_clusters <- suppressWarnings(
    as.numeric(cluster_levels)
  )
  
  if (all(!is.na(numeric_clusters))) {
    cluster_levels <- cluster_levels[order(numeric_clusters)]
  } else {
    cluster_levels <- sort(cluster_levels)
  }
  
  df$Cluster <- factor(
    df$Cluster,
    levels = cluster_levels
  )
  
  # Aggregate by gene ID before applying display labels.
  df_avg <- dplyr::summarise(
    dplyr::group_by(df, .data$Gene, .data$Cluster),
    Average_Expression = mean(
      .data$Expression,
      na.rm = TRUE
    ),
    .groups = "drop"
  )
  
  df_avg$Gene <- factor(
    df_avg$Gene,
    levels = selected_ids,
    labels = display_labels
  )
  
  expression_label <- if (isTRUE(scaleNormExp)) {
    "Mean gene z-score"
  } else {
    "Mean expression"
  }
  
  dot_plot <- ggplot2::ggplot(
    df_avg,
    ggplot2::aes(
      x = .data$Cluster,
      y = .data$Gene,
      size = .data$Average_Expression,
      color = .data$Average_Expression
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_viridis_c() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 90,
        hjust = 1
      )
    ) +
    ggplot2::labs(
      title = "Dot Plot of Average Gene Expression",
      x = "Cluster",
      y = if (geneLabel == "symbol") "Gene symbol" else "Gene ID",
      size = expression_label,
      color = expression_label
    )
  
  # If not NULL, plot the theme of users choice. 
  if(!is.null(dot_theme)) {
    dot_plot <- dot_plot + dot_theme
  }
  
  # PDF save (only if save_pdf = TRUE).
  if(save_pdf){
    if (!grepl("\\.pdf$", filename, ignore.case = TRUE)) {
      filename <- paste0(filename, ".pdf")
    }
    
    ggsave(filename, plot = dot_plot, width = 8, height = 6)
  }
  return(dot_plot)

  }
