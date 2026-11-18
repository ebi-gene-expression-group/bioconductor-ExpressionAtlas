# NORMAL USE BULK  ----

# test for seacrhAtlasExperiment to return Data Frame:
#needs internet connection: 
test_that("searchAtlasExperiment returns data frame", {
  skip_if_offline()
  
  search_data <- searchAtlasExperiments( query = "lung" , detailed = FALSE )
  
  expect_s4_class(search_data, "DataFrame")
})

# test for getAtlasExperiement to return expected output: 
test_that("getAtlasExperiment returns a non-empty SimpleList", {
  skip_if_offline()
  
  simple_list <- getAtlasExperiment("E-MTAB-3007")
  
  expect_true("SimpleList" %in% class(simple_list))

  
  expect_true(nrow(simple_list$rnaseq) > 0)  
  expect_true(ncol(simple_list$rnaseq) > 0)
  
})


# test for getAtlasData to return expected output: 
test_that("getAtlasData returns a non-empty SimpleList", {
  skip_if_offline()
  
  simple_list_data <- getAtlasData(c("E-GEOD-11175", "E-MTAB-1625", "E-MTAB-1624"))
  
  expect_s4_class(simple_list_data, "SimpleList")
  expect_gt(length(simple_list_data), 0)
})


# test for getNormalizedAtlasExperiment to return data frame: 
test_that("getNormalisedAtlasExpression return data frame", {
  skip_if_offline()
  
  data_frame_tpm <- getNormalisedAtlasExpression("E-MTAB-1733", "tpm")
  data_frame_cpm <- getNormalisedAtlasExpression("E-MTAB-1733", "cpm")
  
  expect_s3_class(data_frame_tpm, 'data.frame')
  expect_true(nrow(data_frame_tpm) > 0)
  expect_true(ncol(data_frame_tpm) > 0)
  
  expect_s3_class(data_frame_cpm, 'data.frame')
  expect_true(nrow(data_frame_cpm) > 0)
  expect_true(ncol(data_frame_cpm) > 0)
})


# test for getAnalyticsDifferentialAtlasExpression() to return data frame: 
test_that("getAnalyticsDifferentialAtlasExpression returns data frame for differential bulk data", {
  skip_if_offline()
  
  data_analytics <-  getAnalyticsDifferentialAtlasExpression( "E-MTAB-10104" )
  
  expect_s3_class(data_analytics, 'data.frame')
  expect_true(nrow(data_analytics) > 0)
  expect_true(ncol(data_analytics) > 0)
  
})


# NORMAL USE SINGLE CELL ----

# test for  searchSCAtlasExperiments() to return Data Frame:
test_that("searchSCAtlasExperiments returns data frame", {
  skip_if_offline()
  
  search_single_cell <- searchSCAtlasExperiments(query = "endoderm")
  expect_s4_class(search_single_cell, "DataFrame")
  
})

# test for getAtlasSCExperiment() to return Data Frame:
test_that("getAtlasSCExperiment() returns SingleCellExperiment object", {
  skip_if_offline()
  
  single_cell_data <- getAtlasSCExperiment("E-ENAD-19")
  expect_s4_class(single_cell_data,'SingleCellExperiment')
})

# test for getAtlasSCExperiment() to return not NULL: 
# nrow for genes 
# ncol for cells 
test_that("getAtlasSCExperiment() SingleCellExperiment does not return NULL", {
  skip_if_offline()
  
  single_cell_exp <- getAtlasSCExperiment("E-ENAD-19")
  
  expect_false(is.null(single_cell_exp))
  expect_true(nrow(single_cell_exp) > 0) 
  expect_true(ncol(single_cell_exp) > 0) 
  
})


# ERROR HANDLING BULK ----

#test how searchAtclasexperiment handles empty query
test_that("searchAtlasExperiments rejects a missing query", {
  expect_error(
    searchAtlasExperiments(),
    "Please provide at least one query search term."
  )
})

#test searchAtlasExperiment rejects a non-charater query 
test_that("searchAtlasExperiments rejects a non-character query", {
  expect_error(
    searchAtlasExperiments(query = 123),
    "Please provide query search term(s) as a character vector.",
    fixed=TRUE
  )
})

#testing getAtlasexperiment accession input 
test_that("Checking valid asseccion for getAtlasExperiment function", {
  expect_warning( 
    expect_error(getAtlasExperiment("EMTAB9090"),
                 "Experiment accession not valid"),
    "does not look like an ArrayExpress"
  )
})


#testing getAtlasexperiment accession input, nothing 
test_that("Checking empty input for getAtlasExperiment function", {
  expect_warning( 
    expect_error(getAtlasExperiment(),
                 "Experiment accession not valid. Cannot continue"),
    "Accession missing. Cannot validate"
  )
})

#testing getAtlasData data missing input 
test_that("Checking valid vector input for getAtlasData", {
  expect_error(getAtlasData(),
               "Please provide a vector of experiment accessions to download.")
})

#testing getAtlasData will reject invalid accessions 
test_that("Checking getAtladData fucntion rejects fake asseccions", {
  fake_accessions <- c('E1', 'E2')
  expect_warning(
    expect_error(getAtlasData(fake_accessions), 
                 "None of the accessions passed are valid ArrayExpress/BioStudies accessions. Cannot continue."),
    "does not look like an ArrayExpress")
})

#test getNormalisedAtlasExpression will reject empty input 
test_that("getNormalisedAtlasExpression return data frame", {
  expect_warning(
    expect_error(
    getNormalisedAtlasExpression(),
    "Experiment accession not valid. Cannot continue."
    ),
    "Accession missing. Cannot validate.")
 
})

# test getAnalyticsDifferentialAtlasExpression will reject empty query
test_that("getAnalyticsDifferentialAtlasExpression rejects empty query", {
  expect_warning(
    expect_error(
    getAnalyticsDifferentialAtlasExpression(),
    "Experiment accession not valid. Cannot continue."),
    "Accession missing. Cannot validate.")
})



# ERROR HANDLING SINGLE CELL ----

# test searchSCAtlasExperiments will reject empty string 
test_that("searchSCAtlasExperiments returns data frame", {
  expect_error(
    searchSCAtlasExperiments(),
    "Please provide at least one query search term."
  )
})

# test searchSCAtlasExperiments will reject non-character query 
test_that("searchSCAtlasExperiments rejects a non-character query", {
  expect_error(
    searchSCAtlasExperiments(query = 123),
    "Please provide query search term(s) as a character vector.",
    fixed=TRUE
  )
})

# test searchSCAtlasExperiments will reject non-character secondaryFilter 
test_that("searchSCAtlasExperiments rejects a secondFilter as non-character vector", {
  expect_error(
    searchSCAtlasExperiments(query = 'lung', secondaryFilter =  123),
    "Please provide secondaryFilter as a character vector."
  )
})

# test getAtlasSCExperiment will reject empty input
test_that("Checking valid asseccion for getAtlasSCExperiment function", {
  expect_warning( 
    expect_error(getAtlasSCExperiment(),
                 "Experiment accession not valid"),
    "Accession missing. Cannot validate."
  )
})

# test getAtlasSCExperiment will reject wrong accession number 
test_that("Checking valid asseccion for getAtlasSCExperiment function", {
  expect_warning( 
    expect_error(getAtlasSCExperiment("E-ENAD_000"),
                 "Experiment accession not valid"),
    "does not look like an ArrayExpress"
  )
})



