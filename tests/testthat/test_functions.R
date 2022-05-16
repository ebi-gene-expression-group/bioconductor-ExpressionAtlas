
test_that( "Accession validation returns true or false at the right times", {
    
    expect_true( .isValidExperimentAccession( "E-MTAB-3007" ) )
    expect_false( .isValidExperimentAccession( "DRP000391" ) )
    expect_false( .isValidExperimentAccession( ) )
})

test_that("check test data", {

    data( "atlasRes" )
    data( "rnaseqExps" )
    expect_equal( nrow(atlasRes), 3 )
    expect_equal( ncol(atlasRes), 4 )
    expect_identical( names(rnaseqExps), "E-MTAB-1625" )
})

# the following tests require internet connection

check_api <- function() {
    if ( is.character(getURL("www.ebi.ac.uk/arrayexpress/")) == FALSE ) {
        skip("API not available")
    }
}

test_that("Download the experiment summary for E-GEOD-11175", {

    check_api()
    geod11175 <- getAtlasExperiment( "E-GEOD-11175" )
    expect_identical( names( geod11175 ), "A-AFFY-126" )
})

test_that("Search for cancer datasets in human", {

    check_api()
    cancerRes <- searchAtlasExperiments( properties = "cancer", species = "human"  )
    expect_false( ( nrow(cancerRes) == 0 ) ) 
})
