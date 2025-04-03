
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

test_that("Download Expression Atlas (bulk) data for E-MTAB-1624", {

    skip_if_offline()
    expect_identical( names(getAtlasData( "E-MTAB-1624" )), "E-MTAB-1624" )
})

test_that("Search for bulk cancer datasets in human", {

    skip_if_offline()
    cancer_res <- searchAtlasExperiments( query = "cancer", secondaryFilter = "human"  )
    expect_false( ( nrow(cancer_res) == 0 ) ) 
})

test_that("Search for arabidopsis datasets at Single-Cell Expression Atlas", {

    skip_if_offline()
    arabidopsis_res <- searchSCAtlasExperiments( query = "Arabidopsis"  )
    expect_false( ( nrow(arabidopsis_res) == 0 ) ) 
})
