#valid experiment accession check: 
test_that( "Valid accessions are accepted", {
    
    expect_true( .isValidExperimentAccession( "E-MTAB-3007" ) )
    expect_true( .isValidExperimentAccession( "E-GEOD-36552" ) )

})

test_that("Accession validation rejects an invalid accession", {
  expect_warning(
    result <- .isValidExperimentAccession("EMTAB3007"),
    "does not look like an ArrayExpress"
  )
  
  expect_false(result)
})

test_that("Accession validation rejects empty accession",{
  expect_warning(
    result <- .isValidExperimentAccession(),
    "Accession missing"
  )
  expect_false(result)
} )


#checking whether data is valid from atlaRes from vignette: 
test_that("check test data", {

    data( "atlasRes" )
    data( "rnaseqExps" )
    expect_equal( nrow(atlasRes), 3 )
    expect_equal( ncol(atlasRes), 4 )
    expect_identical( names(rnaseqExps), "E-MTAB-1625" )
})

# the following tests require internet connection

test_that("Download data for E-MTAB-1624", {

    skip_if_offline()
    expect_identical( names(getAtlasData( "E-MTAB-1624" )), "E-MTAB-1624" )
})

test_that("Search for cancer datasets in human", {

    skip_if_offline()
    cancer_res <- searchAtlasExperiments( properties = "cancer", species = "human"  )
    expect_false( ( nrow(cancer_res) == 0 ) ) 
})
