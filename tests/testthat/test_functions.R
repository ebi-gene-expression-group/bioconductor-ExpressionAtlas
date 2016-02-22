
context( "Testing functions from ExpressionAtlas package" )

test_that( "Accession validation returns true or false at the right times", {
    
    expect_true( .isValidExperimentAccession( "E-MTAB-3007" ) )
    expect_false( .isValidExperimentAccession( "DRP000391" ) )
    expect_false( .isValidExperimentAccession( ) )

} )
