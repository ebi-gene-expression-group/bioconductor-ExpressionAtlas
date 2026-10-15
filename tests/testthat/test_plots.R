# BULK PLOTS ----

# test if returns heatmap for tpm 
test_that("heatmapAtlasExperiment returns a Heatmap for TPM", {
  skip_if_offline()
  
  heatmap_tpm <- getNormalisedAtlasExpression("E-GEOD-74322", "tpm")
  
  result <- heatmapAtlasExperiment( heatmap_tpm, "tpm-heatmap.pdf", FALSE, TRUE,  "viridis", 20, TRUE )
  
  expect_s4_class(result, "Heatmap")
})

# test if returns heatmap for cpm
test_that("heatmapAtlasExperiment returns a Heatmap for CPM", {
  skip_if_offline()
  
  heatmap_cpm <- getNormalisedAtlasExpression("E-GEOD-74322", "cpm")
  
  result <- heatmapAtlasExperiment( heatmap_cpm, "tpm-heatmap.pdf", FALSE, TRUE,  "viridis", 20, TRUE )
  
  expect_s4_class(result, "Heatmap")
})

# test if returns volcano plot

test_that("volcanoDifferentialAtlasExperiment returns a non-empty ggplot", {
  skip_if_offline()
  
  ggplot_test <- getAnalyticsDifferentialAtlasExpression("E-MTAB-10104")
  
  result_ggplot <- volcanoDifferentialAtlasExperiment(
    ggplot_test, "E-MTAB-10104",
    FALSE, FALSE, "Gray", "Blue", 1, TRUE
  )
  
  expect_true("ggplot" %in% class(result_ggplot))
  expect_false(is.null(result_ggplot$data))
  expect_gt(nrow(result_ggplot$data), 0)
})


# SINGLE CELL PLOTS ----

# test if returns scatter plot 
test_that("plotDimRedSCAtlasExperiment returns a ggplot", {
  skip_if_offline()
  
  single_cell_test <- getAtlasSCExperiment("E-GEOD-36552")
  
  reductions <- SingleCellExperiment::reducedDimNames(single_cell_test)
  metadata_columns <- colnames(
    SummarizedExperiment::colData(single_cell_test)
  )
  
  expect_gt(length(reductions), 0)
  expect_gt(length(metadata_columns), 0)
  
  result <- plotDimRedSCAtlasExperiment(
    single_cell_test,
    dimRed = reductions[[1]],
    colorby = metadata_columns[[1]]
  )
  
  expect_s3_class(result, "ggplot")
  expect_equal(nrow(result$data), ncol(single_cell_test))
})

# test if returns heatmap

test_that("heatmapSCAtlasExperiment returns heatmap for single cell experiment", {
  skip_if_offline()
  
  heatmap_single_test <- getAtlasSCExperiment( "E-GEOD-36552" )
  
  heatmap_single_result <- heatmapSCAtlasExperiment( heatmap_single_test, genes=NULL, sel.K=NULL, 
                                                     scaleNormExp=TRUE, show_row_names=FALSE ) 
  expect_s4_class(heatmap_single_result, "Heatmap")
})

# test if returns dotPlotSCAtlasExperiment

test_that("dotPlotSCAtlasExperiment returns dot plot for single cell experiment", {
  skip_if_offline()
  
  dotplot_exp <- getAtlasSCExperiment( "E-GEOD-36552" )
  
  result_dotplot <- dotPlotSCAtlasExperiment(dotplot_exp, genes=c('ENSG00000166681','ENSG00000178928', 
                                                                'ENSG00000142182' , 'ENSG00000160282' ),
                                            sel.K=4, scaleNormExp=TRUE) + theme_classic()
  expect_true("ggplot" %in% class(result_dotplot))
  
})

# test if dotPlotSCAtlasExperiment contains data: 
test_that("dotPlotSCAtlasExperiment returns data", {
  skip_if_offline()
  
  dotplot_exp <- getAtlasSCExperiment( "E-GEOD-36552" )
  
  result_dotplot <- dotPlotSCAtlasExperiment(dotplot_exp, genes=c('ENSG00000166681','ENSG00000178928', 
                                                                  'ENSG00000142182' , 'ENSG00000160282' ),
                                             sel.K=4, scaleNormExp=TRUE) + theme_classic()
  
  expect_false(is.null(result_dotplot$data))
  expect_gt(nrow(result_dotplot$data), 0)
})






