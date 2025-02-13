\name{ plotDimRedSCAtlasExperiment }
\alias{ plotDimRedSCAtlasExperiment }
\title{ Plot Dimensionality reductions of a Single-Cell Expression Atlas experiment }
\description{
    This function plots dimensionality reduction (PCA, t-sne, UMAP) of 
    a SingleCellExperiment object representing a Single-Cell Expression Atlas experiment.
}
\usage{
    plotDimRedSCAtlasExperiment( sceObject, dimRed, colorby )
}

\arguments{
    \item{ sceObject }{ Single-cell experiment object, representing a Single-Cell Expression Atlas experiment, obtained with the function getAtlasSCExperiment(). }
    \item{ dimRed }{ Dimensionality reduction to plot. Options are: "X_pca", "X_tsne_perplexity_10", "X_umap_neighbors_20", etc. }
    \item{ colorby }{ Column name in colData to color the cells by. }
}

\value{
    A ggplot2 geom point visualisation, which can be further customised by the user.

}

\examples{
    # Example of plotting PCA dimensionality reduction colored by age

    # Download the experiment summary for E-ENAD-19
    enad19 <- getAtlasSCExperiment( "E-ENAD-19" )

    plotDimRedSCAtlasExperiment(enad19, dimRed = "X_pca", colorby = "age" )
}
