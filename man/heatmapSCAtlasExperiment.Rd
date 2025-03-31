\name{ heatmapSCAtlasExperiment }
\alias{ heatmapSCAtlasExperiment }
\title{ Visualise heatmap for Single-Cell Expression Atlas experiment }
\description{
    This function plots a heatmap from a SingleCellExperiment object coming
    from a Single-Cell Expression Atlas experiment.
}
\usage{
    heatmapSCAtlasExperiment( singleCellExperiment, genes=NULL, sel.K=NULL, scaleNormExp=FALSE, show_row_names=FALSE )
}
\arguments{
    \item{singleCellExperiment}{
        A SingleCellExperiment object from Single Cell Expression Atlas
    }
    \item{genes}{
        Vector of gene IDs to include in the heatmap. If NULL, all genes will be shown for the default clustering. Default is NULL.
    }
    \item{sel.K}{
        Number of clusters for Single Cell Expression Atlas Experiment. If NULL, a gene list must be provided. If both \code{sel.K} and \code{genes} are \code{NULL}, the function will use default clustering for all markers genes for this experiment.
    }
    \item{scaleNormExp}{
        Logical indicating whether to scale normalized expression values, which adjusts the data to have a mean of zero and a standard deviation of one. Default is FALSE.
    }
    \item{show_row_names}{
        Logical indicating whether to show gene names in rows. Default is FALSE.
    }
}
\value{
    A Heatmap-class object from the ComplexHeatmap package. This object can be further customized or combined with other heatmaps using the ComplexHeatmap package's functions.
}
\examples{
    
    egeod6552 <- getAtlasSCExperiment( "E-GEOD-36552" )

    heatmap6552 <- heatmapSCAtlasExperiment( egeod6552, genes=NULL, sel.K=NULL, scaleNormExp=FALSE, show_row_names=FALSE )

    # Example with a specific gene list
    heatmap6552_genes <- heatmapSCAtlasExperiment( egeod6552, genes=c('ENSG00000151611','ENSG00000020577', 'ENSG00000188869' ), sel.K=NULL, scaleNormExp=FALSE, show_row_names=TRUE )
}
