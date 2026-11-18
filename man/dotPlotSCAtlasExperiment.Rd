\name{ dotPlotSCAtlasExperiment }
\alias{ dotPlotSCAtlasExperiment }
\title{ Dot-plot for a Single-Cell Expression Atlas experiment }
\description{
    This function generates a Dot-Plot from a SingleCellExperiment object coming
    from a Single-Cell Expression Atlas experiment.
}
\usage{
    dotPlotSCAtlasExperiment( singleCellExperiment, genes, sel.K=NULL, scaleNormExp=FALSE, geneLabel = c('name', 'id'), filename = 'sc-dotplot', save_pdf = FALSE, dot_theme = NULL )
}
\arguments{
    \item{singleCellExperiment}{
        A SingleCellExperiment object from Single Cell Expression Atlas
    }
    \item{genes}{
        Vector of gene IDs to include in the dot-plot.
    }
    \item{sel.K}{
        Number of clusters for Single Cell Expression Atlas Experiment. If both \code{sel.K} is \code{NULL}, the function will use default clustering for this experiment.
    }
    \item{scaleNormExp}{
        Logical indicating whether to scale normalized expression values, which adjusts the data to have a mean of zero and a standard deviation of one. Default is \code{FALSE}.
    }
    \item{geneLabel}{
      Whether to plot with Gene IDs = 'id' or Gene names = 'name'. Default is 'name'.
    }
    \item{ filename }{
         File name or path for the heatmap PDF. Defaults to
        \code{"sc-dotplot"}, which is saved as \code{sc-dotplot.pdf}.
    }
    \item{ save_pdf  }{
        Save PDF file with the heatmap. Default is \code{FALSE}.
    }
    \item{ dot_theme }{
        Allows user to put a list of ggplot2 components, such as themes, to be added before saving and returning the dot plot. 
        e.g. dot_theme = list(ggplot2::theme_classic())
    }
}
\value{
    A Heatmap-class object from the ComplexHeatmap package. This object can be further customized or combined with other heatmaps using the ComplexHeatmap package's functions.
}
\examples{

    egeod6552 <- getAtlasSCExperiment( "E-GEOD-36552" )

    dotPlotSCAtlasExperiment(egeod6552, genes=c('ENSG00000166681','ENSG00000178928', 'ENSG00000142182' , 'ENSG00000160282' ), sel.K=4, save_pdf = FALSE, dot_theme = NULL)

}
