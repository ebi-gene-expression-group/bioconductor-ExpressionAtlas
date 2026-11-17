\name{ heatmapAtlasExperiment }
\alias{ heatmapAtlasExperiment }
\title{ Generate HeatMap from normalised expression data for an Expression Atlas experiment }
\description{
    This function generates a heatmap for an Expression Atlas experiment, from normalised expression data.
}
\usage{
    heatmapAtlasExperiment( df, filename = "heatmap", save_pdf = FALSE, show_plot = TRUE, palette = "viridis", markers = TRUE, gene_list = NULL, top_n = 100, scaled=FALSE,
    show_heatmap_title = TRUE, geneLabel = c('name', 'id'))
}
\arguments{
    \item{ df }{
        Dataframe in the output format of getNormalisedAtlasExpression.
    }
    \item{ filename }{
        Filename for the output heatmap file. Default is "heatmap.pdf".
    }
    \item{ save_pdf  }{
        Save PDF file with the heatmap. Default is \code{FALSE}.
    }
    \item{ show_plot }{
        Show plot in the R console. Default is \code{TRUE}.
    }
    \item{ palette }{
        Viridis colour palette for the heatmap. Default is "viridis".
        Other options are: magma, plasma, inferno, cividis, rocket, mako, turbo
    }
    \item{ markers }{
        Top ranked genes per tissue. Default is \code{TRUE}. 
    }
    \item{ gene_list }{
        Gene list of user's choise. Can be ploted with top ranked genes when markers = TRUE.
        If markers = FALSE, will be plotting only genes from gene list. 
    }
    \item{ top_n }{
        Number of top most variable genes to be shown in the heatmap. Default is 100.
    }
    \item{ scaled }{
        Scale the normalised gene expression data. Default is \code{FALSE}.
    }
    \item{ show_heatmap_title }{
        Whether to show a title in the heatmap. Default is \code{TRUE}.
    }
    \item{ geneLabel }{
        Whether to plot with GeneID = 'id' or Gene names = 'name'. Default is 'name'.
    }
    
    }
    

\value{
    A PDF file with the generated heatmap.
}
\examples{
    
    # Download the TPM Expression values for experiment E-MTAB-1733
    mtab1733_tpm <- getNormalisedAtlasExpression( "E-MTAB-1733", "tpm" )

    # Generate heatmap with default parameters
    heatmapAtlasExperiment( mtab1733_tpm, "heatmap", TRUE, FALSE, "viridis", TRUE, NULL, 100, FALSE, TRUE, 'name' )

    # Check if the PDF file was created
    file.exists("heatmap.pdf")

    # Download the CPM Expression values for experiment E-MTAB-1733 with custom parameters
    mtab1733_cpm <- getNormalisedAtlasExpression( "E-MTAB-1733", "cpm" )

    # Generate heatmap with default parameters
    heatmapAtlasExperiment( mtab1733_cpm, "E-MTAB-1733-cpm-heatmap", TRUE, FALSE, "viridis", TRUE, NULL, 100, FALSE, TRUE, 'name' )

    # Check if the PDF file was created
    file.exists("E-MTAB-1733-cpm-heatmap.pdf")

}
