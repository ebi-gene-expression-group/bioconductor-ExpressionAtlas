\name{ heatmapAtlasExperiment }
\alias{ heatmapAtlasExperiment }
\title{ Generate HeatMap from normalised expression data for an Expression Atlas experiment }
\description{
    This function generates a heatmap for an Expression Atlas experiment, from normalised expression data.
}
\usage{
    heatmapAtlasExperiment( df, filename = "heatmap", save_pdf = FALSE, show_plot = TRUE, heatmap_color = "Blues", top_n = 100, show_heatmap_title = TRUE )
}
\arguments{
    \item{ df }{
        Dataframe in the output format of getNormalisedAtlasExpression.
    }
    \item{ filename }{
        Filename for the output heatmap file. Default is "heatmap.pdf".
    }
    \item{ save_pdf  }{
        Save PDF file with the heatmap. Default is FALSE.
    }
    \item{ show_plot }{
        Show plot in the R console. Default is TRUE.
    }
    \item{ heatmap_color }{
        Colour palette for the heatmap. Default is "Blues".
    }
    \item{ top_n }{
        Number of top most variable genes to be shown in the heatmap. Default is 100.
    }
    \item{ show_heatmap_title }{
        Whether to show a title in the heatmap. Default is TRUE.
    }
}
\value{
    A PDF file with the generated heatmap.
}
\examples{
    
    # Download the TPM Expression values for experiment E-MTAB-4045
    mtab4045_tpm <- getNormalisedAtlasExpression( "E-MTAB-4045", "tpm" )

    # Generate heatmap with default parameters
    heatmapAtlasExperiment( mtab4045_tpm, "heatmap", TRUE, FALSE, "Blues", 100, TRUE )

    # Check if the PDF file was created
    file.exists("heatmap.pdf")


    # Download the CPM Expression values for experiment E-MTAB-4045 with custom parameters
    mtab4045_cpm <- getNormalisedAtlasExpression( "E-MTAB-4045", "cpm" )

    # Generate heatmap with default parameters
    heatmapAtlasExperiment( mtab4045_cpm, "E-MTAB-4045-cpm-heatmap", TRUE, FALSE, "Blues", 100, TRUE )

    # Check if the PDF file was created
    file.exists("E-MTAB-4045-cpm-heatmap.pdf")


    # Prints out the following:
    # TRUE

}
