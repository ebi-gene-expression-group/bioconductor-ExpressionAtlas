\name{ heatmapAtlasExperiment }
\alias{ heatmapAtlasExperiment }
\title{ Generate HeatMap from normalised expression data for an Expression Atlas experiment }
\description{
    This function generates a heatmap for an Expression Atlas experiment, from normalised expression data.
}
\usage{
    heatmapAtlasExperiment( df, filename = "heatmap.pdf", heatmap_color = "Blues", top_n = 100, show_heatmap_title = TRUE )
}
\arguments{
    \item{ df }{
        Dataframe in the output format of getNormalisedAtlasExpression.
    }
    \item{ filename }{
        Filename for the output heatmap file. Default is "heatmap.pdf".
    }
    \item{ heatmap_color }{
        Colour palette for the heatmap. Default is "Blues".
    }
    \item{ top_n }{
        Number of top genes to be shown in the heatmap. Default is 100.
    }
    \item{ show_heatmap_title }{
        Whether to show a title in the heatmap. Default is TRUE.
    }
}
\value{
    A PDF file with the generated heatmap.
}
\examples{
    
    # Download the Expression values for experiment E-MTAB-4045
    mtab4045_tpm <- getNormalisedAtlasExpression( "E-MTAB-4045", "tpm" )

    # Generate heatmap with default parameters
    heatmapAtlasExperiment( mtab4045_tpm )

    # Check if the PDF file was created
    file.exists("heatmap.pdf")

    # Prints out the following:
    # TRUE

    # Generate heatmap with custom parameters
    heatmapAtlasExperiment( mtab4045_tpm, "E-MTAB-4045-tpm-heatmap.pdf", "Blues", 100, TRUE )

    # Check if the PDF file was created
    file.exists("E-MTAB-4045-tpm-heatmap.pdf")

    # Prints out the following:
    # TRUE

}
