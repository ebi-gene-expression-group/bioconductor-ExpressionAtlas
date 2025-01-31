\name{ heatmapAtlasExperiment }
\alias{ heatmapAtlasExperiment }
\title{ Generate HeatMap from normalised expression data for an Expression Atlas experiment }
\description{
    This function generates heartmap for a single Expression Atlas experiment, from normalised expression data.
}
\usage{
    heatmapAtlasExperiment( df, filename, heatmap_color, top_n, show_heatmap_title )
}
\arguments{
    \item{ df }{
        Dataframe in the output format of getNormalisedAtlasExpression
    }
    \item{ filename }{
        Filename for output heatmap file. Default  
    }
    \item{ heatmap_color }{
        Colour pallate for heatmap. i.e. "Blues"
    }
    \item{ top_n }{
        No of top genes to be shown in Heatmap i.e. 100
    }
    \item{ show_heatmap_title}{
        Show title in Heatmap i.e. TRUE or FALSE
    }
}
\value{
    A PDF file with heatmap.
}
\examples{
    
    # Download the Expression values for experiment E-MTAB-4045
    mtab4045_tpm <- getNormalisedAtlasExpression( "E-MTAB-4045", "tpm" )

    heatmapAtlasExperiment( mtab4045_tpm,  "E-MTAB-4045-tpm-heatmap", "Blues", 100, TRUE )
    
    # See if PDF file created
    file.exists("E-MTAB-4045-tpm-heatmap.pdf")

    # Prints out the following:
    # TRUE

}
