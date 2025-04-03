\name{ volcanoDifferentialAtlasExperiment }
\alias{ volcanoDifferentialAtlasExperiment }
\title{ Generate Volcano plots from analytics data for an Expression Atlas experiment }
\description{
    This function generates a volcano plots for an Expression Atlas experiment, from analytics data.
}
\usage{
    volcanoDifferentialAtlasExperiment( df, filename_prefix = "volcano-plot", save_pdf = FALSE, show_plot = TRUE, low_fc_colour = "gray", high_fc_colour = "blue", cutoff = 1, show_volcanoplot_title = TRUE )
}
\arguments{
    \item{ df }{
        Dataframe in the output format of \code{getAnalysticsDifferentialAtlasExpression}.
    }
    \item{ filename_prefix }{
        Filename for the output volcano plot images. Default is "volcano".
    }
    \item{ save_pdf  }{
        Save PDF file with the heatmap. Default is \code{FALSE}.
    }
    \item{ show_plot }{
        Show plot in the R console. Default is \code{TRUE}.
    }
    \item{ low_fc_colour }{
        Colour palette for the low foldchange. Default is "gray".
    }
    \item{ high_fc_colour }{
        Colour palette for the high foldchange. Default is "blue".
    }
    \item{ cutoff }{
        cutoff for the foldchange. Default is 1.
    }
    \item{ show_volcanoplot_title }{
        Whether to show a title in the volcano plot. Default is \code{TRUE}.
    }
}
\value{
    Image file(s) with the generated volcano plots.
}
\examples{
    
    # Download the differential expression analytics data for experiment E-MTAB-10104
    mtab10104_dea <- getAnalysticsDifferentialAtlasExpression( "E-MTAB-10104" )

    # Generate volcano plots with default parameters
    volcanoDifferentialAtlasExperiment( mtab10104_dea )

    # Check if the image files were created
    # file.exists("heatmap.pdf")

    # Generate volcano plots with custom parameters
    volcanoDifferentialAtlasExperiment( mtab10104_dea, "E-MTAB-10104", FALSE, TRUE, "Gray", "Blue", 1, TRUE )

    # Check if the image files were created
    # file.exists("E-MTAB-4045-cpm-heatmap.pdf")

}
