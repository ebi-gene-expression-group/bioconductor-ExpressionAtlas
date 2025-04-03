\name{ getAtlasExperiment }
\alias{ getAtlasExperiment }
\title{ Download data for an Expression Atlas experiment }
\description{
    This function downloads and returns a SimpleList object representing a
    single Expression Atlas experiment, based on the ArrayExpress accession of
    the experiment.
}
\usage{
    getAtlasExperiment( experimentAccession )
}
\arguments{
    \item{ experimentAccession }{
        ArrayExpress experiment accession e.g. "E-GEOD-11175"
    }
}
\value{
    A SimpleList object representing a single Expression Atlas experiment. The
    SimpleList contains one entry per platform used in the experiment. For
    sequencing experiments, there is a single entry in the list. For microarray
    experiments, there is one entry per array design used. Currently Expression
    Atlas does not support multi-technology (e.g. microarray and RNA-seq)
    experiments.

    For a single-channel microarray experiment, each entry of the list is an
    ExpressionSet object. For a sequencing experiment, the single entry is a
    SummarizedExperiment object.  Please refer to the relevant documentation on
    these classes for more information about them.
    
    \bold{ RNA-seq data }

    Each SummarizedExperiment object contains the following:

    \itemize{
        \item{ 
            Matrix of raw counts (not normalized), in the assays slot, in a
            counts element. 
            }
        \item{ Sample annotations, in the colData slot. }
        \item{ 
            Brief outline of methods, from QC of FASTQ files to production of
            raw counts, in the exptData slot. 
            }
    }

    \bold{ Single-channel microarray data }

    Each ExpressionSet object contains the following:

    \itemize{
        \item{ 
            Matrix of normalized intensity values, in the assayData, accessed
            via: exprs( expressionSet ) 
            }
        \item{ 
            Sample annotations, in the phenoData, accessed via: pData(
            expressionSet ) 
            }
        \item{ 
            Brief outline of normalization method applied, in the
            experimentData slot, accessed via: preproc( experimentData(
            expressionSet ) )
            }
    }
}
\examples{
    
    # Download the experiment summary for E-GEOD-11175
    geod11175 <- getAtlasExperiment( "E-GEOD-11175" )

    # See the entries available (in this case array design accessions)
    names( geod11175 )
    # Prints out the following:
    # [1] "A-AFFY-126"

    # Get the only ExpressionSet object from this experiment.
    eset <- geod11175[[ "A-AFFY-126" ]]
}
