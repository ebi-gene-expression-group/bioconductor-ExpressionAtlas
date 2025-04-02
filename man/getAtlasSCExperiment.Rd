\name{ getAtlasSCExperiment }
\alias{ getAtlasSCExperiment }
\title{ Download data for a Single-Cell Expression Atlas experiment }
\description{
    This function downloads and returns a SingleCellExperiment object representing a
    Single-Cell Expression Atlas experiment.
}
\usage{
    getAtlasSCExperiment( experimentAccession )
}
\arguments{
    \item{ experimentAccession }{
        Single-cell experiment accession e.g. "E-ENAD-19"
    }
}
\value{
    A SingleCellExperiment object representing a single-Cell Expression Atlas 
    experiment.
    
    Each SCE object contains the following:

    \itemize{
        \item{ 
            Assays, containing the raw counts (X), filtered counts (filtered), and normalised
            and filtered counts (normalised).
            }
        \item{ Column data (cell-level metadata) }
        \item{ Different dimensionality reductions (e.g. PCA, t-SNE, UMAP) }
    }

}
\examples{
    
    # Download the experiment summary for E-GEOD-11175
    enad19 <- getAtlasSCExperiment( "E-ENAD-19" )

    enad19
    # Prints out the following:
    #
    # class: SingleCellExperiment 
    # dim: 3644 26 
    # metadata(0):
    # assays(3): X filtered normalised
    # rownames: NULL
    # rowData names(1): PCs
    # colnames(26): ERR2192535 ERR2192536 ... ERR2192559 ERR2192560
    # colData names(27): age cell_type ... louvain_resolution_1.0
    #   louvain_resolution_2.0
    # reducedDimNames(21): X_pca X_tsne_perplexity_1 ...
    #   X_umap_neighbors_n_neighbors_5 X_umap_neighbors_n_neighbors_50
    # mainExpName: NULL
    # altExpNames(0):

}
