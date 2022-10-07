#' grab gene regions from ensembl - using biomaRt
#'
#' @param positions - positions of markers for which to search the
#'    regions in the mart, given as a data.frame with columns named as
#'    filters that are given below
#' @param mart - which database to use - check description in biomaRt::getBM
#' @param attribs - which attributes to fetch - check description in biomaRt::getBM
#' @param filters - which filters to use - check description in biomaRt::getBM
#' @param asGRanges - whether to return the results as GRanges (check
#'    regioneR::GRanges for help) - boolean, defaults to FALSE
#' @param genom_ver - which genome version? default: NULL - newest version
#' 

grabGenesEnsembl <- function(
  positions,
  attribs,
  filters,
  asGRanges = FALSE,
  genome_ver = NULL
){
  if(!requireNamespace('biomaRt')){
    stop("Couldn't find 'biomaRt' package - please install and try once more.",
         call. = FALSE)
  }
  requireNamespace('biomaRt')
  
  if(!all(names(positions) %in% filters)){
    stop("Columns of 'positions' data.frame should be the same as filters!",
         call. = FALSE)
  }

  if(is.null(genome_ver)){
    ensembl_gene <- useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl"
    )
  } else {
    ensembl_gene <- useEnsembl(
      biomart = "ensembl",
      dataset = "hsapiens_gene_ensembl",
      GRCh = genome_ver
    )
  }

  cur_genes_ensembl <- grabRegulRegions(
    positions,
    mart = ensembl_gene,
    attribs = attribs,
    filters = filters
  )
  cur_genes_ensembl <- as_tibble(cur_genes_ensembl) %>%
    filter(uniprotswissprot != "")
  
  if(asGRanges){
    if(nrow(cur_genes_ensembl) == 0){
      return(NULL)
    }
    cur_genes_ensembl <- regioneR::toGRanges(
      as.data.frame(
        cur_genes_ensembl %>%
          mutate(seqnames = paste0("chr", chromosome_name)) %>%
          dplyr::select(
            seqnames,
            start = start_position,
            end = end_position,
            everything()
          ) %>%
          distinct()
      )
    )
  }

  return(cur_genes_ensembl)
}
