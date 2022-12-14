#'  Find genes positioned in the given region; uses karyoploteR
#'  TxDb.Hsapiens.UCSC.hg19.knownGene, org.Hs.eg.db, and regioneR
#'
#' @param positions - data.frame with 'chromosome_name", 'start', and 'end'
#'   columns (all other columns will be ignored)
#' @param asGRangesList - whether to return the results as GRanges (check
#'    regioneR::GRanges for help) - boolean, defaults to FALSE
#'
#' @export

grabGenes <- function(
  positions,
  asGRangesList = TRUE
){
  requireNamespace('org.Hs.eg.db', quietly = TRUE)
  requireNamespace('AnnotationDbi', quietly = TRUE)
  requireNamespace('Biobase', quietly = TRUE)

  # check the format of 'chromosome_name'
  positions_cp <- positions %>%
    mutate(chromosome_name = ifelse(
      grepl("chr", chromosome_name),
      yes = chromosome_name,
      no = paste0("chr", chromosome_name)
    ))

  genes_data_list <- purrr::map(
  	unique(positions_cp$chromosome_name),
  	function(chrom){
	    cur_data <- positions_cp %>%
	    	filter(chromosome_name == chrom)
	    positions_ranges <- regioneR::toGRanges(
	      data.frame(cur_data)
	    )
	    # where to focus?
	    zoom_range <- range(c(cur_data$start, cur_data$end))
	    zoom_reg <- toGRanges(data.frame(
				chrom,
				zoom_range[1],
				zoom_range[2]
			))

	    pdf(NULL)
	    kp_tmp <- karyoploteR::plotKaryotype(
	      chromosomes = cur_row$chromosome_name,
	      zoom = zoom_reg
	    )
	    genes.data <- suppressMessages(karyoploteR::makeGenesDataFromTxDb(
	      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
	      karyoplot = kp_tmp
	    ))
	    if(!is.null(genes.data)){
	      genes.data <- karyoploteR::addGeneNames(genes.data)
	      out <- karyoploteR::mergeTranscripts(genes.data)
	    } else {
	      out <- NULL
	    }
	    dev.off()
	    return(out)
  })

  if(!asGRangesList){
    genes_df <- purrr::map(genes_data_list, function(d){
      if(length(d$genes) == 0){
        return(NULL)
      }
      return(
        as_tibble(d$genes) %>%
          dplyr::select(chr = seqnames, start, end, strand, gene_name = name)
      )
    }) %>% bind_rows() %>%
      distinct()
      return(genes_df)
  }
  return(genes_data_list)
}
