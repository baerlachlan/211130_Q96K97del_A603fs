# Use with objects from de.Rmd loaded

lapply(1:25, function(chr){
  genes <- topTables_chosen$`MPS-IIIB` %>%
    dplyr::filter(DE, chromosome == chr) %>%
    pull(gene_id)
  cat(paste0(chr, "\n"))
  cat(genes, sep="\n")
})

topTables_chosen$`MPS-IIIB` %>%
  dplyr::filter(DE) %>%
  pull(chromosome) %>%
  table()
