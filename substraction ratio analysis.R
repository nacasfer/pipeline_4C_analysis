# Function to compute the subtraction of scores for overlapping regions
subtract_tracks <- function(gr1, gr2) {
  # Find overlapping regions
  overlaps <- findOverlaps(gr1, gr2)
  
  # Extract the overlapping ranges
  gr1_overlap <- gr1[queryHits(overlaps)]
  gr2_overlap <- gr2[subjectHits(overlaps)]
  
  # Ensure the intervals are identical before computing subtraction
  start(gr1_overlap) <- pmax(start(gr1_overlap), start(gr2_overlap))
  end(gr1_overlap) <- pmin(end(gr1_overlap), end(gr2_overlap))
  start(gr2_overlap) <- start(gr1_overlap)
  end(gr2_overlap) <- end(gr1_overlap)
  
  # Compute the subtraction of the scores
  Subtract <- mcols(gr1_overlap)$score - mcols(gr2_overlap)$score 
  
  # Create a new GRanges object with the overlapping intervals and the subtracted score
  result_gr <- GRanges(seqnames = seqnames(gr1_overlap),
                       ranges = IRanges(start = start(gr1_overlap), end = end(gr1_overlap)),
                       score = Subtract)
  
  return(result_gr)
}

# Function to subtract scores between two BEDGraph files and save output
subtract_tracks_bedgraph <- function(file1, file2, output_file) {
  # Load the two BEDGraph files
  gr1 <- file1
  gr2 <- file2
  
  # Compute subtraction for overlapping regions
  gr_subtracted <- subtract_tracks(gr1, gr2)
  
  # Export the result to a new BEDGraph file
  rtracklayer::export(gr_subtracted, con = output_file, format = "BEDGraph")
  
  message("Subtracted BEDGraph file saved to: ", output_file)
}

### Mean of normalized profiles
bedgraph_dir <- choose.dir()
bedgraph_files <- list.files(bedgraph_dir, pattern = "\\.bedgraph$", full.names = TRUE)
bedgraph_files

# Import each BEDGraph file into a GRanges object and store in a list
Mean_bedgraph <- lapply(bedgraph_files, function(file) {
  import(file, format = "BEDGraph")
})
# Assign names to the list based on the file names
names(Mean_bedgraph) <- gsub(".bedgraph","" ,basename(bedgraph_files))
Mean_bedgraph


### Call the function  for each comparison

subtract_tracks_bedgraph(file1 = Mean_bedgraph$Mean_Pat1_DBT_corm, 
                         file2 = Mean_bedgraph$Mean_Ctrl1_DBT_corm, 
                         output_file = "Subtract_Pat1_Ctrl1.bedgraph")

subtract_tracks_bedgraph(file1 = Mean_bedgraph$Mean_Pat1_DBT_corm, 
                         file2 = Mean_bedgraph$Mean_Ctrl2_DBT_corm, 
                         output_file = "Subtract_Pat1_Ctrl2.bedgraph")

subtract_tracks_bedgraph(file1 = Mean_bedgraph$Mean_Pat1_DBT_corm, 
                         file2 = Mean_bedgraph$Mean_Pat2_DBT_corm, 
                         output_file = "Subtract_Pat1_Pat2.bedgraph")

