library("rtracklayer")
library("GenomicRanges")
library("ggplot2")
library(dplyr)


## define regions to estimate background and the exluded proximity region
background_3=GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=103760001, end=110040000))# region corresponding to 3 consecutive TADs located on 3' of the region to be quantified
background_5=GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=93920001, end=96880000))# region corresponding to 3 consecutive TADs located on 5' of the region to be quantified
Target_region=GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=96960001, end=103600000))# reference TAD based on the dataset
exluded_reg=GRanges(seqnames=Rle("chr1"), ranges=IRanges(start=100245000, end=100256300))# region were proximity effects are too high (defined arbitrary based on plot)

Background_tot=c(background_3, background_5)
Target_seq=GenomicRanges::setdiff(Target_region, exluded_reg)

## Imprt bedgraph files into a granges list
# Define directory and files
bedgraph_dir <- choose.dir()
bedgraph_files <- list.files(bedgraph_dir, pattern = "\\.bedgraph$", full.names = TRUE)
bedgraph_files

# Import each BEDGraph file into a GRanges object and store in a list
granges_list <- lapply(bedgraph_files, function(file) {
  import(file, format = "BEDGraph")
})
# Assign names to the list based on the file names
names(granges_list) <- gsub("_1.bedgraph","" ,basename(bedgraph_files))
granges_list

granges_list <- lapply(granges_list, function(gr) {
  resize(gr, width(gr) + 50, fix = "start")
})

### normalize GRanges objects
# Function to normalize GRanges objects
normalize_granges <- function(gr_obj, background_reg=Background_tot, reg_to_norm=Target_seq) {
  # Step 1: Calculate background cutoff and quantify max scores
  backgr_cutoff =subsetByOverlaps(gr_obj, background_reg)$score
  backgr_cutoff<-backgr_cutoff[backgr_cutoff<= quantile(backgr_cutoff, 0.90)]
  backgr_cutoff=max(backgr_cutoff)
  print(paste0("the background value for the dataset ", deparse(substitute(gr_obj)), "is ", backgr_cutoff))

  # Step 2: Subset the GRanges object based on overlaps with target sequence
  gr_norm <- subsetByOverlaps(gr_obj, reg_to_norm)
  
  # Step 3: Normalize scores by subtracting the background value
  gr_norm$score <- gr_norm$score - backgr_cutoff
  
  # Step 4: Replace negative scores with 0
  gr_norm$score[gr_norm$score < 0] <- 0
  
  # Step 5: Calculate total signal after background normalization
  total_signal <- sum(gr_norm$score)
  print(paste0("the total signal value for the dataset ", deparse(substitute(gr_obj)), "is ", total_signal))
  # Step 6: Normalize scores as a percentage of the total signal
  gr_norm$score <- gr_norm$score / total_signal * 100
  gr_norm$cum_score=cumsum(gr_norm$score)
  
  return(gr_norm)
}


# Apply the normalization function to each GRanges object in the list
normalized_granges_list <- lapply(granges_list, normalize_granges,
                                  background_reg=Background_tot, 
                                  reg_to_norm=Target_seq)



### Export the normalized files

output_dir=choose.dir()

# Function to export a GRanges object to a BEDGraph file with a track header
export_granges<- function(gr_obj, file_name, output_dir, 
                                                   track_height = "80", 
                                                   view_limits = "0:0.1")
  {
  # Create the full path for the output file
  output_file <- file.path(output_dir, paste0(file_name, "_norm.bedgraph"))
  
  # Create the track header line with additional details, including maxHeightPixels, viewLimits, and scale
  header <- paste0("track type=bedGraph name=\"", file_name, 
                   "\" visibility=full autoScale=off graphType=bar ",
                   "windowingFunction=maximum smoothingWindow=4 ",
                   "maxHeightPixels=", track_height, " viewLimits=", view_limits
                   )
  # Open the file for writing
  file_conn <- file(output_file, "w")
  
  # Write the track header
  writeLines(header, con = file_conn)
  
  # Export the GRanges object as BEDGraph, appending to the file after the header
  rtracklayer::export(gr_obj, file_conn, format = "BEDGraph")
  
  # Close the file
  close(file_conn)
  
  message("Exported: ", output_file)
}

# Apply the export function to each GRanges object in the list, including a track header
lapply(names(normalized_granges_list), function(file_name) {
  gr_obj <- normalized_granges_list[[file_name]]
  export_granges(gr_obj, file_name, output_dir)
})




### extracting info and Plotting the cumsum values


# Extracting the end coordinates and scores from each GRanges object of the normalized list

# Function to extract relevant data from a GRanges object
extract_data <- function(gr_obj, file_name) {
  data.frame(
    end = end(gr_obj),         # X-axis: end coordinate of the intervals
    score = gr_obj$cum_score,      # Y-axis: normalized score
    file = file_name           # File name to differentiate between files
  )
}


# Combine the extracted data into a single data frame
plot_data <- do.call(rbind, lapply(names(normalized_granges_list), function(file_name) {
  gr_obj <- normalized_granges_list[[file_name]]
  extract_data(gr_obj, file_name)
}))

# Plot the data using ggplot2
condition_colors=c("Ct1_DBT_1_1"="darkolivegreen",
                   "Ct1_DBT_2_1"="#6E8B3D",
                   "Ct2_DBT_1_1"="#F0E68C",
                   "Ct2_DBT_2_1" ="#FFF68F",
                   "Pt1_DBT_1_1"="indianred4",
                   "Pt1_DBT_2_1"="#CD2626",
                   "Pt2_DBT_1_1"="deepskyblue2",
                   "Pt2_DBT_2_1"="#8DEEEE")

plot_all=ggplot(plot_data, aes(x = end, y = score, color = file)) +
  geom_line() +  # Use lines to connect points for each file
  labs(title = "Cumulative Sum of Signal Scores vs End Coordinate", 
       x = "End Coordinate (bp)", 
       y = "Cumulative Sum of Signal Score (%)") +
  theme_minimal() +
  theme(legend.title = element_blank())+
  scale_color_manual(values = condition_colors)  # Apply the specific colors

plot_all

# Define the range of X (end coordinates) you want to plot
x_min <- 100179000  # Start of the range
x_max <- 100969000  # End of the range
# Define the X-axis range for the VP shaded area
shade_xmin <- 100245000  # Start of the shaded area
shade_xmax <- 100256300  # End of the shaded area
# Define the X-axis range for shading the area of interest #1
shade_xmin2 <- 100847000  # Start of the shaded area
shade_xmax2 <- 100948000  # End of the shaded area
# Define the X-axis range for shading the area of interest #1
shade_xmin3 <- 100344000  # Start of the shaded area
shade_xmax3 <- 100356261  # End of the shaded area


# Filter the plot_data to include only the desired range of end coordinates
filtered_plot_data <- plot_data %>%
  filter(end >= x_min & end <= x_max)
# Plot the filtered data
plot_zoom=ggplot(filtered_plot_data, aes(x = end, y = score, color = file)) +
  geom_rect(aes(xmin = shade_xmin, xmax = shade_xmax, ymin = -Inf, ymax = Inf), 
            fill = "ivory2", alpha = 0.2, inherit.aes = FALSE)+  # Add shaded area
  geom_rect(aes(xmin = shade_xmin2, xmax = shade_xmax2, ymin = -Inf, ymax = Inf), 
            fill = "antiquewhite", alpha = 0.8, inherit.aes = FALSE)+  # Add shaded area
  geom_rect(aes(xmin = shade_xmin3, xmax = shade_xmax3, ymin = -Inf, ymax = Inf), 
            fill = "antiquewhite", alpha = 0.8, inherit.aes = FALSE)+ 
  geom_line() +  # Use lines to connect points for each file
  labs(title = "Cumulative Sum of Signal Scores vs End Coordinate", 
       x = "End Coordinate (bp)", 
       y = "Cumulative Sum of Signal Score (%)") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = condition_colors)  # Apply the specific colors

plot_zoom


library("plotly")
ggplotly(plot_zoom)



# Function to extract relevant data from a GRanges object
extract_data2 <- function(gr_obj, file_name) {
  data.frame(
    end = end(gr_obj),         # X-axis: end coordinate of the intervals
    score = gr_obj$score,      # Y-axis: normalized score
    file = file_name           # File name to differentiate between files
  )
}
df<- do.call(rbind, lapply(names(normalized_granges_list), function(file_name) {
  gr_obj <- normalized_granges_list[[file_name]]
  extract_data2(gr_obj, file_name)
}))



#
distal_reg_data1=df[which(df$end> 100027007 & df$end< 100142095),]

df_distal_reg_data1=data.frame(Coord=distal_reg_data1[which(grepl("Ct1_DBT_1", distal_reg_data1$file)),1],
                            Ctrl1_rep1=distal_reg_data1[which(grepl("Ct1_DBT_1", distal_reg_data1$file)), 2], 
                            Ctrl1_rep2=distal_reg_data1[which(grepl("Ct1_DBT_2", distal_reg_data1$file)), 2],
                            Ctrl2_rep1=distal_reg_data1[which(grepl("Ct2_DBT_1", distal_reg_data1$file)), 2],
                            Ctrl2_rep2=distal_reg_data1[which(grepl("Ct2_DBT_2", distal_reg_data1$file)), 2],
                            Pat1_rep1=distal_reg_data1[which(grepl("Pt1_DBT_1", distal_reg_data1$file)), 2],
                            Pat1_rep2=distal_reg_data1[which(grepl("Pt1_DBT_2", distal_reg_data1$file)), 2],
                            Pat2_rep1=distal_reg_data1[which(grepl("Pt2_DBT_1", distal_reg_data1$file)), 2],
                            Pat2_rep2=distal_reg_data1[which(grepl("Pt2_DBT_2", distal_reg_data1$file)),2]
                            )

write.table(df_distal_reg_data1, file = "DBT_distal_reg1.txt", append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
#
distal_reg_data2=df[which(df$end> 100141977 & df$end< 100245119),]

df_distal_reg_data2=data.frame(Coord=distal_reg_data2[which(grepl("Ct1_DBT_1", distal_reg_data2$file)),1],
                               Ctrl1_rep1=distal_reg_data2[which(grepl("Ct1_DBT_1", distal_reg_data2$file)), 2], 
                               Ctrl1_rep2=distal_reg_data2[which(grepl("Ct1_DBT_2", distal_reg_data2$file)), 2],
                               Ctrl2_rep1=distal_reg_data2[which(grepl("Ct2_DBT_1", distal_reg_data2$file)), 2],
                               Ctrl2_rep2=distal_reg_data2[which(grepl("Ct2_DBT_2", distal_reg_data2$file)), 2],
                               Pat1_rep1=distal_reg_data2[which(grepl("Pt1_DBT_1", distal_reg_data2$file)), 2],
                               Pat1_rep2=distal_reg_data2[which(grepl("Pt1_DBT_2", distal_reg_data2$file)), 2],
                               Pat2_rep1=distal_reg_data2[which(grepl("Pt2_DBT_1", distal_reg_data2$file)), 2],
                               Pat2_rep2=distal_reg_data2[which(grepl("Pt2_DBT_2", distal_reg_data2$file)),2]
)

write.table(df_distal_reg_data2, file = "DBT_distal_reg2.txt", append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)
#
distal_reg_data3=df[which(df$end> 100259335 & df$end< 100298553),]

df_distal_reg_data3=data.frame(Coord=distal_reg_data3[which(grepl("Ct1_DBT_1", distal_reg_data3$file)),1],
                               Ctrl1_rep1=distal_reg_data3[which(grepl("Ct1_DBT_1", distal_reg_data3$file)), 2], 
                               Ctrl1_rep2=distal_reg_data3[which(grepl("Ct1_DBT_2", distal_reg_data3$file)), 2],
                               Ctrl2_rep1=distal_reg_data3[which(grepl("Ct2_DBT_1", distal_reg_data3$file)), 2],
                               Ctrl2_rep2=distal_reg_data3[which(grepl("Ct2_DBT_2", distal_reg_data3$file)), 2],
                               Pat1_rep1=distal_reg_data3[which(grepl("Pt1_DBT_1", distal_reg_data3$file)), 2],
                               Pat1_rep2=distal_reg_data3[which(grepl("Pt1_DBT_2", distal_reg_data3$file)), 2],
                               Pat2_rep1=distal_reg_data3[which(grepl("Pt2_DBT_1", distal_reg_data3$file)), 2],
                               Pat2_rep2=distal_reg_data3[which(grepl("Pt2_DBT_2", distal_reg_data3$file)),2]
)

write.table(df_distal_reg_data3, file = "DBT_distal_reg3.txt", append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#
distal_reg_data4=df[which(df$end> 100337197 & df$end< 100372674),]

df_distal_reg_data4=data.frame(Coord=distal_reg_data4[which(grepl("Ct1_DBT_1", distal_reg_data4$file)),1],
                               Ctrl1_rep1=distal_reg_data4[which(grepl("Ct1_DBT_1", distal_reg_data4$file)), 2], 
                               Ctrl1_rep2=distal_reg_data4[which(grepl("Ct1_DBT_2", distal_reg_data4$file)), 2],
                               Ctrl2_rep1=distal_reg_data4[which(grepl("Ct2_DBT_1", distal_reg_data4$file)), 2],
                               Ctrl2_rep2=distal_reg_data4[which(grepl("Ct2_DBT_2", distal_reg_data4$file)), 2],
                               Pat1_rep1=distal_reg_data4[which(grepl("Pt1_DBT_1", distal_reg_data4$file)), 2],
                               Pat1_rep2=distal_reg_data4[which(grepl("Pt1_DBT_2", distal_reg_data4$file)), 2],
                               Pat2_rep1=distal_reg_data4[which(grepl("Pt2_DBT_1", distal_reg_data4$file)), 2],
                               Pat2_rep2=distal_reg_data4[which(grepl("Pt2_DBT_2", distal_reg_data4$file)),2]
)

write.table(df_distal_reg_data4, file = "DBT_distal_reg4.txt", append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#
distal_reg_data5=df[which(df$end> 100843237 & df$end< 100911206),]

df_distal_reg_data5=data.frame(Coord=distal_reg_data5[which(grepl("Ct1_DBT_1", distal_reg_data5$file)),1],
                               Ctrl1_rep1=distal_reg_data5[which(grepl("Ct1_DBT_1", distal_reg_data5$file)), 2], 
                               Ctrl1_rep2=distal_reg_data5[which(grepl("Ct1_DBT_2", distal_reg_data5$file)), 2],
                               Ctrl2_rep1=distal_reg_data5[which(grepl("Ct2_DBT_1", distal_reg_data5$file)), 2],
                               Ctrl2_rep2=distal_reg_data5[which(grepl("Ct2_DBT_2", distal_reg_data5$file)), 2],
                               Pat1_rep1=distal_reg_data5[which(grepl("Pt1_DBT_1", distal_reg_data5$file)), 2],
                               Pat1_rep2=distal_reg_data5[which(grepl("Pt1_DBT_2", distal_reg_data5$file)), 2],
                               Pat2_rep1=distal_reg_data5[which(grepl("Pt2_DBT_1", distal_reg_data5$file)), 2],
                               Pat2_rep2=distal_reg_data5[which(grepl("Pt2_DBT_2", distal_reg_data5$file)),2]
)

write.table(df_distal_reg_data5, file = "DBT_distal_reg5.txt", append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE,
            col.names = TRUE)