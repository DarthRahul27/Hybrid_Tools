#!/usr/bin/env Rscript

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(optparse)

# Define a function to print usage statement
print_usage <- function() {
  cat("
Usage: ./plot_genomic_data.R -d <datafile>

Options:
  -d, --data     Path to the input data file. The file should be a tab-delimited file
                 containing columns: Chromosome, Running_Position, First_Parent_Count,
                 Second_Parent_Count.

Example:
  ./plot_genomic_data.R -d /path/to/your/datafile.tab
")
}

# Define command-line options
option_list <- list(
  make_option(c("-d", "--data"), type="character", default=NULL,
              help="Path to the input data file", metavar="character")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if the file path is provided
if (is.null(opt$data)) {
  print_usage()
  stop("Input data file is required.", call.=FALSE)
}

# Read the data
data <- read.table(opt$data, header=TRUE, sep="\t")

# Convert Chromosome to factor with levels ordered numerically
data <- data %>%
  mutate(Chromosome = as.factor(Chromosome),
         Chromosome_num = as.numeric(sub("chr", "", Chromosome))) %>%
  arrange(Chromosome_num)

# Calculate cumulative positions for each chromosome
data <- data %>%
  group_by(Chromosome) %>%
  mutate(chr_len = max(Running_Position)) %>%
  ungroup() %>%
  mutate(total_len = cumsum(lag(chr_len, default = 0)),
         cumulative_position = Running_Position + total_len)

# Get chromosome boundaries and midpoints
chromosome_boundaries <- data %>%
  group_by(Chromosome) %>%
  summarize(chr_len = max(Running_Position), 
            start = min(cumulative_position), 
            end = max(cumulative_position)) %>%
  mutate(mid = (start + end) / 2) %>%
  arrange(as.numeric(sub("chr", "", Chromosome))) # Ensure order by numeric values

# Plot using ggplot2
p <- ggplot(data, aes(x = cumulative_position)) +
  geom_point(aes(y = First_Parent_Count, color = "First Parent"), size = 1) +
  geom_point(aes(y = Second_Parent_Count, color = "Second Parent"), size = 1) +
  scale_x_continuous(labels = function(x) x / 1e6, name = "Genomic Position (Mb)") +
  scale_y_continuous(name = "Count") +
  geom_vline(data = chromosome_boundaries, aes(xintercept = end), linetype = "dashed", color = "black") +
  scale_color_manual(values = c("First Parent" = "#A6C6F4", "Second Parent" = "#F4A6A6")) +
  labs(color = "Parent") +
  theme_minimal() +
  theme(
    plot.margin = unit(c(1, 5, 1, 5), "cm"),
    axis.text.x = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("Parent Counts Across Genome") +
  geom_text(data = chromosome_boundaries, aes(x = mid, y = Inf, label = Chromosome), 
            angle = 90, vjust = 1, hjust = 1, size = 3, color = "black")

# Open a PDF device
pdf("genomic_plot.pdf", height = 12, width = 18)
print(p)
dev.off()
