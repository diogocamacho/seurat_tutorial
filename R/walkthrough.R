#' libraries et al
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(here)
library(EnsDb.Hsapiens.v75)
library(AnnotationDbi)
library(BSgenome.Hsapiens.UCSC.hg19)

#' load data
count_data <- Seurat::Read10X(data.dir = here::here("./data/pbmck3/filtered_gene_bc_matrices/hg19/"), gene.column = 1)


#' Using the `Seurat` package to load the 10X genomics data here, while keeping the Ensembl IDs. Let's now update the gene names using the `EnsDb.Hsapiens.v79` package.

#' gene symbols
gene_symbols <- AnnotationDbi::mapIds(x = EnsDb.Hsapiens.v75, keys = count_data@Dimnames[[1]], column = c("SYMBOL", "GENENAME"), keytype = "GENEID")


#' gene filtering
# this will tell you in how many cells a given gene is found
presence_cells <- apply(count_data, 1, function(y) length(which(y != 0)))

# this will tell you how many genes are found in a given cell
genes_per_cell <- apply(count_data, 2, function(y) length(which(y != 0)))

#' filter data
#' we will filter the data based on the number of genes found in a given cell and the number of cells a given gene is found in. we will also account for the percentage of mitochondrial genes, selecting only those cells that have a small percentage of mitochondrial genes.
gids <- which(presence_cells >= 3)
sids <- which(genes_per_cell > 200 & genes_per_cell < 2500)

new_counts <- count_data[gids, sids]
new_symbols <- gene_symbols[gids]

#' mitochondrial genes
mito_genes <- grep(pattern = "^MT-", x = new_symbols)


num_reads_sample <- apply(new_counts, 2, sum)
num_mit_reads <- apply(new_counts[mito_genes, ], 2, sum)
percent_mit <- apply(cbind(num_mit_reads, num_reads_sample), 1, function(y) y[1] / y[2])


#' violin plot of features
df1 <- tibble::tibble(number_genes = apply(new_counts, 2, function(y) length(which(y != 0))),
                      mito_gene_percent = percent_mit,
                      name = "project")

p1 <- df1 %>% 
  ggplot() + 
  geom_violin(aes(x = name, y = number_genes, fill = name)) + 
  geom_jitter(aes(x = name, y = number_genes), shape = 21, size = 1, fill = "gray50", color = "black") +
  scale_fill_manual(values = c("darkred")) + 
  labs(title = NULL, x = NULL, y = "Number of genes per cell") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")

p2 <- df1 %>% 
  ggplot() + 
  geom_violin(aes(x = name, y = mito_gene_percent, fill = name)) + 
  geom_jitter(aes(x = name, y = mito_gene_percent), shape = 21, size = 1, fill = "gray50", color = "black") +
  scale_fill_manual(values = c("orange")) + 
  labs(title = NULL, x = NULL, y = "Percentage of mitochondrial reads per cell") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")


p1 + p2

# plot total read count versus mitochondrial read count
df2 <- tibble::tibble(project_name = "project",
                      cell_id = seq(1, ncol(new_counts)),
                      total_counts = num_reads_sample,
                      mitochondrial_counts = num_mit_reads, 
                      mitochondrial_percentage = percent_mit,
                      number_genes = apply(new_counts, 2, function(y) length(which(y != 0))),
                      mit_call = "bad") %>%
  dplyr::mutate(., mit_call = replace(mit_call, mitochondrial_percentage < 0.05, "good"))

 
p1 <- df2 %>% 
  ggplot() + 
  geom_violin(aes(x = project_name, y = number_genes, fill = project_name)) + 
  geom_jitter(aes(x = project_name, y = number_genes), shape = 21, size = 1, fill = "gray50", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("darkred")) + 
  labs(title = NULL, x = NULL, y = "Number of genes per cell") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")

p2 <- df2 %>% 
  ggplot() + 
  geom_violin(aes(x = project_name, y = mitochondrial_percentage, fill = project_name)) + 
  geom_jitter(aes(x = project_name, y = mitochondrial_percentage), shape = 21, size = 1, fill = "gray50", color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("orange")) + 
  labs(title = NULL, x = NULL, y = "Percentage of mitochondrial reads per cell") + 
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")


p3 <- df2 %>%
  ggplot() + 
  geom_point(aes(x = total_counts, y = mitochondrial_percentage, fill = mit_call), shape = 21, size = 1, color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0.05, color = "red", lty = 2) + 
  scale_fill_manual(values = c("darkred", "darkgreen")) + 
  labs(title = NULL, x = "Total read counts per cell", y = "Percentage of mitochondrial reads per cell") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")

p4 <- df2 %>%
  ggplot() + 
  geom_point(aes(x = total_counts, y = number_genes, fill = project_name), shape = 21, size = 1, color = "black", alpha = 0.5) +
  scale_fill_manual(values = c("darkblue")) + 
  labs(title = NULL, x = "Total read counts per cell", y = "Number of genes per cell") + 
  theme_bw() +
  theme(axis.title.x = element_text(size = 12, color = "black"), 
        axis.text.x = element_text(size = 12, color = "black"), 
        axis.title.y = element_text(size = 12, color = "black"), 
        axis.text.y = element_text(size = 12, color = "black"),
        legend.position = "none")


(p1 + p2) / (p3 + p4)

# log normalization
log_normalization <- function(count_matrix, scaling_factor = 10000) {
  if (missing(scaling_factor)) scaling_factor <- 10000
  
  total_counts <- apply(count_matrix, 2, sum)
  
  count_normalized <- sweep(count_matrix, 2, total_counts, "/")
  scaled_normalized <- count_normalized * scaling_factor
  log_normalized <- log1p(scaled_normalized)
  
  return(log_normalized)
}


# feature selection (gene filtering)
# here, the idea is to keep the genes that show high variance across all samples
# seurat uses a variance stabilizing transformatipn, modeling the variance as a function of the loess distribution for each gene
