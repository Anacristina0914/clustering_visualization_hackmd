# Clustering Visualization

[![hackmd-github-sync-badge](https://hackmd.io/hd9_lbvER5CzTDiSZbVXpw/badge)](https://hackmd.io/hd9_lbvER5CzTDiSZbVXpw)

## Table of contents
- [UMAP Description](#UMAP_Description)
    - [Number of neighbors](#Number_of_neighbors)
    - [Minimum distance](#Minimum_Distance)
    - [Distance metric](#Distance_metric)
- [Data Preprocessing](#Data_pre-processing)
- [Clustering](#Clustering)
- [UMAP implementation](#UMAP_Implementation)
- [Heatmap](#Heatmap)
- [Line plot](#Line_plot)
## UMAP_Description
When working with biological data, which often includes many variables (e.g., proteins, genes, or autoantibodies, as in our case), methods like UMAP can simplify these data into visual representations. These visualizations help us understand patterns or clusters in the data that might otherwise be hidden.

For techniques like UMAP or t-SNE, it is possible to adjust the parameters that influence how the algorithm learns the structure of the data. Three key parameters significantly impact UMAP visualizations: the number of neighbors, the minimum distance, and the distance metric.

### Number_of_neighbors
This parameter controls how UMAP balances local versus global structure in the data. It determines the size of the local neighborhood UMAP considers when learning the data’s manifold structure. Low values focus on very local structures, while higher values concentrate on a broader neighborhood, potentially losing finer details to capture a more global perspective.

### Minimum_Distance
This parameter controls how closely UMAP can pack points together in the reduced dimensional space. Lower values result in tighter, "clumpier" clusters, which can highlight local or detailed structure—useful for identifying clusters. Higher values spread the points out, focusing on preserving the overall topology of the dataset, and may highlight broader structures.

### Distance_metric
This defines how the distances between points are calculated in the input space. For binary data, we can use metrics like Hamming, Jaccard, Dice, RussellRao, or Kulczynski. These metrics assume the data is raw (i.e., not pre-processed). However, if the distances have already been calculated, one can use a pre-computed distance matrix.

In our case, since we have already used the Gower distance to measure similarity or dissimilarity between observations, we focus on optimizing the other two parameters (number of neighbors and minimum distance).

Below are some R code snippets that can be used to test different combinations of these values and evaluate which ones best represent your data. If and Rmd is used these can be pasted in a new Rmd.

```{r setup, include=FALSE}
# Import needed libraries 
library(readxl)
library(tidyverse)
library(reshape2)
library(psych)
library(cluster)
library(Rtsne)
library(umap) # to install run: install.packages(umap)
library(ggplot2)
library(fmsb) # to install run: install.packages(fmsb)
```

## Data_pre-processing
```{r proj_data}
# Define location and import data
data_dir = "/path/to/autoantibody_data"
# Where DataBase-Clusters.xlsx must include autoantibody and cluster information as columns
db_clusters = readxl::read_excel(paste0(data_dir,"/DataBase-Clusters.xlsx"), sheet = "your_excel_sheetname")
# Make another df with only subset of data
autoant_info = db_clusters
```

First 16 columns contain information on ID's, cohort and autoantibodies. Extract this information and re-format as factor. This may be different from your data so adjust the name of the autoantibodies, and the column numbers you wish to extract.
```{r subset_data}
# Make all autoantibodies 0 or 1 
binary_autoant_info <- autoant_info %>%
  mutate(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52,Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), ~ ifelse(. == "yes", 1, 0)))
# Convert 0 and 1 to factors
binary_autoant_info[c(1:16)] <- lapply(binary_autoant_info[c(1:16)], factor)
# Summary info of dataset
summary(binary_autoant_info)
```
Calculate distance matrix by:
1. Removing the individuals with all autoantibodies negative.
2. Calculate matrix using the rest.
```{r dist_matrix}
# Filter out rows that contain all 0's for autoantibodies
autoant_info_noneg <- binary_autoant_info %>%
  filter(rowSums(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52,Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), ~ . == 1)) >= 1) %>%
  # Filter only autoantibodies columns
  select(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52,Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG)

# Get only the IIDs of non-negative autoant
autoant_iid_noneg <- binary_autoant_info %>%
  filter(rowSums(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52,Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), ~ . == 1)) >= 1) %>% select(IID)
autoant_iid_noneg <- c(autoant_iid_noneg)$IID

gower_autoant_info_noneg <- daisy(autoant_info_noneg, metric="gower")
gower_autoant_info_noneg_matrix <-as.matrix(gower_autoant_info_noneg)
summary(gower_autoant_info_noneg)
```
Calculate silhouette scores as described by Lina
```{r}
 sil_width <- c(NA)
 n_clusters <- 10 # Adjust to include more clusters, e.g. 15, 20.
for(i in 2:n_clusters){
	pam_fit <- pam(gower_autoant_info_noneg,diss=TRUE,k=i)
	sil_width[i]<-pam_fit$silinfo$avg.width
	}

sillplot <- plot(1:10, sil_width,
  xlab="Number of clusters",
  ylab="Silhouette Width")
lines(1:10,sil_width)

sillplot
```
![Example_silhouette_score](https://hackmd.io/_uploads/Hktz41yaA.png)

## Clustering
```{r}
pam_fit<- pam(gower_autoant_info_noneg,diss=TRUE,k=3)
fviz_silhouette(pam_fit, palette = "jco", ggtheme = theme_classic()) #This function is from the factoextra package
```

## UMAP_Implementation
```{r umap_modular}
# Define functions to calculate and plot UMAP
plot_outpath = "/path/to/output/umap_plots" # Modify to your own path

# Define different values for the number of neighbors and min_distance (Here you can add different parameters of your choice)
nneight <- c(10L, 15L, 20L, 25L, 30L, 35L, 40L, 45L, 50L) # Discrete values
min_dist <- c(0.1, 0.25, 0.3, 0.5, 0.75, 0.99) # Value > 0 < 1

# Function to run UMAP
run_custom_umap <- function(data, neighbors, min_dist, n_epochs = 200) {
  custom.config <- umap.defaults
  custom.config$n_neighbors <- neighbors
  custom.config$min_dist <- min_dist
  custom.config$n_epochs <- n_epochs
  
  umap_result <- umap(as.matrix(data), input = "dist", config = custom.config)
  umap_coords <- as.data.frame(umap_result$layout)
  colnames(umap_coords) <- c("UMAP1", "UMAP2")
  
  return(umap_coords)
}

# Function to plot UMAP
plot_custom_umap <- function(umap_coords, clusters, neighbors, min_dist, pt.size = 0.5) {
  umap_coords$Cluster <- as.factor(clusters)
  
  ggplot_umap <- ggplot(umap_coords, aes(x = UMAP1, y = UMAP2, color = Cluster)) +
    geom_point(alpha = 0.7, size = pt.size) +
    theme_minimal() +
    labs(
      title = sprintf("UMAP Plot nneigh %s, min_dist %s", neighbors, min_dist),
      x = "UMAP Dimension 1", 
      y = "UMAP Dimension 2"
    )
  
  return(ggplot_umap)
}
```

Use the previously defined funcitons to run UMAP across different combinations of neighbors and min distance values.

```{r plot_umaps, fig.width=7, fig.width=6}
# Loop over parameter combinations and save results
umap_plots <- list()

for (neigh in nneight) { # Loop over neighbor values
  for (dist in min_dist) { # Loop over distance values
    message(sprintf("Computing UMAP for nneigh = %s, min_dist = %s", neigh, dist))
    
    umap_coords <- run_custom_umap(gower_autoant_info_noneg, neighbors = neigh, min_dist = dist)
    ggplot_umap <- plot_custom_umap(umap_coords, pam_fit$clustering, neigh, dist)
    
    umap_plots[[as.character(neigh)]][[as.character(dist)]] <- ggplot_umap
    ggsave(
      filename = sprintf("%s/umap_neigh%s_dist%s.png", plot_outpath, as.character(neigh), dist),
      plot = ggplot_umap, bg = "white"
    )
    
    print(ggplot_umap)
  }
}
```

This generates a list containing all the possible UMAP plots, it also prints it to the console in Rstudio. Here are some examples on how the UMAP can change based on the parameters:
![umap_neigh50_dist0.99](https://hackmd.io/_uploads/ByNqBy1TR.png)
![umap_neigh50_dist0.1](https://hackmd.io/_uploads/Hy_nHJ1TA.png)
![umap_neigh45_dist0.3](https://hackmd.io/_uploads/BJ9RB1kpC.png)
![umap_neigh25_dist0.25](https://hackmd.io/_uploads/HkBkUkyaC.png)
![umap_neigh15_dist0.5](https://hackmd.io/_uploads/HkzlLJkaR.png)
![umap_neigh10_dist0.25](https://hackmd.io/_uploads/r1Ax8k1a0.png)

## Heatmap

Alternatively a heatmap plot could be used to represent the distance metric (gower in this case), but the relationships of the clusters may be less obvious. See code and example below:

```{r heatmap, fig.width=8, fig.height=8}
# Cluster information
clusters <- as.factor(pam_fit$clustering)

heatmap_gowerdist <- pheatmap(gower_autoant_info_noneg, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         display_numbers = FALSE, 
         annotation_row = annotation_row, 
         main = "Gower Distance Heatmap with Cluster Annotations")
```

![cluster_heatmap](https://hackmd.io/_uploads/rJjX_kk6C.png)

## Line_plot
Pre-process data to include only data from autoantibodies and cluster information.
```{r}
#Select only autoantibody variables
binary_data_matrix <- binary_autoant_info %>% 
  filter(IID %in% autoant_iid_noneg) %>%
  select(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52,Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG)
  
# Add cluster information to matrix
binary_data_matrix$Cluster <- as.factor(pam_fit$clustering)

#Convert binary columns to numeric
binary_data_matrix <- binary_data_matrix %>%
  mutate(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52, Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), 
                ~ as.numeric(as.character(.))))

# Group by cluster and count how many individuals have positive values for each autoantobody
grouped_sums <- binary_data_matrix %>%
  group_by(Cluster) %>%  # Replace 'Cluster' with your actual group column name
  summarise(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52, Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), sum))

# Count total number of individuals per group
cluster_sizes <- binary_data_matrix %>%
  group_by(Cluster) %>%
  summarise(cluster_size = n())

# Divide the sum of 1's by the number of individuals per group and multiply by 100 for percentage
cluster_freq <- grouped_sums %>%
  left_join(cluster_sizes, by = "Cluster") %>%
  mutate(across(c(aDNAMx, aNukelosomer, Ribo_P, RNP68, RNPA, Sm, SmRNP, Ro52, Ro60Mx, SSBMx, OaCLIgG, OaCLIgM, OB2GP1IgG), 
                ~ . / cluster_size *100))
                
# Remove cluster information no longer needed
cluster_freq$Group_Size <- NULL
# Change format for plotting
cluster_freq_long <- melt(cluster_freq, id.vars = "Cluster", variable.name = "Autoantibody", value.name = "Frequency")

# Plot
ggplot(cluster_freq_long, aes(x = Autoantibody, y = Frequency, group = Cluster, color = Cluster)) +
  geom_line(size = 1.5) + 
  geom_point(size = 2) +
  theme_minimal() +
  labs(title = "Frequency of autoantibody values by Cluster", x = "Autoantibody", y = "Percentage of Positivity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

The result should look like this:
![line_plot_clusters](https://hackmd.io/_uploads/HyEujmypA.png)












