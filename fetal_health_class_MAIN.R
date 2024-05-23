####            MATEUSZ PINDYK              ####
####                PROJEKT                 ####
####  Wnioskowanie w Warunkach Niepewnosci  ####
####      Fetal Health Classification       ####


#setwd("")
library(readr)
library(dplyr)
fetal_health <- read_csv("fetal_health.csv")  #load dataset

####  declaring names for health conditions ####
fetal_health$fetal_health <- factor(
  fetal_health$fetal_health,
  levels = c(1, 2, 3),
  labels = c("NORMAL", "SUSPECT", "PATHOLOGICAL")
)

####  deleting spare columns  ####
fetal_health <- fetal_health %>% select(-one_of(
  c(
    'prolongued_decelerations',
    'accelerations',
    'light_decelerations',
    'severe_decelerations',
    'histogram_tendency',
    "histogram_number_of_peaks",
    "histogram_number_of_zeroes",
    "histogram_mode",
    "histogram_median",
    "histogram_variance"
  )
))

#### deleting minimal and maximal values ####
nrow_fetal_bef <- nrow(fetal_health)

columns_to_clean <- c(
  "baseline value",
  "mean_value_of_short_term_variability",
  "mean_value_of_long_term_variability",
  "histogram_width",
  "histogram_min",
  "histogram_max",
  "histogram_mean"
)
# Function to remove rows with min and max values in a specific column
clean_column <- function(data, column) {
  min_val <- min(data[[column]], na.rm = TRUE)
  max_val <- max(data[[column]], na.rm = TRUE)
  
  # Remove rows with the min and max values
  data <-
    data[data[[column]] != min_val & data[[column]] != max_val, ]
  return(data)
}

# Apply the cleaning function to each column of interest
for (column in columns_to_clean) {
  fetal_health <- clean_column(fetal_health, column)
}

cat(
  "Deleted" ,
  nrow_fetal_bef - nrow(fetal_health) ,
  "rows containing minimal or maximal values."
)


####  grouping data   ####
data_group <- function(data_select) {
  min_value <- min(data_select, na.rm = TRUE)
  max_value <- max(data_select, na.rm = TRUE)
  diff_value <- max_value - min_value
  
  data_groups <- cut(
    data_select,
    breaks = c(
      min_value - 1,
      min_value + 0.2 * diff_value,
      min_value + 0.4 * diff_value,
      min_value + 0.6 * diff_value,
      min_value + 0.8 * diff_value,
      max_value + 1
    ),
    #labels = c("Very Low", "Low", "Medium", "High", "Very High"),
    include.lowest = TRUE
  )
  return(data_groups)
}

#fetal_health$baseline_value_group <- data_group(fetal_health$`baseline value`)

# List of column names to process (excluding the 12th column)
names_of_columns <- names(fetal_health)
names_of_columns <- names_of_columns[-12]

# Loop over the columns and create new grouped columns
for (i in names_of_columns) {
  new_column_name <- paste(i, 'group', sep = "_")
  fetal_health[[new_column_name]] <-
    as.factor(data_group(fetal_health[[i]])) #as factors
}

# New data frame with groups
fetal_health_grouped <- fetal_health[, c(13:23, 12)]
names(fetal_health_grouped)[1] <- "baseline_value_group"

####  independence tests ####

library(bnlearn)

# Create an empty data frame
independence_matrix <-
  data.frame(matrix(
    ncol = ncol(fetal_health_grouped),
    nrow = ncol(fetal_health_grouped)
  ))

# Assign the column names from fetal_health_grouped to the empty data frame
colnames(independence_matrix) <- colnames(fetal_health_grouped)
rownames(independence_matrix) <- colnames(fetal_health_grouped)

fetal_health_grouped <-
  as.data.frame(fetal_health_grouped) #conv to df due it only works with bnlearn tests

for (i in names(independence_matrix)) {
  for (j in names(independence_matrix)) {
    if (i != j) {
      test_result <-
        ci.test(i, j, test = "x2", data = fetal_health_grouped)
      p_value <- test_result$p.value
      independence_matrix[i, j] <- p_value
    }
  }
}

independence_matrix_less <- independence_matrix < 0.05

#### HC Network ####

library(igraph)
hc_network <- hc(fetal_health_grouped)
# Convert the Bayesian network to an igraph object
bn_igraph <- as.igraph(hc_network)

# Plot the igraph object
plot(
  bn_igraph,
  vertex.size = 30,
  vertex.label.cex = 0.8,
  edge.arrow.size = 0.5
)

# connect the missing node
hc_network <-
  hc(fetal_health_grouped, whitelist = matrix(c("fetal_movement_group", "fetal_health"), ncol = 2))

bn_igraph <- as.igraph(hc_network)

plot(
  bn_igraph,
  vertex.size = 30,
  vertex.label.cex = 0.7,
  edge.arrow.size = 0.2
)

#### IAMB Network ####
iamb_network <- iamb(fetal_health_grouped)
bn_igraph <- as.igraph(iamb_network)

plot(
  bn_igraph,
  vertex.size = 30,
  vertex.label.cex = 0.7,
  edge.arrow.size = 0.2
)
# connect the missing nodes

iamb_network <-
  iamb(fetal_health_grouped, whitelist = matrix(
    c(
      "fetal_movement_group",
      "fetal_health",
      "uterine_contractions_group" ,
      "fetal_health",
      "percentage_of_time_with_abnormal_long_term_variability_group",
      "fetal_health",
      "mean_value_of_long_term_variability_group" ,
      "mean_value_of_short_term_variability_group",
      "abnormal_short_term_variability_group",
      "fetal_health",
      "fetal_health",
      "histogram_mean_group"
    ),
    ncol = 2,
    byrow = T
  ))

iamb_network <- set.arc(iamb_network,       "histogram_width_group",
                        "histogram_max_group")
iamb_network <- set.arc(iamb_network,       "histogram_min_group",
                        "histogram_max_group")
iamb_network <- set.arc(iamb_network,       "histogram_min_group",
                        "histogram_width_group")

bn_igraph <- as.igraph(iamb_network)

plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)


#### Fast.IAMB Network ####
fast_iamb_network <- fast.iamb(fetal_health_grouped)

bn_igraph <- as.igraph(fast_iamb_network)

plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

fast_iamb_network <-
  fast.iamb(fetal_health_grouped, whitelist = matrix(
    c(
      "mean_value_of_long_term_variability_group" ,
      "mean_value_of_short_term_variability_group",
      "mean_value_of_short_term_variability_group",
      "histogram_min_group",
      "histogram_max_group",
      "histogram_mean_group",
      "baseline_value_group",
      "histogram_mean_group",
      "uterine_contractions_group" ,
      "fetal_health",
      "percentage_of_time_with_abnormal_long_term_variability_group",
      "fetal_health",
      "abnormal_short_term_variability_group",
      "fetal_health",
      "fetal_movement_group",
      "fetal_health",
      "fetal_health",
      "histogram_mean_group"
    ),
    ncol = 2,
    byrow = T
  ))

bn_igraph <- as.igraph(fast_iamb_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)


#### TABU Network ####
tabu_network <- tabu(fetal_health_grouped)

bn_igraph <- as.igraph(tabu_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

tabu_network <-
  tabu(fetal_health_grouped, whitelist = matrix(c("fetal_movement_group", "fetal_health"), ncol = 2))
bn_igraph <- as.igraph(tabu_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)


#### PC.STABLE Network ####
pc_network <- pc.stable(fetal_health_grouped)


bn_igraph <- as.igraph(pc_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

pc_network <-
  pc.stable(fetal_health_grouped, whitelist = matrix(
    c(
      "fetal_movement_group",
      "fetal_health",
      "mean_value_of_long_term_variability_group" ,
      "mean_value_of_short_term_variability_group",
      "mean_value_of_short_term_variability_group",
      "histogram_min_group",
      "uterine_contractions_group" ,
      "fetal_health",
      "percentage_of_time_with_abnormal_long_term_variability_group",
      "fetal_health",
      "abnormal_short_term_variability_group",
      "fetal_health",
      "abnormal_short_term_variability_group",
      "mean_value_of_short_term_variability_group"
      
    ),
    ncol = 2,
    byrow = T
  ))
pc_network <-
  set.arc(pc_network, "histogram_mean_group", "histogram_max_group")
pc_network <-
  set.arc(pc_network, "baseline_value_group", "histogram_mean_group")
bn_igraph <- as.igraph(pc_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

#### GS Network ####
gs_network <- gs(fetal_health_grouped)

bn_igraph <- as.igraph(gs_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

gs_network <-
  gs(fetal_health_grouped, whitelist = matrix(
    c(
      "histogram_width_group",
      "histogram_max_group",
      "baseline_value_group",
      "histogram_min_group",
      "abnormal_short_term_variability_group",
      "fetal_health",
      "percentage_of_time_with_abnormal_long_term_variability_group",
      "fetal_health",
      "histogram_max_group",
      "histogram_mean_group",
      "fetal_movement_group",
      "fetal_health",
      "uterine_contractions_group",
      "fetal_health",
      "mean_value_of_long_term_variability_group",
      "mean_value_of_short_term_variability_group",
      "abnormal_short_term_variability_group",
      "mean_value_of_short_term_variability_group"
    ),
    ncol = 2,
    byrow = T
  ))

gs_network <-
  set.arc(gs_network,
          "abnormal_short_term_variability_group",
          "baseline_value_group")

bn_igraph <- as.igraph(gs_network)
plot(
  bn_igraph,
  vertex.size = 25,
  vertex.label.cex = 0.6,
  edge.arrow.size = 0.2
)

#### Choosing network ####

# list of network names
networks <- c(
  "hc_network",
  "iamb_network",
  "fast_iamb_network",
  "tabu_network",
  "pc_network",
  "gs_network"
)

# a data frame to store the scores
network_scores <- as.data.frame(matrix(nrow = length(networks), ncol = 2))
colnames(network_scores) <- c("Network", "Score")

# a loop to compute the scores for each network
for (i in 1:length(networks)) {
  network_name <- networks[i]
  network_object <- get(network_name)  # Get the actual network object
  network_scores[i, 1] <- network_name
  network_scores[i, 2] <- score(network_object, data = fetal_health_grouped, type = "bic")
}

cat("Sieci z najwiekszym wynikiem to :\n")
network_scores[network_scores$Score == max(network_scores$Score),]

