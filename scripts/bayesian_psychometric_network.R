# ==============================================================================
# STEP 1: DATA LOADING AND PSYCHOMETRIC PREPARATION
# ==============================================================================
library(MPsychoR)
library(dplyr)

#Loading the raw Rogers dataset (Adult Sample) from MPsychoR package
data("Rogers",package = "MPsychoR")

# ------------------------------------------------------------------------------
# 1.1 handling Psychometric Redundancies (QIDS-SR Filter Questions)
# ------------------------------------------------------------------------------
# According to QIDS-SR scoring guidelines, sleep and weight symptoms are 
# captured via multiple alternative questions. To prevent artificial negative 
# correlations in the network, we aggregate these items by taking the maximum 
# score (highest severity) for each psychological domain.

# Aggregate Sleep Disturbances: Max of Initial, Middle, Terminal Insomnia, and Hypersomnia
sleep_composite <- apply(Rogers[, 1:4], 1, max, na.rm = TRUE)

# Aggregate Appetite/Weight Issues: Max of Weight Decrease/Increase and Appetite Decrease/Increase
weight_composite <- apply(Rogers[, 6:9], 1, max, na.rm = TRUE)

# ------------------------------------------------------------------------------
# 1.2 Constructing the Analytical Dataset
# ------------------------------------------------------------------------------
# We select the independent single-item symptoms and bind them with our composites.
# QIDS-SR remaining items: 
#   Item 5: Mood
#   Item 10: Concentration
#   Item 11: Self-View (Guilt)
#   Item 12: Suicidal Ideation
#   Item 13: Involvement (Interest)
#   Item 14: Energy/Fatigue
#   Item 15: Psychomotor Slowing
#   Item 16: Psychomotor Agitation

depression_items <- cbind(
  Sleep_Disturbance = sleep_composite,
  Weight_Appetite = weight_composite,
  Mood            = Rogers[, 5],
  Concentration   = Rogers[, 10],
  Guilt =   Rogers[, 11],
  Suicide   = Rogers[, 12],
  Interest   = Rogers[, 13],
  Fatigue    = Rogers[, 14],
  Slowing    = Rogers[, 15],
  Agitation    = Rogers[, 16]
)

# Y-BOCS-SR items: 10 items assessing OCD symptoms (Obsessions and Compulsions)
# Columns 17 to 26 in the original Rogers dataset

ocd_items <- Rogers[, 17:26]
# Rename OCD items for better readability in the network plot

colnames(ocd_items) <- c(
  "Time_Obs", "Interf_Obs", "Distress_Obs", "Resist_Obs", "Contr_Obs",
  "Time_Comp", "Interf_Comp", "Distress_Comp", "Resist_Comp", "Contr_Comp"
)

# Combine Depression and OCD into a single clinical data frame
clinical_network_data <- as.data.frame(cbind(depression_items, ocd_items)
                                       )
# ------------------------------------------------------------------------------
# 1.3 Handling Missing Data (Listwise Deletion)
# ------------------------------------------------------------------------------
# For psychometric network models, listwise deletion is standard when the missingness
# is minimal, ensuring stable matrix estimations for both Frequentist and Bayesian samplers.

cleaned_network_data <-na.omit(clinical_network_data)

# Check the final sample size after cleaning
cat("Original sample size:", nrow(Rogers), "\n")
cat("Cleaned sample size (complete cases):", nrow(cleaned_network_data), "\n")

# ==============================================================================
# STEP 2: FREQUENTIST NETWORK ESTIMATION (EBICGLASSO)
# ==============================================================================

#Load the libraries for frequentist network analysis

library(bootnet)
library(qgraph)

# 2.1 Estimate the Network Structural Matrix using EBICglasso
# We treat the data as ordinal using 'corMethod = "cor_auto"' (highly recommended 
# for clinical Likert scales) and apply the LASSO regularization filter.

frequentist_network <- estimateNetwork(
  data = cleaned_network_data,
  default = "EBICglasso",
  corMethod = "cor_auto", # Automatically detects ordinal variables and uses polychoric correlations
  tuning = 0.5   # standard EBIC tuning parameter (gamma)
)  

# ------------------------------------------------------------------------------
# 2.2 Constructing Visual Elements and Plotting Directly on Screen
# ------------------------------------------------------------------------------

# Define 20 short abbreviations to keep the plot clean and readable
node_labels <- c(
  "DP1_Slp", "DP2_Wgh", "DP3_Moo", "DP4_Cnc", "DP5_Glt", 
  "DP6_Suc", "DP7_Int", "DP8_Fat", "DP9_Slw", "DP10_Agt",
  "OCD1_TO", "OCD2_IO", "OCD3_DO", "OCD4_RO", "OCD5_CO",
  "OCD6_TC", "OCD7_IC", "OCD8_DC", "OCD9_RC", "OCD10_CC"
)

# Group the 20 nodes into Depression (1 to 10) and OCD (11 to 20) domains
node_groups <- list(
  Depression = 1:10,
  OCD        = 11:20
)

# Plot the network using screen device (without file generation)
# Click 'Zoom' in your RStudio plot panel and expand into a perfect square!
frequentist_plot <- plot(
  frequentist_network,
  groups   = node_groups,
  color    = c("#ff9999", "#99ccff"), # Light red (Depression) and light blue (OCD)
  labels   = node_labels,             # Enforces our clean short abbreviations
  legend   = TRUE,                    
  borders  = TRUE,                    
  vsize    = 6,                       
  theme    = "classic",               
  layout   = "spring"                 
)

# Lock the spatial coordinates for the exact matching with the Bayesian model later
saved_layout <- frequentist_plot$layout

# ------------------------------------------------------------------------------
# 2.3 Calculating Node Centrality Indices
# ------------------------------------------------------------------------------

# Generate a clean centrality plot for Node Strength
# This bypasses the empty data frame issue and shows the hierarchy immediately
qgraph::centralityPlot(
  frequentist_network, 
  include = "Strength", 
  orderBy = "Strength"
)

# ==============================================================================
# STEP 3: BAYESIAN GRAPHICAL MODELING (EASYBGM)
# ==============================================================================

# Load the modern Bayesian psychometric network library
library(easybgm)


# 3.1 Estimate the Bayesian Network Structure
# We define type = "ordinal" since our symptoms are Likert-scale items.
# This will run a high-precision MCMC sampler behind the scenes.
# Note: This command might take up to a minute to finish computing.

bayesian_network <- easybgm(
  data = cleaned_network_data,
  type = "ordinal"
)

# Print the basic summary of the Bayesian estimation
print(bayesian_network)

# ------------------------------------------------------------------------------
# 3.2 Visualizing the Bayesian Network (Targeting the Parameters Matrix)
# ------------------------------------------------------------------------------

# ==============================================================================
# STEP 3: BAYESIAN NETWORKS - EXPORTING THE POSTERIOR PARAMETERS GRAPH
# ==============================================================================

# Reset any hanging or corrupted RStudio graphical devices
graphics.off()

# Initialize a high-resolution PNG device inside the working directory
# Using standard academic dimensions (2000x2000 pixels at 300 DPI)
png("bayesian_network.png", width = 2000, height = 2000, res = 300)

# Render the Bayesian Psychometric Network using the posterior parameters matrix
# Spatial coordinates are calculated via the Fruchterman-Reingold algorithm
qgraph(
  input    = bayesian_network$parameters, 
  layout   = "spring",                
  labels   = node_labels,             
  groups   = node_groups,             
  color    = c("#ff9999", "#99ccff"), 
  borders  = TRUE,
  vsize    = 6,
  theme    = "classic"
)

# Safely close the file device and write the image to disk
dev.off()


