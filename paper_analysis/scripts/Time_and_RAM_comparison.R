library(ggplot2)
library(dplyr)
library(scales)
project_ram <- function(data) {
  # Identify all unique Tool/Thread combinations
  combinations <- data %>% select(Tool, Threads) %>% distinct()
  projections <- data.frame()
  
  # The full set of nodes we want in the final benchmark
  full_node_set <- c(10000, 20000, 40000, 60000, 80000, 100000)
  
  for (i in 1:nrow(combinations)) {
    t <- combinations$Tool[i]
    thr <- combinations$Threads[i]
    
    # Filter the observed data for this specific tool/thread
    subset_df <- data %>% filter(Tool == t, Threads == thr)
    
    # Check which nodes are missing
    existing_nodes <- unique(subset_df$Nodes)
    missing_nodes <- setdiff(full_node_set, existing_nodes)
    
    # If no nodes are missing or not enough data to fit a curve (need 3+ points)
    if(length(missing_nodes) == 0 || nrow(subset_df) < 3) next
    
    # Fit quadratic model: Peak_RAM ~ Nodes + Nodes^2
    model_ram <- lm(Peak_RAM_GB ~ poly(Nodes, 2, raw = TRUE), data = subset_df)
    
    # Predict values for the missing nodes
    ram_pred <- predict(model_ram, data.frame(Nodes = missing_nodes))
    
    # Store results
    projections <- rbind(projections, data.frame(
      Tool = t,
      Threads = thr,
      Nodes = missing_nodes,
      Peak_RAM_GB = as.numeric(ram_pred),
      Type = "Projected"
    ))
  }
  return(projections)
}

project_time <- function(data) {
  combinations <- data %>% select(Tool, Threads) %>% distinct()
  projections <- data.frame()
  
  full_node_set <- c(10000, 20000, 40000, 60000, 80000, 100000)
  
  for (i in 1:nrow(combinations)) {
    t <- combinations$Tool[i]
    thr <- combinations$Threads[i]
    
    subset_df <- data %>% filter(Tool == t, Threads == thr)
    
    existing_nodes <- unique(subset_df$Nodes)
    missing_nodes <- setdiff(full_node_set, existing_nodes)
    
    if(length(missing_nodes) == 0 || nrow(subset_df) < 3) next
    
    # Fit quadratic model for time complexity
    model_time <- lm(Runtime_Sec ~ poly(Nodes, 2, raw = TRUE), data = subset_df)
    
    # Predict values for the missing nodes
    time_pred <- predict(model_time, data.frame(Nodes = missing_nodes))
    
    projections <- rbind(projections, data.frame(
      Tool = t,
      Threads = thr,
      Nodes = missing_nodes,
      Runtime_Sec = as.numeric(time_pred),
      Type = "Projected"
    ))
  }
  return(projections)
}
# 1. Load and Prepare Data (same as previous)
df <- read.delim("~/T-ChroNet/paper_analysis/data/banchmark/final_performance_comparison.csv", sep = ",")
df$Threads <- as.factor(df$Threads)



ram_projections <- project_ram(df)
full_ram_data <- df %>%
  mutate(Type = "Observed") %>%
  select(Tool, Threads, Nodes, Peak_RAM_GB, Type) %>%
  bind_rows(ram_projections) %>%
  mutate(Nodes = as.numeric(Nodes))

# Create the Single Row Plot
ram_plot_final <- ggplot(full_ram_data, aes(x = Nodes, y = Peak_RAM_GB, color = Threads, group = Threads)) +
  geom_line(linetype = "dashed", size = 0.7, alpha = 0.4) +
  geom_line(data = filter(full_ram_data, Type == "Observed"), size = 1.1) +
  geom_point(data = filter(full_ram_data, Type == "Observed"), size = 2) +
  geom_point(data = filter(full_ram_data, Type == "Projected"), size = 1.8, shape = 1) +
  
  scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3), breaks = c(10000, 60000, 100000)) +
  scale_y_continuous(labels = unit_format(unit = "GB")) +
  
  # nrow = 1 forces all plots onto a single horizontal line
  facet_wrap(~Tool, nrow = 1) + 
  
  theme_bw(base_size = 12) +
  labs( y = "Memory (GB)", x = "Nodes") +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey95"))

pdf("~/T-ChroNet/paper_analysis/data/banchmark/Supplementary_Figure_RAM.pdf",
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(ram_plot_final)
dev.off()


time_projections <- project_time(df)
full_time_data <- df %>%
  mutate(Type = "Observed") %>%
  select(Tool, Threads, Nodes, Runtime_Sec, Type) %>%
  bind_rows(time_projections) %>%
  mutate(Runtime_Min = Runtime_Sec / 60, Nodes = as.numeric(Nodes))

# Create the Single Row Plot
time_plot_final <- ggplot(full_time_data, aes(x = Nodes, y = Runtime_Min, color = Threads, group = Threads)) +
  geom_line(linetype = "dashed", size = 0.7, alpha = 0.4) +
  geom_line(data = filter(full_time_data, Type == "Observed"), size = 1.1) +
  geom_point(data = filter(full_time_data, Type == "Observed"), size = 2) +
  geom_point(data = filter(full_time_data, Type == "Projected"), size = 1.8, shape = 1) +
  
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-3), breaks = c(10000, 60000, 100000)) +
  scale_y_continuous(labels = unit_format(unit = "min")) +
  
  facet_wrap(~Tool, nrow = 1) + 
  
  theme_bw(base_size = 12) +
  labs( y = "Time (Minutes)", x = "Nodes") +
  theme(legend.position = "top", strip.background = element_rect(fill = "grey95"))


pdf("~/T-ChroNet/paper_analysis/data/banchmark/Supplementary_Figure_Time.pdf",
  height = 5, 
  width = 9,
  family = "ArialMT",
  useDingbats = FALSE)
print(time_plot_final)
dev.off()