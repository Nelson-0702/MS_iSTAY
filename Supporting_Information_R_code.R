
# ==============================================================================
# 0. Load Required Libraries
# ==============================================================================
library(tidyverse)    
library(mvtnorm)      
library(knitr)        
library(RColorBrewer)
library(gtools)       
# Check if gridExtra is installed, install if not, then load it
if(!require(gridExtra)) install.packages("gridExtra")
library(gridExtra)
library(iSTAY)
library(lmerTest)

# --- Core Function: Calculate Stability (S_q) ---
S_q <- function(data, q){
  # Convert input to matrix and remove rows with sum = 0
  data <- as.matrix(data)
  data <- data[rowSums(data) > 0, , drop = FALSE]
  
  # Get dimensions (Time steps and Number of components)
  T_step <- ncol(data)
  K <- nrow(data)
  
  if(K >= 1){
    # Calculate column sums (total abundance per time step)
    gamma_data <- colSums(data)
    Z_gamma <- sum(gamma_data)
    
    # Return NAs if total sum is 0
    if (Z_gamma == 0) {
      return(list(S_alpha=NA, S_beta=NA, max_S_beta=NA, S_gamma=NA))
    }
    
    # Calculate relative proportions (probabilities)
    p_t_gamma <- gamma_data / Z_gamma
    p_t_gamma <- p_t_gamma[p_t_gamma > 0]
    
    # Calculate Diversity Index (Hill numbers)
    # Handle the special case where q = 1 (limit approaches Shannon entropy)
    if(q != 1){
      I_gamma <- sum(p_t_gamma^q)^(1/(1-q))
    } else {
      I_gamma <- exp(-sum(p_t_gamma * log(p_t_gamma)))
    }
    
    # Normalize to get the Stability metric
    if (T_step == 1) {
      S_gamma <- 1
    } else {
      S_gamma <- (I_gamma - 1) / (T_step - 1)
    }
    
    # Return results (S_gamma is the main focus here)
    if (K == 1) {
      return(list(S_alpha=S_gamma, S_beta=0, max_S_beta=0, S_gamma=S_gamma))
    }
    
    return(list(S_gamma=S_gamma)) 
  }
}

# --- Wrapper Function: Generate Profile ---
# Calculates stability for a range of q values
calculate_stability_profile <- function(series_data, series_name, q_range) {
  data_matrix <- matrix(series_data, nrow = 1)
  # Initialize a dataframe to store results
  profile <- data.frame(q = q_range, Stability = numeric(length(q_range)), Series = series_name)
  
  # Loop through each q value
  for (i in 1:length(q_range)) {
    q_val <- q_range[i]
    stability_result <- S_q(data_matrix, q = q_val)
    profile$Stability[i] <- stability_result$S_gamma
  }
  return(profile)
}

# --- Data Preparation ---
# Define the data series
series_A_data <- c(0.999, 0.001)
series_B_data <- c(rep(0.3049, 3), rep(0.0122, 6), 0.0121)
series_C_data <- c(0.25, 0.3, 0.45)
series_D_data <- c(0.8, rep(0.02, 2), rep(0.01, 16))

# Define range for q (Order)
q_sequence <- seq(0, 2, by = 0.01)

# Calculate profiles for each series
profile_A <- calculate_stability_profile(series_A_data, "A", q_sequence)
profile_B <- calculate_stability_profile(series_B_data, "B", q_sequence)
profile_C <- calculate_stability_profile(series_C_data, "C", q_sequence)
profile_D <- calculate_stability_profile(series_D_data, "D", q_sequence)

# Combine data for plotting
data_for_plot_a <- rbind(profile_A, profile_B)
data_for_plot_b <- rbind(profile_C, profile_D)

# --- Visualization Settings ---
# Define a common theme for consistent styling
custom_theme <- theme_bw() +
  theme(
    plot.title = element_text(size = 20, face = "plain", hjust = 0),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_line(color = "grey95")
  )

# --- Plot (a): Series A & B ---
plot_a <- ggplot(data_for_plot_a, aes(x = q, y = Stability, color = Series, linetype = Series)) +
  geom_line(size = 1.2) +
  # Custom colors: A (Orange), B (Greenish)
  scale_color_manual(values = c("A" = "darkorange", "B" = "#99CC00")) +
  scale_linetype_manual(values = c("A" = "solid", "B" = "dashed")) +
  ylim(0, 1.05) +
  # Add text annotations directly on the plot
  annotate("text", x = 1.5, y = 0.08, label = "Series A", color = "darkorange", size = 7) +
  annotate("text", x = 1.5, y = 0.45, label = "Series B", color = "#99CC00", size = 7) +
  labs(
    title = "(a)",
    x = "Order q", 
    y = "Stability"
  ) +
  custom_theme

# --- Plot (b): Series C & D ---
plot_b <- ggplot(data_for_plot_b, aes(x = q, y = Stability, color = Series, linetype = Series)) +
  geom_line(size = 1.2) +
  # Custom colors: C (Red), D (Blue)
  scale_color_manual(values = c("C" = "red", "D" = "#337ab7")) +
  scale_linetype_manual(values = c("C" = "solid", "D" = "dashed")) +
  ylim(0, 1.05) +
  # Add text annotations
  annotate("text", x = 1.5, y = 0.82, label = "Series C", color = "red", size = 7) +
  annotate("text", x = 1.4, y = 0.22, label = "Series D", color = "#337ab7", size = 7) +
  labs(
    title = "(b)",
    x = "Order q", 
    y = "Stability"
  ) +
  custom_theme

# --- Final Output ---
# Arrange the two plots side-by-side
grid.arrange(plot_a, plot_b, ncol = 2)



########################################################################################################################################################################




# ==============================================================================
# Unified Analysis Script: LDM Properties and Synchrony Metrics Comparison
# 
# Description:
# This script consolidates analysis for:
# 1. The relationship between LDM (Local Dependence Measure) and correlation (rho).
# 2. Comparison of LDM against entropy-based stability metrics (Hill numbers q=1, q=2)
#    using simulation and heatmaps.
# ==============================================================================


# Set global seed for initial consistency (specific functions may override this)
set.seed(2026)

# ==============================================================================
# Reference: Appendix S3_Figure S1.R
# ==============================================================================

cat("\n--- Starting Part 1: LDM vs Rho Analysis ---\n")

# 1.1 Define Functions for Part 1
# ------------------------------------------------------------------------------

# Function: Calculate LDM for two vectors x and y
phi_lm <- function(x, y){
  vx <- var(x)
  vy <- var(y)
  cxy <- cov(x, y)
  # Denominator: (sigma_x + sigma_y)^2
  den <- vx + vy + 2 * sqrt(vx * vy) 
  # Numerator: Var(X) + Var(Y) + 2Cov(X,Y) = Var(X+Y)
  num <- vx + vy + 2 * cxy
  return(num / den)
}

# Function: Generate bivariate normal samples
# Note: Added method = "svd" to handle cases where rho = 1 or -1
sample_xy <- function(n, mu, s1, s2, rho){
  Sigma <- matrix(c(s1^2, rho*s1*s2, rho*s1*s2, s2^2), 2, 2, byrow = TRUE)
  mvtnorm::rmvnorm(n, mean = mu, sigma = Sigma, method = "svd")
}

# Function: Format ratios (display integers without decimals)
fmt_ratio <- function(x){
  ifelse(abs(x - round(x)) < 1e-10, as.character(as.integer(round(x))), as.character(x))
}

# 1.2 Simulation Setup
# ------------------------------------------------------------------------------
n  <- 10          # Sample size per iteration
B  <- 10000      # Number of bootstrap/replicates
rhos <- seq(-1, 1, by = 0.1) # Correlation sequence

sigma1_fixed <- 1
sigma2_vec <- c(1, 2, 4, 8, 16, 32, 5.8)

# Create parameter grid
sig_pairs <- tibble(
  sigma1 = sigma1_fixed,
  sigma2 = sigma2_vec
) |>
  mutate(
    ratio = sigma2 / sigma1,
    label = fmt_ratio(ratio)
  )

grid <- tidyr::crossing(rho = rhos, sig_pairs)

# 1.3 Execute Simulation
# ------------------------------------------------------------------------------
cat("Executing LDM simulation (Part 1)...\n")

sim <- purrr::map_dfr(1:B, function(b){
  grid |>
    mutate(
      ph = purrr::pmap_dbl(
        list(s1 = sigma1, s2 = sigma2, r = rho),
        function(s1, s2, r) {
          xy <- sample_xy(n, c(0,0), s1, s2, r)
          phi_lm(xy[,1], xy[,2])
        }
      ),
      rep = b
    )
})

# Aggregate results (Mean and Standard Error)
sumry <- sim |>
  group_by(ratio, label, sigma1, sigma2, rho) |>
  summarise(
    ph_mean  = mean(ph),
    ph_se    = sd(ph)/sqrt(n()),
    .groups = "drop"
  )

# 1.4 Prepare Plotting Data
# ------------------------------------------------------------------------------

# (A) Label data: Position labels at the leftmost point of each line
lab_df <- sumry |>
  group_by(ratio, label) |>
  slice_min(rho, n = 1, with_ties = FALSE) |>
  ungroup() |>
  mutate(
    ratio_chr = as.character(label),
    lab = paste0("sigma[2]/sigma[1]==", ratio_chr)
  )

# (B) Calculate intersection where LDM = 0.5 (using linear interpolation)
cross_df <- sumry |>
  group_by(ratio, label) |>
  arrange(rho, .by_group = TRUE) |>
  mutate(
    y = ph_mean,
    y0 = y - 0.5,
    y0_prev = lag(y0),
    rho_prev = lag(rho),
    y_prev = lag(y)
  ) |>
  filter(
    !is.na(y0_prev) &
      (y0 == 0 | y0_prev == 0 | (y0 * y0_prev < 0))
  ) |>
  mutate(
    rho_cross = if_else(
      y0 == 0,
      rho,
      rho_prev + (0.5 - y_prev) * (rho - rho_prev) / (y - y_prev)
    ),
    y_cross = 0.5
  ) |>
  ungroup()

# Select the rho with the smallest absolute value for the intersection (usually unique)
cross_df_one <- cross_df |>
  group_by(ratio, label) |>
  slice_min(abs(rho_cross), n = 1, with_ties = FALSE) |>
  ungroup()

# 1.5 Generate Plot (Simulation Result)
# ------------------------------------------------------------------------------
x_min <- min(sumry$rho, na.rm = TRUE)
x_max <- max(sumry$rho, na.rm = TRUE)

p1 <- ggplot(sumry, aes(x = rho, y = ph_mean, group = label, color = label)) +
  # Dashed line at y=0.5
  geom_segment(
    aes(x = -1, xend = 1, y = 0.5, yend = 0.5), 
    linetype = "dashed", 
    linewidth = 0.7, 
    color = "gray50",
    inherit.aes = FALSE
  ) +
  geom_line(linewidth = 1) +
  # Red intersection points
  geom_point(
    data = cross_df_one,
    aes(x = rho_cross, y = y_cross),
    inherit.aes = FALSE,
    color = "red",
    size = 2.5
  ) +
  # Line labels
  geom_text(
    data = lab_df,
    aes(label = lab),
    nudge_x = -0.05,
    hjust = 1,
    size = 3.5,
    show.legend = FALSE,
    parse = TRUE
  ) +
  coord_cartesian(
    xlim = c(x_min - 0.25, x_max),
    ylim = c(0, 1)
  ) +
  labs(
    title = expression(paste("Simulation: ", phi[LM], " vs. ", rho)),
    x = expression(rho),
    y = expression(phi[LM])
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

print(p1)


#################################################################################################################


iSTAY_Multiple = function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, 
                           end_T = NULL) 
{
  NA_num <- sum(is.na(data))
  if (NA_num != 0) {
    stop("There are some NA in the data.")
  }
  if (sum(which(order.q <= 0))) {
    stop("Order q need larger than 0.")
  }
  if (Alltime != TRUE) {
    if (is.null(start_T) | is.null(end_T)) {
      stop("Need to set the length of time series for calculating.")
    }
    if ((start_T >= end_T) | length(start_T) != 1 | length(end_T) != 
        1) {
      stop("Starting and ending time need to be a number, and ending time needs larger than starting time.")
    }
  }
  Stabillity_Multiple <- function(ZZ, q) {
    K <- ncol(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    
    ww <- 1/nrow(ZZ)
    
    if (nrow(ZZ) != 1) p_i <- as.data.frame(apply(ZZ, 2, function(w) w/z_iplus)) else p_i = matrix(apply(ZZ, 2, function(w) w/z_iplus), nrow = 1)
    ZZ = p_i
    z_plusk <- apply(ZZ, 2, sum)
    z_plusplus <- sum(ZZ)
    p_pool <- z_plusk/z_plusplus
    # ww <- z_iplus/z_plusplus
    
    if (q == 1) {
      p_i_new <- p_i
      p_pool_new <- p_pool
      alpha <- ( exp(-sum(ww * p_i_new[p_i_new > 0] * log(p_i_new[p_i_new > 0])) ) - 1)/(K - 1)
      gamma <- (exp(-sum(p_pool[p_pool > 0] * log(p_pool_new[p_pool_new > 
                                                               0]))) - 1)/(K - 1)
    }
    else {
      alpha <- (sum(ww * (p_i)^q)^(1/(1 - q)) - 1)/(K - 1)
      gamma <- (sum(p_pool^q)^(1/(1 - q)) - 1)/(K - 1)
    }
    return(c(gamma, alpha, (gamma/alpha), (gamma - alpha)))
  }
  # Synchrony <- function(ZZ, q) {
  #   M <- nrow(ZZ)
  #   K <- ncol(ZZ)
  #   if (M == 1) {
  #     value <- 1
  #   }
  #   else {
  #     z_iplus <- apply(ZZ, 1, sum)
  #     if (length(which(z_iplus == 0)) != 0) {
  #       ZZ <- ZZ[-which(z_iplus == 0), ]
  #     }
  #     ZZ <- as.matrix(ZZ)
  #     z_iplus <- apply(ZZ, 1, sum)
  #     if (length(z_iplus) <= 1) {
  #       value <- NA
  #     }
  #     else {
  #       z_plusk <- apply(ZZ, 2, sum)
  #       z_plusplus <- sum(ZZ)
  #       p_i <- apply(ZZ, 2, function(w) w/z_iplus)
  #       ww <- z_iplus/z_plusplus
  #       if (q == 1) {
  #         J <- exp(-sum(ZZ[ZZ > 0]/z_plusplus * log(ZZ[ZZ > 
  #                                                        0]/z_plusplus)))
  #         J <- ifelse(J < K, J, J/M + K * (M - 1)/M)
  #         pool <- z_plusk/z_plusplus
  #         pool[which(pool == 0)] <- 10^(-15)
  #         G <- exp(-sum(pool * log(pool)))
  #         A <- exp(-sum(ZZ[ZZ > 0]/z_plusplus * log(p_i[p_i > 
  #                                                         0])))
  #         value <- (J - G)/(J - A)
  #       }
  #       else {
  #         J <- sum((ZZ/z_plusplus)^q)^(1/(1 - q))
  #         J <- ifelse(J < K, J, J/M + K * (M - 1)/M)
  #         G <- sum((z_plusk/z_plusplus)^q)^(1/(1 - q))
  #         A <- sum(apply(ZZ, 2, function(w) (w/z_iplus)^q * 
  #                          ww))^(1/(1 - q))
  #         value <- (J - G)/(J - A)
  #       }
  #     }
  #   }
  #   return(value)
  # }
  
  Synchrony <- function(ZZ, q) {
    
    
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    
    if(sum(rowSums(ZZ)>0)!=1  ){
      
      K = nrow(ZZ)
      Time = ncol(ZZ)
      beta = Stabillity_Multiple(ZZ,q)[4]
      gamma = Stabillity_Multiple(ZZ,q)[1]
      
      beta_max = (1 - 1/K)/(Time - 1)+(1 - 1/K)*gamma
      1 - beta/beta_max
      
    }else{
      
      1
      
    } 
    
    
  }
  
  if (is.data.frame(data) | is.matrix(data)) {
    if (Alltime == TRUE) {
      subdata <- data
    }
    else {
      subdata <- data[, c(start_T:end_T)]
    }
    out <- as.matrix(sapply(order.q, function(qq) c(Stabillity_Multiple(subdata, 
                                                                        q = qq), Synchrony(subdata, q = qq))))
    result <- data.frame(Site = rep(1, length(order.q)), 
                         Order_q = order.q, t(out))
    colnames(result)[3:7] <- c("Stab_Gamma", "Stab_Alpha", 
                               "Stab_Beta_multiplicative", "Stab_Beta_additive", 
                               "Synchrony")
    result <- result[, -5]
  }
  else {
    out <- lapply(order.q, function(qq) {
      cal <- lapply(data, function(ZZ) {
        if (Alltime == T) {
          subZZ <- ZZ
        }
        else {
          subZZ <- ZZ[, c(start_T:end_T)]
        }
        outout <- c(Stabillity_Multiple(subZZ, q = qq), 
                    Synchrony(subZZ, q = qq))
        result <- data.frame(Order_q = qq, t(outout))
        colnames(result)[2:6] <- c("Stab_Gamma", "Stab_Alpha", 
                                   "Stab_Beta_multiplicative", "Stab_Beta_additive", 
                                   "Synchrony")
        result <- result[, -4]
        return(result)
      })
      cal2 <- do.call(rbind, cal)
      calcal <- data.frame(Site = names(data), cal2)
      return(calcal)
    })
    result <- do.call(rbind, out)
  }
  colnames(result) <- c("Dataset", "Order_q", "Gamma", "Alpha", 
                        "Beta", "Synchrony")
  return(result)
}




# Jena LMM ----------------------------------------------------------------

data("Data_Jena_hierarchical_structure")
data("Data_Jena_462_populations")
data("Data_Jena_20_metacommunities")
data("Data_Jena_76_metapopulations")

multi_gamma_diffk <- function(data, blocksowndiv) {
  
  # List of functions to iterate
  stability_functions <- list(
    E3_Routledge = iSTAY_Multiple
  )
  
  # Initialize empty data.frames for final results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  combined_plotdata <- data.frame()
  
  for (method_name in names(stability_functions)) {
    method_func <- stability_functions[[method_name]]
    
    # Stability data for current method
    # stability_data <- method_func(data, order.q = c(0.3, 0.7, 1, 2))
    stability_data <- method_func(data, order.q = c(0.5, 1, 2))
    
    # Base plotdata for this method
    plotdata <- data.frame(
      block = rep(blocksowndiv$block, 3),
      sowndiv = rep(blocksowndiv$sowndiv, 3),
      gamma = c(stability_data[, 3]),
      version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(blocksowndiv))),
      Class = gsub("_.*", "", method_name),
      Method = gsub(".*_", "", method_name)
    )
    
    # Loop over each version to fit LMMs separately
    for (v in levels(plotdata$version)) {
      subdata <- filter(plotdata, version == v)
      
      # Linear mixed-effects model
      model <- lmerTest::lmer(gamma ~ 1 + sowndiv + (1 + sowndiv | block), data = subdata)
      summary_model <- summary(model)
      
      # Add predictions and significance
      subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
      subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
      subdata$sign <- ifelse(summary_model$coefficients[2, 1] > 0, "positive", "negative")
      
      # Add slope and intercepts for blocks
      slopes <- coef(model)$block
      slope_intercept <- data.frame(
        block = rownames(slopes),
        intercept = slopes[, 1],
        slope = slopes[, 2],
        x_min = tapply(subdata$sowndiv, subdata$block, min),
        x_max = tapply(subdata$sowndiv, subdata$block, max),
        version = rep(v,4),
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Add slope text for visualization
      slope_text_metric <- data.frame(
        block = c(rownames(slopes), "Total"),
        slope_text = c(paste0("slope = ", sprintf("%.4f", slopes[, 2])), paste0("slpoe = ", sprintf("%.4f", summary_model$coefficients[2, 1]))),
        version = v,
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Append to combined results
      combined_total_fit <- rbind(combined_total_fit, subdata)
      combined_part_fit <- rbind(combined_part_fit, slope_intercept)
      combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
    }
    
    # Append the overall plotdata for this method
    combined_plotdata <- rbind(combined_plotdata, plotdata)
  }
  
  # Return combined results as a list of 4 data.frames
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = combined_plotdata
  ))
}

multi_alpha_diffk <- function(data, blocksowndiv) {
  
  # List of functions to iterate
  stability_functions <- list(
    E3_Routledge = iSTAY_Multiple
  )
  
  # Initialize empty data.frames for final results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  combined_plotdata <- data.frame()
  
  for (method_name in names(stability_functions)) {
    method_func <- stability_functions[[method_name]]
    
    # Stability data for current method
    stability_data <- method_func(data, order.q = c(0.5, 1, 2))
    
    # Base plotdata for this method
    plotdata <- data.frame(
      block = rep(blocksowndiv$block, 3),
      sowndiv = rep(blocksowndiv$sowndiv, 3),
      alpha = c(stability_data[, 4]),
      version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(blocksowndiv))),
      Class = gsub("_.*", "", method_name),
      Method = gsub(".*_", "", method_name)
    )
    
    # Loop over each version to fit LMMs separately
    for (v in levels(plotdata$version)) {
      subdata <- filter(plotdata, version == v)
      
      # Linear mixed-effects model
      model <- lmerTest::lmer(alpha ~ 1 + sowndiv + (1 + sowndiv | block), data = subdata)
      summary_model <- summary(model)
      
      # Add predictions and significance
      subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
      subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
      subdata$sign <- ifelse(summary_model$coefficients[2, 1] > 0, "positive", "negative")
      
      # Add slope and intercepts for blocks
      slopes <- coef(model)$block
      slope_intercept <- data.frame(
        block = rownames(slopes),
        intercept = slopes[, 1],
        slope = slopes[, 2],
        x_min = tapply(subdata$sowndiv, subdata$block, min),
        x_max = tapply(subdata$sowndiv, subdata$block, max),
        version = rep(v, 4),
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Add slope text for visualization
      slope_text_metric <- data.frame(
        block = c(rownames(slopes), "Total"),
        slope_text = c(paste0("slope = ", sprintf("%.4f", slopes[, 2])), paste0("slpoe = ", sprintf("%.4f", summary_model$coefficients[2, 1]))),
        version = v,
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Append to combined results
      combined_total_fit <- rbind(combined_total_fit, subdata)
      combined_part_fit <- rbind(combined_part_fit, slope_intercept)
      combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
    }
    
    # Append the overall plotdata for this method
    combined_plotdata <- rbind(combined_plotdata, plotdata)
  }
  
  # Return combined results as a list of 4 data.frames
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = combined_plotdata
  ))
}

multi_beta_diffk <- function(data, blocksowndiv) {
  
  # List of functions to iterate
  stability_functions <- list(
    E3_Routledge = iSTAY_Multiple
  )
  
  # Initialize empty data.frames for final results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  combined_plotdata <- data.frame()
  
  for (method_name in names(stability_functions)) {
    method_func <- stability_functions[[method_name]]
    
    # Stability data for current method
    stability_data <- method_func(data, order.q = c(0.5, 1, 2))
    
    # Base plotdata for this method
    plotdata <- data.frame(
      block = rep(blocksowndiv$block, 3),
      sowndiv = rep(blocksowndiv$sowndiv, 3),
      beta = c(stability_data[, 5]),
      version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(blocksowndiv))),
      Class = gsub("_.*", "", method_name),
      Method = gsub(".*_", "", method_name)
    )
    
    # Loop over each version to fit LMMs separately
    for (v in levels(plotdata$version)) {
      subdata <- filter(plotdata, version == v)
      
      # Linear mixed-effects model
      model <- lmerTest::lmer(beta ~ 1 + sowndiv + (1 + sowndiv | block), data = subdata)
      summary_model <- summary(model)
      
      # Add predictions and significance
      subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
      subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
      subdata$sign <- ifelse(summary_model$coefficients[2, 1] > 0, "positive", "negative")
      
      # Add slope and intercepts for blocks
      slopes <- coef(model)$block
      slope_intercept <- data.frame(
        block = rownames(slopes),
        intercept = slopes[, 1],
        slope = slopes[, 2],
        x_min = tapply(subdata$sowndiv, subdata$block, min),
        x_max = tapply(subdata$sowndiv, subdata$block, max),
        version = rep(v, 4),
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Add slope text for visualization
      slope_text_metric <- data.frame(
        block = c(rownames(slopes), "Total"),
        slope_text = c(paste0("slope = ", sprintf("%.4f", slopes[, 2])), paste0("slpoe = ", sprintf("%.4f", summary_model$coefficients[2, 1]))),
        version = v,
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Append to combined results
      combined_total_fit <- rbind(combined_total_fit, subdata)
      combined_part_fit <- rbind(combined_part_fit, slope_intercept)
      combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
    }
    
    # Append the overall plotdata for this method
    combined_plotdata <- rbind(combined_plotdata, plotdata)
  }
  
  # Return combined results as a list of 4 data.frames
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = combined_plotdata
  ))
}

multi_syn_diffk <- function(data, blocksowndiv) {
  
  # List of functions to iterate
  stability_functions <- list(
    E3_Routledge = iSTAY_Multiple
  )
  
  # Initialize empty data.frames for final results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  combined_plotdata <- data.frame()
  
  for (method_name in names(stability_functions)) {
    method_func <- stability_functions[[method_name]]
    
    # Stability data for current method
    stability_data <- method_func(data, order.q = c(0.5, 1, 2))
    
    # Base plotdata for this method
    plotdata <- data.frame(
      block = rep(blocksowndiv$block, 3),
      sowndiv = rep(blocksowndiv$sowndiv, 3),
      Synchrony = c(stability_data[, 6]),
      version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(blocksowndiv))),
      Class = gsub("_.*", "", method_name),
      Method = gsub(".*_", "", method_name)
    )
    
    # Loop over each version to fit LMMs separately
    for (v in levels(plotdata$version)) {
      subdata <- filter(plotdata, version == v)
      
      # Linear mixed-effects model
      model <- lmerTest::lmer(Synchrony ~ 1 + sowndiv + (1 + sowndiv | block), data = subdata)
      summary_model <- summary(model)
      
      # Add predictions and significance
      subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
      subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
      subdata$sign <- ifelse(summary_model$coefficients[2, 1] > 0, "positive", "negative")
      
      # Add slope and intercepts for blocks
      slopes <- coef(model)$block
      slope_intercept <- data.frame(
        block = rownames(slopes),
        intercept = slopes[, 1],
        slope = slopes[, 2],
        x_min = tapply(subdata$sowndiv, subdata$block, min),
        x_max = tapply(subdata$sowndiv, subdata$block, max),
        version = rep(v,4),
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Add slope text for visualization
      slope_text_metric <- data.frame(
        block = c(rownames(slopes), "Total"),
        slope_text = c(paste0("slope = ", sprintf("%.4f", slopes[, 2])), paste0("slpoe = ", sprintf("%.4f", summary_model$coefficients[2, 1]))),
        version = v,
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Append to combined results
      combined_total_fit <- rbind(combined_total_fit, subdata)
      combined_part_fit <- rbind(combined_part_fit, slope_intercept)
      combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
    }
    
    # Append the overall plotdata for this method
    combined_plotdata <- rbind(combined_plotdata, plotdata)
  }
  
  # Return combined results as a list of 4 data.frames
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = combined_plotdata
  ))
}

multi_asyn_diffk <- function(data, blocksowndiv) {
  
  # List of functions to iterate
  stability_functions <- list(
    E3_Routledge = iSTAY_Multiple
  )
  
  # Initialize empty data.frames for final results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  combined_plotdata <- data.frame()
  
  for (method_name in names(stability_functions)) {
    method_func <- stability_functions[[method_name]]
    
    # Stability data for current method
    stability_data <- method_func(data, order.q = c(0.5, 1, 2))
    
    # Base plotdata for this method
    plotdata <- data.frame(
      block = rep(blocksowndiv$block, 3),
      sowndiv = rep(blocksowndiv$sowndiv, 3),
      Asynchrony = 1 - c(stability_data[, 6]),
      version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(blocksowndiv))),
      Class = gsub("_.*", "", method_name),
      Method = gsub(".*_", "", method_name)
    )
    
    # Loop over each version to fit LMMs separately
    for (v in levels(plotdata$version)) {
      subdata <- filter(plotdata, version == v)
      
      # Linear mixed-effects model
      model <- lmerTest::lmer(Asynchrony ~ 1 + sowndiv + (1 + sowndiv | block), data = subdata)
      summary_model <- summary(model)
      
      # Add predictions and significance
      subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
      subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
      subdata$sign <- ifelse(summary_model$coefficients[2, 1] > 0, "positive", "negative")
      
      # Add slope and intercepts for blocks
      slopes <- coef(model)$block
      slope_intercept <- data.frame(
        block = rownames(slopes),
        intercept = slopes[, 1],
        slope = slopes[, 2],
        x_min = tapply(subdata$sowndiv, subdata$block, min),
        x_max = tapply(subdata$sowndiv, subdata$block, max),
        version = rep(v,4),
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Add slope text for visualization
      slope_text_metric <- data.frame(
        block = c(rownames(slopes), "Total"),
        slope_text = c(paste0("slope = ", sprintf("%.4f", slopes[, 2])), paste0("slpoe = ", sprintf("%.4f", summary_model$coefficients[2, 1]))),
        version = v,
        Class = gsub("_.*", "", method_name),
        Method = gsub(".*_", "", method_name)
      )
      
      # Append to combined results
      combined_total_fit <- rbind(combined_total_fit, subdata)
      combined_part_fit <- rbind(combined_part_fit, slope_intercept)
      combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
    }
    
    # Append the overall plotdata for this method
    combined_plotdata <- rbind(combined_plotdata, plotdata)
  }
  
  # Return combined results as a list of 4 data.frames
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = combined_plotdata
  ))
}


# 10yr window -------------------------------------------------------------

struct_plot <- data.frame(plot = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 1),
                          block = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 2),
                          sowndiv = as.numeric(map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 3)))

blocksowndiv_spe <- struct_plot |> select(block, sowndiv)
colnames(blocksowndiv_spe) <- c("block", "sowndiv")
blocksowndiv_spe$sowndiv <- log2(as.numeric(blocksowndiv_spe$sowndiv))

# 10yr 3set ---------------------------------------------------------------

struct_plot <- data.frame(plot = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 1),
                          block = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 2),
                          sowndiv = as.numeric(map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 3)))

blocksowndiv_spe <- struct_plot |> select(block, sowndiv)
colnames(blocksowndiv_spe) <- c("block", "sowndiv")
blocksowndiv_spe$sowndiv <- log2(as.numeric(blocksowndiv_spe$sowndiv))

# Community Stability(76) ----------------------------------------------------

valid_start_years <- c(2005, 2010, 2015)
year_windows <- lapply(valid_start_years, function(start_year) {
  years <- start_year:(start_year + 9)
  if (start_year == 2003) {
    years <- c(2003, 2005:2013)
  } else if (2004 %in% years) {
    return(NULL)
  }
  return(as.character(years))
})
year_windows <- Filter(Negate(is.null), year_windows)

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
    sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
    sub_df
  })
  names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
  out <- multi_gamma_diffk(Data_Jena_76_metapopulations_10, blocksowndiv_spe)
  return(out)
})

year_labels <- sapply(year_windows, function(yrs){
  paste0(min(yrs), "-", max(yrs))
})

plotdata_text1 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$slope_text |>
    filter(block == "Total") |>
    mutate(year = year_labels[i])
}))

plotdata1 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$total_fit |>
    mutate(year = year_labels[i])
}))

plotdata1$significance <- factor(plotdata1$significance, levels = c("significant", "non-significant"))

if(!any(plotdata1$significance == "non-significant")) {
  dummy <- plotdata1[1, ]  # 以第一列為範例
  dummy$sowndiv <- NA      # 可以設定為 NA 或適當的數值
  dummy$predicted <- NA
  dummy$gamma <- NA
  dummy$significance <- "non-significant"
  plotdata1 <- rbind(plotdata1, dummy)
}

ggplot() +
  geom_point(data = plotdata1,
             aes(x = sowndiv, y = gamma, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = plotdata1,
            aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
  geom_text(data = plotdata_text1,
            aes(x = 2.5, y = 0.25, label = slope_text, color = year,
                hjust = rep(-c(0,0,0), 3),
                vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("#EA0000","#0066CC","purple")) +
  scale_linetype_manual(values = c("solid", "dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), label = c(1, 2, 4, 8, 16)) +
  facet_grid(~ version) +
  guides(linetype = guide_legend(keywidth = 3.1)) +
  labs(linetype="", x = "Number of species (log2 scale)", y = "Community stability") + 
  theme_bw() +
  theme(legend.position = "bottom", 
        axis.title = element_text(size = 13),
        legend.box.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        axis.text = element_text(size = 16),
        text = element_text(size = 14), 
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14))

# ggsave("Community stability 0512.png", width = 9.5, height = 5)


# Population Stability(76) -----------------------------------------------------

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
    sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
    sub_df
  })
  names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
  out <- multi_alpha_diffk(Data_Jena_76_metapopulations_10, blocksowndiv_spe)
  return(out)
})

plotdata_text2 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$slope_text |>
    filter(block == "Total") |>
    mutate(year = year_labels[i])
}))

plotdata2 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$total_fit |>
    mutate(year = year_labels[i])
}))

plotdata2 <- plotdata2 |> mutate(sign = factor(sign, levels = c("positive", "negative")))
plotdata2$significance <- factor(plotdata2$significance, levels = c("significant", "non-significant"))

ggplot()+
  geom_point(data = plotdata2,
             aes(x = sowndiv, y = alpha, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = plotdata2, aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
  geom_text(data = plotdata_text2,
            aes(x = 2.5, y = 0.8, label = slope_text, color = year,
                hjust = rep(-c(0,0,0), 3),
                vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("#EA0000","#0066CC","purple")) +
  scale_linetype_manual(values = c("solid","dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels = c(1, 2, 4, 8, 16)) +
  facet_grid(~ version) +
  guides(linetype = guide_legend(keywidth = 3.1)) +
  labs(linetype = "", x = "Number of species (log2 scale)", y = "Population stability") + theme_bw()+
  theme(legend.position="bottom", 
        # axis.text=element_text(size=10), 
        axis.title=element_text(size=13),
        legend.box.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        text = element_text(size = 14), 
        axis.text = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14))

# ggsave("Population stability 0512.png", width = 9.5, height = 5)

# Beta(76) -----------------------------------------------------------

names <- struct_plot |>
  filter(sowndiv != 1) |>
  mutate(names = paste0(plot, "_", block, "_", sowndiv)) |>
  select(names) |> unlist()

blocksowndiv_spe2 <- struct_plot |> filter(sowndiv != 1) |> select(block, sowndiv)
colnames(blocksowndiv_spe2) <- c("block", "sowndiv")
blocksowndiv_spe2$sowndiv <- log2(as.numeric(blocksowndiv_spe2$sowndiv))


# Jena_10yr_window_result <- lapply(year_windows, function(years){
#   Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
#     sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
#     sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
#     sub_df
#   })
#   names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
#   out <- multi_beta_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
#   return(out)
# })
# 
# plotdata_text3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |>
#     mutate(year = year_labels[i])
# }))
# 
# plotdata3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$total_fit |>
#     mutate(year = year_labels[i])
# }))
# 
# plotdata3 <- plotdata3 |> mutate(sign = factor(sign, levels = c("positive", "negative")))
# plotdata3$significance <- factor(plotdata3$significance, levels = c("significant", "non-significant"))
# 
# 
# ggplot() +
#   geom_point(data = plotdata3, 
#              aes(x = sowndiv, y = beta, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
#   geom_line(data = plotdata3, aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
#   geom_text(data = plotdata_text3,
#             aes(x = 2.55, y = 0.545, label = slope_text, color = year,
#                 hjust = rep(-c(0,0,0), 3),
#                 vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
#   scale_color_manual(values = c("#EA0000","#0066CC","purple")) +
#   scale_linetype_manual(values = c("solid","dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
#   scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(2, 4, 8, 16)) +
#   facet_grid(~ version) +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
#   labs(linetype = "", x = "Number of species (log2 scale)", y = "Beta stability") + 
#   theme_bw() +
#   theme(legend.position = "bottom", 
#         # axis.text=element_text(size=10),
#         axis.title=element_text(size=13),
#         legend.box.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
#         legend.text = element_text(size = 12),
#         legend.title = element_blank(), 
#         # legend.margin = margin(0, 0, 0, 0),
#         # legend.box.margin = margin(-10, -10, -5, -10), 
#         text = element_text(size = 14), 
#         # plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
#         axis.text = element_text(size = 16),
#         strip.text.x = element_text(size = 14),
#         strip.text.y = element_text(size = 14),
#         axis.title.x = element_text(hjust = 0.5, size = 14),
#         axis.title.y = element_text(hjust = 0.5, size = 14))

# ggsave("Population Beta 0512.png", width = 9.5, height = 5)


# Synchroy (76) -----------------------------------------------------------

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      
    sub_df[sub_df == 0] <- 1e-15             
    sub_df
  })
  names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
  out <- multi_syn_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
  return(out)
})

plotdata_text3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$slope_text |>
    filter(block == "Total") |>
    mutate(year = year_labels[i])
}))

plotdata3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$total_fit |>
    mutate(year = year_labels[i])
}))

plotdata3 <- plotdata3 |> mutate(sign = factor(sign, levels = c("positive", "negative")))
plotdata3$significance <- factor(plotdata3$significance, levels = c("significant", "non-significant"))


ggplot() +
  geom_point(data = plotdata3, 
             aes(x = sowndiv, y = Synchrony, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = plotdata3, aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
  geom_text(data = plotdata_text3,
            aes(x = 2.5, y = 0.135, label = slope_text, color = year,
                hjust = rep(-c(0,0,0), 3),
                vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("#EA0000","#0066CC","purple")) +
  scale_linetype_manual(values = c("solid","dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(2, 4, 8, 16)) +
  facet_grid(~ version) +
  guides(linetype = guide_legend(keywidth = 3.1)) +
  labs(linetype = "", x = "Number of species (log2 scale)", y = "Population synchrony") + 
  theme_bw() +
  theme(legend.position = "bottom", 
        # axis.text=element_text(size=10),
        axis.title=element_text(size=13),
        legend.box.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin = margin(-10, -10, -5, -10), 
        text = element_text(size = 14), 
        # plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        axis.text = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14))

# ggsave("Population Synchrony 0512.png", width = 9.5, height = 5)

# Asynchrony(76) ----------------------------------------------------------

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
    sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
    sub_df
  })
  names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
  out <- multi_asyn_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
  return(out)
})

plotdata_text3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$slope_text |>
    filter(block == "Total") |>
    mutate(year = year_labels[i])
}))

plotdata3 <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
  Jena_10yr_window_result[[i]]$total_fit |>
    mutate(year = year_labels[i])
}))

plotdata3 <- plotdata3 |> mutate(sign = factor(sign, levels = c("positive", "negative")))
plotdata3$significance <- factor(plotdata3$significance, levels = c("significant", "non-significant"))


ggplot() +
  geom_point(data = plotdata3, 
             aes(x = sowndiv, y = Asynchrony, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = plotdata3, aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
  geom_text(data = plotdata_text3,
            aes(x = 3, y = 0.115, label = slope_text, color = year,
                hjust = rep(-c(0,0,0), 3),
                vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("#EA0000","#0066CC","purple")) +
  scale_linetype_manual(values = c("solid","dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(2, 4, 8, 16)) +
  facet_grid(~ version) +
  guides(linetype = guide_legend(keywidth = 3.1)) +
  labs(linetype = "", x = "Number of species (log2 scale)", y = "Population asynchrony") + 
  theme_bw() +
  theme(legend.position = "bottom", 
        # axis.text=element_text(size=10),
        axis.title=element_text(size=13),
        legend.box.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.text = element_text(size = 12),
        legend.title = element_blank(), 
        # legend.margin = margin(0, 0, 0, 0),
        # legend.box.margin = margin(-10, -10, -5, -10), 
        text = element_text(size = 14), 
        # plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        axis.text = element_text(size = 16),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.title.x = element_text(hjust = 0.5, size = 14),
        axis.title.y = element_text(hjust = 0.5, size = 14))

# ggsave("Population Asynchrony 0512.png", width = 9.5, height = 5)



# Weight vs Stab ----------------------------------------------------------

Jena_tidy <- Data_Jena_462_populations |> 
  mutate(rowname = row.names(Data_Jena_462_populations),
         block = str_sub(rowname, 1, 2),
         number_of_species = as.numeric(sub(".*_([0-9]+)_.*", "\\1", rowname)),
         plot = str_sub(sub(".*_(B[1234]A[0-9]{2})_.*", "\\1", rowname), 3, 5),
         species = str_sub(rowname, -7, -1),
         label = paste(block, plot, number_of_species, sep = "_"))


gg_alpha_wrt_weight_group_by_species <- function(data = Jena_tidy, q = 1, threshold = 1 / 6){
  
  # 這個函數可以輸出每個species的計算stab的斜率和預測值
  fit_slope <- function(q){
    
    temp <- data |> 
      mutate(total_biomass = rowSums(across(contains("20")))) |> 
      group_by(label) |> 
      mutate(weight = total_biomass / sum(total_biomass),
             Stab = iSTAY_Single(data = across(contains("20")), order.q = q)$Stab,
             number_of_species = factor(number_of_species),
             iSTAY_Multiple(data = across(contains("20")), order.q = q)) |> 
      ungroup() |> 
      group_by(number_of_species) 
    
    temp <- split(temp, temp$number_of_species)
    
    result1 <- lapply(temp, function(x) {
      
      summary_model <- lm(Stab ~ weight, data = x) |> summary()
      data.frame(number_of_species = unique(x$number_of_species),
                 Slope = tryCatch({
                   summary_model$coefficients["weight", "Estimate"]
                 }, error = function(e) {
                   NA 
                 }),
                 p.value = tryCatch({
                   summary_model$coefficients["weight", "Pr(>|t|)"]
                 }, error = function(e) {
                   NA 
                 })) |> 
        mutate(Significant = if_else(p.value < 0.05, "Yes", "No"))
      
    }) |> map_df(~.x)
    
    result2 <- lapply(temp, function(x) {
      
      model <- lm(Stab ~ weight, data = x)
      predictions <- predict(model, newdata = x |> select(weight, Stab))
      data.frame(number_of_species = unique(x$number_of_species),
                 weight = x$weight,
                 predicted = predictions)
      
    }) |> map_df(~.x) |> 
      remove_rownames()
    
    list(result1 = result1, result2 = result2)
    
  }
  

  
  temp <- data |> 
    mutate(total_biomass = rowSums(across(contains("20")))) |> 
    group_by(label) |> 
    mutate(weight = total_biomass / sum(total_biomass),
           Stab = iSTAY_Single(data = across(contains("20")), order.q = q)$Stab,
           number_of_species = factor(number_of_species),
           iSTAY_Multiple(data = across(contains("20")), order.q = q)) |>
    ungroup() |> 
    group_by(number_of_species) |> 
    mutate(mean_Alpha = mean(Alpha),
           number_of_species = factor(number_of_species, levels = c(1, 2, 4, 8, 16)))
  
  temp_ono <- temp |> 
    filter(species == "Ono.vic") |> 
    group_by(number_of_species) |> 
    mutate(mean_ono = mean(Stab),
           text = paste0("Mean Ono.Vic = ", sprintf("%.4f", mean_ono)))
  
  temp |> 
    ggplot() +
    geom_point(aes(x = weight, y = Stab), size = 1, alpha = 0.5) +
    geom_hline(data = temp_ono, aes(yintercept = mean_ono), linetype = "dashed", linewidth = 1, color = "#F8766D") +
    # geom_point(data = temp |> filter(species == "Ono.vic"), aes(x = weight, y = Stab, color = species), alpha = 0.5, size = 5) +
    geom_point(data = temp |> mutate(species = ifelse(species == "Ono.vic", "Ono.vic", "Species")) |> 
                 mutate(species = fct_inorder(species)),
               aes(x = weight, y = Stab, color = species, size = species), alpha = 0.5) +
    geom_line(data = fit_slope(q = q)[[2]], aes(x = weight, y = predicted, group = number_of_species), color = "#00BFC4", linewidth = 1) +
    geom_text(data = fit_slope(q = q)[[1]], 
              aes(x = 0.75, y = 0.09, group = number_of_species, label = paste0("Slope = ", sprintf("%.4f", Slope))), color = "#00BFC4",
              size = 4, key_glyph = draw_key_path) +
    geom_text(data = temp_ono, 
              aes(x = 0.61, y = 0.03, group = number_of_species, label = text), color = "#F8766D",
              size = 4, key_glyph = draw_key_path) +
    ylim(c(0, 1.1)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_size_manual(values = c(1, 5)) +
    scale_colour_manual(values = c('black', '#F8766D')) +
    
    # labs(title = paste0("q = ", q), x = "Weight", y = "Stability", color = "Species") +
    # labs(title = paste0("q = ", q), x = "Relative biomass", y = "Stability", color = "Species") +
    labs(title = paste0("q = ", q), x = "Relative biomass", y = "Stability", color = NULL, size = NULL) +
    facet_wrap(~ number_of_species, nrow = 1, labeller = labeller(number_of_species = c(`2` = "2 species", `4` = "4 species", 
                                                                                        `8` = "8 species", `16` = "16 species"))) +
    theme_bw() + 
    theme(legend.position = 'bottom', 
          legend.text = element_text(size = 12),
          # legend.title = element_blank(), 
          legend.margin = margin(0, 0, 0, 0), 
          # legend.box.margin = margin(-10, -10, -5, -10), 
          legend.key.width = unit(1, "cm"),
          text = element_text(size = 14), 
          panel.spacing = unit(0.75, "lines"),
          # plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
          axis.text = element_text(size = 12),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(hjust = 0.5, size = 14),
          axis.title.y = element_text(hjust = 0.5, size = 14)) +
    guides(size = guide_legend(override.aes = list(size = 5)))
  
  
}

gg_alpha_wrt_weight_group_by_species(data = Jena_tidy |> filter(number_of_species != 1), q = 1, threshold = 1 / 6)
gg_alpha_wrt_weight_group_by_species(data = Jena_tidy |> filter(number_of_species != 1), q = 2, threshold = 1 / 6)

