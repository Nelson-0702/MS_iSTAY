library(tidyverse)
library(lmerTest)
library(gridExtra)
library(ggpubr)

load("data/Data_Jena_20_metacommunities.rda")
load("data/Data_Jena_76_community_populations.rda")
load("data/Data_Jena_462_populations.rda")

output_dir <- "figures"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
save_pdf <- function(p, name, width, height) {
  ggplot2::ggsave(file.path(output_dir, name), plot = p, width = width, height = height, bg = "white")
}

# Stability calculation functions-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

#' Calculate information-based temporal stability for a single time series.
#'
#' \code{iSTAY_Single} computes temporal stability measures of order q > 0 for individual time series, accounting for the number of time points.
#'
#' @param data A numeric vector representing a single time series, or a data.frame or matrix with sampling units as rows and time points as columns.
#' @param order.q A numerical vector specifying the orders of stability (q > 0). Default is c(1, 2).
#' @param Alltime Logical indicating whether all time points are used.
#' @param start_T The starting column when Alltime = FALSE.
#' @param end_T The ending column when Alltime = FALSE.
#'
#' @return A data.frame with the following columns: \cr
#' Dataset: identifier of the input time series \cr
#' Order_q: order of stability (q > 0) \cr
#' Stability: information-based temporal stability measure of order q
#' 


iSTAY_Single = function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, end_T = NULL)
{
  if (is.vector(data)) {
    data <- matrix(data, nrow = 1)
  }
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
  stability <- function(vector, q) {
    if (sum(vector != 0) == 0) {
      out <- NA
    }
    else {
      K <- length(vector)
      if (q == 1) {
        H <- sum((vector[vector > 0]/sum(vector)) *
                   log(vector[vector > 0]/sum(vector))) * (-1)
        out <- (exp(H) - 1)/(K - 1)
      }
      else {
        up <- 1 - ((sum((vector[vector > 0]/sum(vector))^q))^(1/(1 - q)))
        out <- up/(1 - K)
      }
    }
    return(out)
  }
  if (Alltime == TRUE) {
    stab <- as.matrix(apply(data, 1, function(vec) sapply(order.q,
                                                          function(qq) stability(vec, q = qq))))
  }
  else {
    subdata <- data[, c(start_T:end_T)]
    stab <- as.matrix(apply(subdata, 1, function(vec) sapply(order.q,
                                                             function(qq) stability(vec, q = qq))))
  }
  result <- data.frame(Assemblage = rep(rownames(as.data.frame(data)),
                                        length(order.q)), Order_q = rep(order.q, each = nrow(data)),
                       Stability = as.vector(t(stab)))
  colnames(result)[1] <- c("Dataset")
  return(result)
}

#' Calculate information-based stability and synchrony for multiple time series.
#'
#' \code{iSTAY_Multiple} computes gamma and alpha stability, multiplicative beta stability, and synchrony of order q > 0 based on the hierarchical decomposition of temporal stability. Both biomass-weighting and equal-weighting schemes are supported.
#'
#' @param data A data.frame or matrix containing multiple time series, with sampling units as rows and time points as columns, or a list of such data frames or matrices.
#' @param order.q A numerical vector specifying the orders of stability and synchrony (q > 0). Default is c(1, 2).
#' @param equal_weights Logical. If TRUE, equal weights are assigned to all sampling units. If FALSE, weights are proportional to their total biomass. Default is TRUE.
#' @param Alltime Logical indicating whether all time points are used.
#' @param start_T The starting column when Alltime = FALSE.
#' @param end_T The ending column when Alltime = FALSE.
#'
#' @return A data.frame with the following columns: \cr
#' Dataset: the input dataset \cr
#' Order_q: order of stability or synchrony \cr
#' Gamma: gamma stability of order q \cr
#' Alpha: alpha stability of order q \cr
#' Beta: multiplicative beta stability (Gamma/Alpha) of order q \cr
#' Synchrony: synchrony measure of order q

iSTAY_Multiple = function (data, order.q = c(1, 2), equal_weights = T, Alltime = TRUE, start_T = NULL, end_T = NULL)
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
  Stabillity_Multiple <- function(ZZ, q){
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    Times <- ncol(ZZ); K <- nrow(ZZ)
    p_i <- as.data.frame(apply(ZZ, 2, function(w) w / z_iplus))
    if(dim(p_i)[2] == 1) p_i = t(p_i)
    if(equal_weights){
      p_pool <-  colMeans(p_i)
      ww <- rep(1/K,K)
    }else{
      z_plusk <- apply(ZZ, 2, sum)
      ww <- z_iplus/sum(ZZ)
      p_pool <- ww %*% as.matrix(p_i)
    }
    if(q == 1){
      alpha <- (exp(sum(ww * apply(p_i, 1, function(p) -sum(p[p > 0] * log(p[p > 0]))))) - 1)/(Times - 1)
      gamma <- (exp(-sum(p_pool[p_pool > 0] * log(p_pool[p_pool > 0]))) - 1)/(Times - 1)
    }else{
      alpha <- (sum(ww * rowSums(p_i^q))^(1/(1-q)) - 1)/(Times - 1)
      gamma <- (sum(p_pool^q)^(1/(1 - q)) - 1)/(Times - 1)
    }
    return( c(gamma, alpha, (gamma/alpha), (gamma - alpha)) )
  }
  Synchrony <- function(ZZ, q){
    z_iplus <- apply(ZZ, 1, sum)
    if (length(which(z_iplus == 0)) != 0) {
      ZZ <- ZZ[-which(z_iplus == 0), ]
    }
    ZZ <- as.matrix(ZZ)
    Times <- ncol(ZZ); K <- nrow(ZZ)
    z_iplus <- apply(ZZ, 1, sum)
    gamma_diversity = Stabillity_Multiple(ZZ,q)[1] * (Times - 1) + 1
    alpha_diversity = Stabillity_Multiple(ZZ,q)[2] * (Times - 1) + 1
    if(sum(rowSums(ZZ)>0)!=1){
      if(dim(ZZ)[1] > dim(ZZ)[2]){
        1 - (1 - alpha_diversity/gamma_diversity)/(1 - 1/gamma_diversity)
      }else{
        p_i <- as.data.frame(apply(ZZ, 2, function(w) w / z_iplus))
        if(equal_weights){
          ww <- rep(1/K, K)
          p_joint <- diag(ww, K) %*% as.matrix(p_i)
        }else{
          ww <- z_iplus/sum(ZZ)
          p_joint <- diag(ww, K) %*% as.matrix(p_i)
        }
        if(q == 1){
          joint_diversity <- exp(-sum(p_joint[p_joint > 0] * log(p_joint[p_joint > 0])) )
        }else{
          joint_diversity <- sum((p_joint)^q)^(1/(1 - q))
        }
        1 - (1 - alpha_diversity/gamma_diversity)/(1 - alpha_diversity/joint_diversity)
      }
    }else{
      1
    }
  }
  if (is.data.frame(data) | is.matrix(data)) {
    if (Alltime == TRUE) {
      subdata <- data
    }else {
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
  }else {
    out <- lapply(order.q, function(qq) {
      cal <- lapply(data, function(ZZ) {
        if (Alltime == T) {
          subZZ <- ZZ
        }else {
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

# Shared plotting functions-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

#' @Function fig_1a is used to plot Figure 1 (a).
#' @param output is an object computed from function 'iSTAY_Multiple' in the package 'iStay'.

fig_1a <- function(output){
  output |>
    ggplot() +
    geom_line(aes(Order_q, Gamma, color = Dataset,linetype = Dataset), linewidth = 2.5)+
    scale_linetype_manual(name="Plot",values = c("solid", "dashed"),labels=c("B1A04_4", "B4A14_2")) +
    scale_color_manual(name="Plot",values = c("red", "blue"),labels=c("B1A04_4", "B4A14_2")) +
    labs(x = "Order.q", y = "Gamma stability") +
    guides(
      linetype = guide_legend(keywidth = 3),
      color = guide_legend(keywidth = 3)
    ) +
    theme_bw() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          text = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(hjust = 0.5, size = 20),
          axis.title.y = element_text(hjust = 0.5, size = 20))+
    theme(panel.border = element_rect(colour="black", size=1.2))
}

#' @Function fig_1b is used to plot Figure 1 (b).
#' @param output is an object computed from function 'iSTAY_Multiple' in the package 'iStay'.
#' @param output2 are objects computed from function '' in the package 'iStay'.

fig_1b <- function(output, output2){
  output |>
    ggplot() +
    geom_line(aes(Order_q, Alpha, color = Dataset, linetype = Dataset), linewidth = 2.5) +
    geom_line(data = output2, aes(Order_q, Stability, color = Site, group = Dataset, linetype = Site), linewidth = 1.5, alpha = 0.3,show.legend = FALSE) +
    scale_color_manual(name="Plot",values = c("red", "blue"),labels=c("B1A04_4", "B4A14_2")) +
    scale_linetype_manual(name = "Plot",values = c("solid","dashed"),labels=c("B1A04_4", "B4A14_2"))+
    labs(x = "Order.q", y = "Alpha stability", color = "Plot") +
    guides(
      linetype = guide_legend(keywidth = 3),
      color = guide_legend(keywidth = 3)
    )+theme_bw() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          text = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(hjust = 0.5, size = 20),
          axis.title.y = element_text(hjust = 0.5, size = 20))+
    theme(panel.border = element_rect(colour="black", size=1.2))
}

#' @Function fig_1c is used to plot Figure 1 (c).
#' @param output is an object computed from function 'iSTAY_Multiple' in the package 'iStay'.

fig_1c <- function(output){
  output |>
    ggplot() +
    geom_line(aes(Order_q, Synchrony, color = Dataset,linetype = Dataset), linewidth = 2.5) +
    labs(x = "Order.q", y = "Synchrony among species", color = "Plot") +
    scale_linetype_manual(name="Plot",values = c("solid", "dashed"),labels=c("B1A04_4", "B4A14_2")) +
    scale_color_manual(name="Plot",values = c("red", "blue"),labels=c("B1A04_4", "B4A14_2")) +
    guides(
      linetype = guide_legend(keywidth = 3),
      color = guide_legend(keywidth = 3)
    )+
    theme_bw() +
    theme(legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
          axis.text = element_text(size = 16),
          text = element_text(size = 14),
          strip.text.x = element_text(size = 14),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(hjust = 0.5, size = 20),
          axis.title.y = element_text(hjust = 0.5, size = 20))+
    theme(panel.border = element_rect(colour="black", size=1.2))
}

#' @Function LMM_2_to_4 is used to calculate linear mixed model results for Figure 2 to Figure4.
#' @param output is an object computed from function 'iSTAY_Multiple' in the package 'iStay'.
#' @param structure A data.frame containing block-level structural information, 
#' typically including columns such as `block` and `log2_sowndiv`. This is used to fit 
#' linear mixed-effects models and to associate each diversity level with its corresponding block.
#' @param metric_name A character string specifying the ecological stability metric 
#' to be analyzed. Acceptable values include column names from the `output` object, 
#' such as `"Gamma"`, `"Alpha"`, or `"Synchrony"`.

LMM_2_to_4 <- function(output, structure, metric_name) {
  plotdata <- data.frame(
    block = rep(structure$block, 3),
    log2_sowndiv = rep(structure$log2_sowndiv, 3),
    value = c(output[[metric_name]]),
    version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(structure)))
  )
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  for (v in levels(plotdata$version)) {
    subdata <- dplyr::filter(plotdata, version == v)
    model <- lmerTest::lmer(value ~ 1 + log2_sowndiv + (1 + log2_sowndiv | block), data = subdata)
    summary_model <- summary(model)
    subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
    subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
    slopes <- coef(model)$block
    slope_intercept <- data.frame(
      block = rownames(slopes),
      intercept = slopes[, 1],
      slope = slopes[, 2],
      x_min = tapply(subdata$log2_sowndiv, subdata$block, min),
      x_max = tapply(subdata$log2_sowndiv, subdata$block, max),
      version = v
    )
    slope_text_metric <- data.frame(
      block = c(rownames(slopes), "Total"),
      slope_text = c(
        paste0("slope = ", sprintf("%.4f", slopes[, 2])),
        paste0("slope = ", sprintf("%.4f", summary_model$coefficients[2, 1]))
      ),
      version = v
    )
    combined_total_fit <- rbind(combined_total_fit, subdata)
    combined_part_fit <- rbind(combined_part_fit, slope_intercept)
    combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
  }
  names(combined_total_fit)[names(combined_total_fit) == "value"] <- metric_name
  names(plotdata)[names(plotdata) == "value"] <- metric_name
  return(list(
    total_fit = combined_total_fit,
    part_fit = combined_part_fit,
    slope_text = combined_slope_text,
    plotdata = plotdata
  ))
}

#' @Function fig_2_or_4 is used to plot Figure 2 or Figure4.
#' @param output is an object computed from function LMM_2_or_4.
#' @param metric_name A character string specifying the ecological stability metric 
#' to be analyzed. Acceptable values include column names from the `output` object, 
#' such as `"Gamma"`, `"Alpha"`, or `"Synchrony"`.

fig2_or_4 <- function(output, metric_name) {
  color_palette <- c("#EA0000", "#4f772d", "#0066CC", "#fb8500")
  total_fit <- output$total_fit
  part_fit <- output$part_fit
  slope_text <- output$slope_text
  total_fit$significance <- factor(total_fit$significance, levels = c("significant", "non-significant"))
  if (!any(total_fit$significance == "non-significant")) {
    dummy <- total_fit[1, ]
    dummy$log2_sowndiv <- NA
    dummy$predicted <- NA
    metric_col <- setdiff(names(total_fit), c("log2_sowndiv", "predicted", "significance", "sign", "block", "version"))
    dummy[[metric_name]] <- NA
    dummy$significance <- "non-significant"
    total_fit <- rbind(total_fit, dummy)
  }
  p <- ggplot()
  if (identical(output, output_fig_2a)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Gamma, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 0.5, y = 0.25, label = slope_text),
                hjust = c(-1.2,-1.2,-1.2), vjust = rep(0,3), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 0.5, y = 0.25, label = slope_text, color = block,
                    hjust = rep(-c(1.2,0,1.2,0), 3),
                    vjust = rep(c(2,2,4,4), 3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Community stability")  +
      facet_wrap(~ version)+
      scale_y_continuous(limits = c(0.15, 1.0))
  } else if (identical(output, output_fig_2b)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Alpha, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 0.5, y = 0.7, label = slope_text),
                hjust = rep(-1.25, 3), vjust = c(-8,-8,-8), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 0.5, y = 0.7, label = slope_text, color = block,
                    hjust = rep(-c(0,1.25,0,1.25), 3),
                    vjust = rep(c(-6, -6, -4, -4),3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Population stability") +
      facet_wrap(~ version)+
      scale_y_continuous(limits = c(0, 1.0))
  } else if (identical(output, output_fig_2c)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Synchrony, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 1.25, y = 0.40, label = slope_text),
                hjust = rep(-1.25, 3), vjust = c(-16,-16,-16), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 1.25, y = 0.40, label = slope_text, color = block,
                    hjust = rep(-c(0, 1.25, 0, 1.25), 3),
                    vjust = rep(c(-14, -14, -12, -12),3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Population synchrony") +
      facet_wrap(~ version) +
      coord_cartesian(ylim = c(0, 1.0))
  } else if (identical(output, output_fig_4a)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Gamma, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 0.5, y = 0.6, label = slope_text),
                hjust = rep(-1.25, 3), vjust = c(7,7,7), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 0.5, y = 0.6, label = slope_text, color = block,
                    hjust = rep(-c(0,1.25,0,1.25), 3),
                    vjust = rep(c(9, 9, 11, 11),3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Gamma stability") +
      facet_wrap(~ version) +
      scale_y_continuous(limits = c(0.35, 1.05))
  } else if (identical(output, output_fig_4b)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Alpha, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 0.5, y = 0.6, label = slope_text),
                hjust = rep(-1.25, 3), vjust = c(8,8,8), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 0.5, y = 0.6, label = slope_text, color = block,
                    hjust = rep(-c(0,1.25,0,1.25), 3),
                    vjust = rep(c(10, 10, 12, 12),3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Alpha stability") +
      facet_wrap(~ version) +
      scale_y_continuous(limits = c(0.3, 1.0))
  } else if (identical(output, output_fig_4c)) {
    p <- p +
      geom_point(data = total_fit, aes(x = log2_sowndiv, y = Synchrony, color = block), size = 3, alpha = 0.5) +
      geom_text(data = slope_text |> filter(block == "Total"),
                aes(x = 0.5, y = 0.8, label = slope_text),
                hjust = rep(-1.25, 3), vjust = rep(6, 3), size=4.5) +
      geom_text(data = slope_text |> filter(block != "Total"),
                aes(x = 0.5, y = 0.8, label = slope_text, color = block,
                    hjust = rep(-c(0, 1.25, 0, 1.25), 3),
                    vjust = rep(c(8, 8, 10, 10),3)),
                size=4.5, key_glyph = draw_key_path) +
      labs(y = "Plot spatial synchrony") +
      facet_wrap(~ version) +
      coord_cartesian(ylim = c(0.7, 1.0))
  }
  p <- p +
    geom_segment(data = part_fit,
                 aes(x = x_min, xend = x_max,
                     y = intercept + slope * x_min,
                     yend = intercept + slope * x_max,
                     color = block)) +
    geom_line(data = total_fit,
              aes(x = log2_sowndiv, y = predicted, linetype = significance), linewidth = 1.5) +
    scale_color_manual(values = color_palette) +
    scale_linetype_manual(values = c("significant" = "solid", "non-significant" = "dashed"),
                          labels = c("significant" = "Significant", "non-significant" = "Non-significant"),
                          drop = FALSE)+
    scale_x_continuous(breaks = c(0, 1, 2, 3, 4), labels = c(1, 2, 4, 8, 16)) +
    labs(linetype = "", x = "Number of species (log2 scale)") +
    guides(linetype = guide_legend(keywidth = 3.1)) +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 13),
          legend.box.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
          legend.text = element_text(size = 14),
          legend.title = element_blank(),
          axis.text = element_text(size = 16),
          text = element_text(size = 14),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 14),
          axis.title.x = element_text(hjust = 0.5, size = 20),
          axis.title.y = element_text(hjust = 0.5, size = 20))+
    theme(panel.border = element_rect(colour="black", size=1.2))
  return(p)
}

#' @Function fig_3_left is used to plot left panel of Figure 3.
#' @param summart_df is a summary result for the output of function 'iSTAY_Multiple' in the package 'iStay'.

fig_3_left <- function(summary_df) {
  color_palette <- c("#F8766D", "#A3A500", "#00BF7D", "#343a40", "#E76BF3")
  shape_palette <- c(20, 15, 18, 17, 19)
  x_breaks <- c(2003, 2005, 2007, 2009, 2011, 2013, 2015)
  base_theme <- theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 13),
      legend.box.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.text = element_text(size = 14),
      text = element_text(size = 14),
      axis.text = element_text(size = 14),
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(hjust = 0.5, size = 20),
      axis.title.y = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )  +
    theme(panel.border = element_rect(colour="black", size=1.2))
  p1 <- ggplot(summary_df) +
    geom_point(aes(x = Start_year, y = mean_gamma, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_gamma, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),
      color    = guide_legend(keywidth = 3)
    ) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    scale_x_continuous(breaks = x_breaks) +
    labs(y = "Community stability", x = "Starting year of sliding 10-year windows", color = "Number of species", shape = "Number of species") +
    base_theme
  p2 <- ggplot(summary_df) +
    geom_point(aes(x = Start_year, y = mean_alpha, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_alpha, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),
      color    = guide_legend(keywidth = 3)
    ) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    scale_x_continuous(breaks = x_breaks) +
    labs(y = "Population stability", x = "Starting year of sliding 10-year windows", color = "Number of species", shape = "Number of species") +
    base_theme
  p3 <- summary_df |>
    filter(log2_sowndiv != 0) |>
    ggplot() +
    geom_point(aes(x = Start_year, y = mean_syn, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_syn, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),
      color    = guide_legend(keywidth = 3)
    ) +
    scale_color_manual(values = color_palette[-1]) +
    scale_shape_manual(values = shape_palette[-1]) +
    scale_x_continuous(breaks = x_breaks) +
    labs(y = "Population synchrony", x = "Starting year of sliding 10-year windows", color = "Number of species", shape = "Number of species") +
    base_theme
  return(list(Gamma_plot = p1, Alpha_plot = p2, Synchrony_plot = p3))
}

#' @Function slope_3 is used to calculate slope over window for Figure 3.
#' @param metric_name A character string specifying the ecological stability metric 
#' to be analyzed. Acceptable values include column names from the `output` object, 
#' such as `"Gamma"`, `"Alpha"`, or `"Synchrony"`.

slope_3 <- function(metric_name) {
  lapply(seq_along(year_windows), function(i) {
    data_10yr <- lapply(Data_Jena_76_community_populations, \(df) df[, year_windows[[i]], drop = FALSE])
    names(data_10yr) <- names(Data_Jena_76_community_populations)
    output <- iSTAY_Multiple(data_10yr, order.q = c(0.5, 1, 2), equal_weights = T)
    structure_window <- structure3
    if(metric_name == "Synchrony") {
      keep <- which(structure3$log2_sowndiv != 0)
      output <- output[unlist(lapply(0:2, function(j) keep + j * nrow(structure3))), ]
      structure_window <- structure3[keep, , drop = FALSE]
    }
    result <- LMM_2_to_4(output, structure_window, metric_name = metric_name)
    result$slope_text |>
      dplyr::filter(block == "Total") |>
      dplyr::mutate(
        starting_year = as.numeric(year_windows[[i]][1]),
        slope_num = readr::parse_number(slope_text)
      )
  }) |> dplyr::bind_rows()
}

#' @Function fig_3_right is used to plot right panel of Figure 3.
#' @param output is an object computed from function slpoe_3.

fig_3_right <- function(output) {
  y_label <- "Effect of species richness"
  default_colors <- c("#386641", "#6a994e", "#a7c957")
  ggplot(output) +
    geom_line(aes(x = starting_year, y = slope_num, color = version, linetype = version), linewidth = 1) +
    geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
    scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
    scale_color_manual(values = default_colors) +
    labs(y = y_label, x = "Starting year of sliding 10-year windows") +
    guides(
      linetype = guide_legend(keywidth = 3),
      color    = guide_legend(keywidth = 3)
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 13),
      legend.box.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm"),
      legend.text = element_text(size = 18),
      legend.title = element_blank(),
      text = element_text(size = 14),
      axis.text = element_text(size = 16),
      strip.text.x = element_text(size = 20),
      strip.text.y = element_text(size = 14),
      axis.title.x = element_text(hjust = 0.5, size = 20),
      axis.title.y = element_text(hjust = 0.5, size = 20)
    )+
    theme(panel.border = element_rect(colour="black", size=1.2))
}

# ============================================================================================================================================================================================================================================================================================================================== #
# 
# This code includes four parts:
# 
# (1) Figure 1. Gamma, alpha and synchrony profiles for q between 0 and 2 within each plot. 
# (2) Figure 2. Biodiversity–stability and biodiversity–synchrony relationships based on 76 plots.
# (3) Figure 3. Temporal effects of species richness on stability and synchrony based on 12 consecutive overlapping 10-year moving window.
# (4) Figure 4. Biodiversity–stability and biodiversity–synchrony relationships based on 20 sets
#
# See "Brief guide" for details. 
# 
# ============================================================================================================================================================================================================================================================================================================================== #


# ============================================================================================================================================================================================================================================================================================================================== #
# Figure 1. Gamma, alpha and synchrony profiles for q between 0 and 2 within each plot.-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

df1 <- list(B4A14_2 = Data_Jena_76_community_populations$B4A14_B4_2,
            B1A04_4 = Data_Jena_76_community_populations$B1A04_B1_4)
df1 <- Map(function(x, nm) {
  rownames(x) <- paste0(nm, rownames(x))
  x
}, df1, names(df1))
output_fig_1 <- iSTAY_Multiple(df1, order.q = seq(0.1, 2, 0.01), equal_weights = F)
output_fig_1b <- list(
  B1A04_4 = iSTAY_Single(df1$B1A04_4, order.q = seq(0.1, 2, 0.01)),
  B4A14_2 = iSTAY_Single(df1$B4A14_2, order.q = seq(0.1, 2, 0.01))
) |>
  bind_rows(.id = "Site")

## Figure 1 (a)
fig_1a = fig_1a(output_fig_1)

## Figure 1 (b)
fig_1b = fig_1b(output_fig_1, output_fig_1b)

## Figure 1 (c)
fig_1c = fig_1c(output_fig_1)

p1a <- annotate_figure( fig_1a, top = text_grob("(a) Plot or gamma stability", size = 18, face = "bold", hjust = 1.05) )
p1b <- annotate_figure( fig_1b, top = text_grob("(b) Species or alpha stability", size = 18, face = "bold", hjust = 0.92) )
row1 <- arrangeGrob( p1a, p1b, ncol = 2, widths = c(8, 8) )
p1c <- annotate_figure( fig_1c, top = text_grob("(c) Synchrony among species", size = 18, face = "bold", hjust = 0.92) )
blank <- grid::nullGrob()
row2 <- arrangeGrob( p1c, blank, ncol = 2, widths = c(8, 8) )
g <- arrangeGrob(row1, row2, grid::nullGrob(), ncol = 1, heights = c(1, 1, 1))
save_pdf(g, "Figure_1.pdf", 16, 20)

# ============================================================================================================================================================================================================================================================================================================================== #
# Figure 2. Biodiversity–stability and biodiversity–synchrony relationships based on 76 plots.-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

split_names2 <- str_split(names(Data_Jena_76_community_populations), "_", simplify = TRUE)
structure2 <- data.frame(
  block = split_names2[, 2],
  log2_sowndiv = log2(as.numeric(split_names2[, 3]))
)
keep2 <- which(structure2$log2_sowndiv != 0)
structure2c <- structure2[keep2, , drop = FALSE]
output2 <- iSTAY_Multiple(Data_Jena_76_community_populations, order.q = c(0.5, 1, 2), equal_weights = F)
output_fig_2a <- LMM_2_to_4(output2, structure = structure2, metric_name = "Gamma")

## Figure 2 (a)
fig_2a = fig2_or_4(output_fig_2a, metric_name = "Gamma")
output_fig_2b <- LMM_2_to_4(output2, structure = structure2, metric_name = "Alpha")

## Figure 2 (b)
fig_2b = fig2_or_4(output_fig_2b, metric_name = "Alpha")
output2c <- output2[unlist(lapply(0:2, function(j) keep2 + j * nrow(structure2))), ]
output_fig_2c <- LMM_2_to_4(output2c, structure = structure2c, metric_name = "Synchrony")

## Figure 2 (c)
fig_2c = fig2_or_4(output_fig_2c, metric_name = "Synchrony")

fig2_all <- arrangeGrob(
  annotate_figure(fig_2a, top = text_grob("(a) Biodiversity–stability relationships for plot/gamma stability", size = 20, face = "bold", just = "left", hjust = 0.62)),
  annotate_figure(fig_2b, top = text_grob("(b) Biodiversity–stability relationships for species/alpha stability", size = 20, face = "bold", just = "left", hjust = 0.60)),
  annotate_figure(fig_2c, top = text_grob("(c) Biodiversity–synchrony relationships for species synchrony within a plot", size = 20, face = "bold", just = "left", hjust = 0.51)),
  ncol = 1)
save_pdf(fig2_all, "Figure_2.pdf", 12, 20)

# ============================================================================================================================================================================================================================================================================================================================== #
# Figure 3. Temporal effects of species richness on stability and synchrony based on 12 consecutive overlapping 10-year moving window.-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

### Left

split_names3 <- str_split(names(Data_Jena_76_community_populations), "_", simplify = TRUE)

# Create 10-year moving windows, excluding those containing 2004

year_windows <- lapply(2003:2015, function(start) {
  yrs <- if (start == 2003) c(2003, 2005:2013) else start:(start + 9)
  if (2004 %in% yrs) return(NULL)
  as.character(yrs)
}) |> compact()
names(year_windows) <- as.character(c(2003, 2005:2015))
output_fig_3_left <- lapply(names(year_windows), function(start_year) {
  window_years <- year_windows[[start_year]]
  biomass_data <- lapply(Data_Jena_76_community_populations, \(df) df[, window_years, drop = FALSE])
  output3_left <- iSTAY_Multiple(biomass_data, order.q = c(0.5, 1, 2), equal_weights = T)
  output3_left |>
    mutate(
      Start_year = as.numeric(start_year),
      log2_sowndiv = rep(log2(as.numeric(split_names3[, 3])), 3)
    )
})
Summary_fig_3_left <- bind_rows(output_fig_3_left) |>
  mutate(Order_q = paste0("q = ", Order_q)) |>
  group_by(Start_year, Order_q, log2_sowndiv) |>
  summarise(
    mean_gamma = mean(Gamma),
    mean_alpha = mean(Alpha),
    mean_syn = mean(Synchrony),
    .groups = "drop"
  )

## Figure 3 (a) left
fig_3a_left = fig_3_left(Summary_fig_3_left)$Gamma_plot

## Figure 3 (b) left
fig_3b_left = fig_3_left(Summary_fig_3_left)$Alpha_plot

## Figure 3 (c) left
fig_3c_left = fig_3_left(Summary_fig_3_left)$Synchrony_plot
structure3 <- data.frame(
  block = split_names3[, 2],
  log2_sowndiv = log2(as.numeric(split_names3[, 3]))
)

## Figure 3 (a) right
output_fig_3a_right <- slope_3(metric_name = "Gamma")
fig_3a_right = fig_3_right(output_fig_3a_right)

## Figure 3 (b) right
output_fig_3b_right <- slope_3(metric_name = "Alpha")
fig_3b_right = fig_3_right(output_fig_3b_right)

## Figure 3 (c) right
output_fig_3c_right <- slope_3(metric_name = "Synchrony")
fig_3c_right = fig_3_right(output_fig_3c_right)

p_a <- annotate_figure(fig_3a_left,
                       top = text_grob("(a) Gamma stability", size = 18, face = "bold",
                                       hjust = 2.3))
p_b <- annotate_figure(fig_3a_right,
                       top = text_grob("(b) Effect of biodiversity on gamma stability", size = 18, face = "bold",
                                       hjust = 0.61))
row1 <- arrangeGrob(p_a, p_b, ncol = 2, widths = c(12, 8))
p_c <- annotate_figure(fig_3b_left,
                       top = text_grob("(c) Alpha stability", size = 18, face = "bold", hjust = 2.55))
p_d <- annotate_figure(fig_3b_right,
                       top = text_grob("(d) Effect of biodiversity on alpha stability", size = 18, face = "bold", hjust = 0.61))
row2 <- arrangeGrob(p_c, p_d, ncol = 2, widths = c(12, 8))
p_e <- annotate_figure(fig_3c_left,
                       top = text_grob("(e) Species-level synchrony", size = 18, face = "bold", hjust = 1.59))
p_f <- annotate_figure(fig_3c_right,
                       top = text_grob("(f) Effect of biodiversity on synchrony", size = 18, face = "bold", hjust = 0.68))
row3 <- arrangeGrob(p_e, p_f, ncol = 2, widths = c(12, 8))
g <- arrangeGrob(row1, row2, row3, ncol = 1, heights = c(1,1,1))
save_pdf(g, "Figure_3.pdf", 20, 18)

# ============================================================================================================================================================================================================================================================================================================================== #
# Figure 4. Biodiversity–stability and biodiversity–synchrony relationships based on 20 sets-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --

split_names4 <- str_split(names(Data_Jena_20_metacommunities), "_", simplify = TRUE)
structure4 <- data.frame(
  block = split_names4[, 1],
  log2_sowndiv = log2(as.numeric(split_names4[, 2]))
)
output4 <- iSTAY_Multiple(Data_Jena_20_metacommunities, order.q = c(0.5, 1, 2), equal_weights = F)

## Figure 4 (a)
output_fig_4a <- LMM_2_to_4(output4, structure = structure4, metric_name = "Gamma")
fig_4a = fig2_or_4(output_fig_4a, metric_name = "Gamma")

## Figure 4 (b)
output_fig_4b <- LMM_2_to_4(output4, structure = structure4, metric_name = "Alpha")
fig_4b = fig2_or_4(output_fig_4b, metric_name = "Alpha")


## Figure 4 (c)
output_fig_4c <- LMM_2_to_4(output4, structure = structure4, metric_name = "Synchrony")
fig_4c = fig2_or_4(output_fig_4c, metric_name = "Synchrony")
fig4_all <- arrangeGrob(
  annotate_figure(fig_4a, top = text_grob("(a) Biodiversity–stability relationships for gamma stability", size = 20, face = "bold", just = "left", hjust = 0.68)),
  annotate_figure(fig_4b, top = text_grob("(b) Biodiversity–stability relationships for alpha stability", size = 20, face = "bold", just = "left", hjust = 0.70)),
  annotate_figure(fig_4c, top = text_grob("(c) Biodiversity-synchrony relationships for spatial synchrony among plots", size = 20, face = "bold", just = "left", hjust = 0.525)),
  ncol = 1)
save_pdf(fig4_all, "Figure_4.pdf", 12, 20)
