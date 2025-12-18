# ========================================================================================================== #

# Install and upload libraries ----
library(ggh4x)
library(tidyverse)
library(lmerTest)
library(ggplot2)

#install.packages(iSTAY)
library(iSTAY) # Just for the Jena dataset. 

data("Data_Jena_20_metacommunities")
data("Data_Jena_76_metapopulations")

#' Calculate stability for a single time series.
#'
#' \code{iSTAY_Single} computes the stability of order q for a single time series.
#'
#' @param data A \code{vector} of time series data, or a \code{data.frame} with sampling units as rows and time points as columns.
#' @param order.q a numerical vector specifying the orders of stability. Default is c(1,2).
#' @param Alltime Logical (\code{TRUE} or \code{FALSE}), indicating whether to use all time points in the data.
#' @param start_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the starting column (time point) for the analysis interval.
#' @param end_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the ending column (time point) for the analysis interval.
#'
#'
#'
#' @return a dataframe with columns: \cr
#' Dataset: the input dataset \cr
#' Order_q: order of stability \cr
#' Stability: stability measures of order q

iSTAY_Single = function (data, order.q = c(1, 2), Alltime = TRUE, start_T = NULL, 
                         end_T = NULL) 
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


#' Calculate stability and synchrony for multiple time series.
#'
#' \code{iSTAY_Multiple} computes gamma, alpha, and beta stability, as well as synchrony, for multiple time-series data.
#'
#' @param data A \code{data.frame} containing multiple time series data, with sampling units as rows and time points as columns, or a \code{list} of \code{data.frames} with each data frame representing multiple time series.
#' @param order.q A numerical vector specifying the orders of stability and synchrony. Default is c(1,2).
#' @param Alltime Logical (\code{TRUE} or \code{FALSE}), indicating whether to use all time points in the data.
#' @param start_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the starting column (time point) for the analysis interval.
#' @param end_T (Applicable only if \code{Alltime = FALSE}) a positive integer specifying the ending column (time point) for the analysis interval.
#'
#'
#' @return a data frame with the following columns: \cr
#'  Dataset: the input dataset \cr
#'  Order_q: order of stability or synchrony \cr
#'  Gamma, Alpha, Beta: stability measures of order q \cr 
#'  Synchrony: synchrony measure of order q

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



#' @Function fig_1a is used to plot Figure 1 (a).
#' @param output is an object computed from function 'iSTAY_Multiple' in the package 'iStay'.
#' 
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
#' 
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
#' 
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
#'

LMM_2_to_4 <- function(output, structure, metric_name) {
  # Build plotdata using the selected metric column from precomputed output
  plotdata <- data.frame(
    block = rep(structure$block, 3),
    log2_sowndiv = rep(structure$log2_sowndiv, 3),
    value = c(output[[metric_name]]),
    version = factor(rep(c("q = 0.5", "q = 1", "q = 2"), each = nrow(structure)))
  )
  
  # Initialize containers for results
  combined_total_fit <- data.frame()
  combined_part_fit <- data.frame()
  combined_slope_text <- data.frame()
  
  # Loop over each q version
  for (v in levels(plotdata$version)) {
    subdata <- dplyr::filter(plotdata, version == v)
    
    # Fit linear mixed-effects model
    model <- lmerTest::lmer(value ~ 1 + log2_sowndiv + (1 + log2_sowndiv | block), data = subdata)
    summary_model <- summary(model)
    
    # Add predictions and significance markers
    subdata$predicted <- predict(model, newdata = subdata, re.form = NA)
    subdata$significance <- ifelse(summary_model$coefficients[2, 5] < 0.05, "significant", "non-significant")
    
    # Get slopes and intercepts per block
    slopes <- coef(model)$block
    slope_intercept <- data.frame(
      block = rownames(slopes),
      intercept = slopes[, 1],
      slope = slopes[, 2],
      x_min = tapply(subdata$log2_sowndiv, subdata$block, min),
      x_max = tapply(subdata$log2_sowndiv, subdata$block, max),
      version = v
    )
    
    # Create slope text for plot annotations
    slope_text_metric <- data.frame(
      block = c(rownames(slopes), "Total"),
      slope_text = c(
        paste0("slope = ", sprintf("%.4f", slopes[, 2])),
        paste0("slope = ", sprintf("%.4f", summary_model$coefficients[2, 1]))
      ),
      version = v
    )
    
    # Combine results
    combined_total_fit <- rbind(combined_total_fit, subdata)
    combined_part_fit <- rbind(combined_part_fit, slope_intercept)
    combined_slope_text <- rbind(combined_slope_text, slope_text_metric)
  }
  
  # Rename "value" column to the actual metric name
  names(combined_total_fit)[names(combined_total_fit) == "value"] <- metric_name
  names(plotdata)[names(plotdata) == "value"] <- metric_name
  
  # Return results
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
#'

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
      scale_y_continuous(limits = c(0, 1.0))
    
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
      scale_y_continuous(limits = c(0.35, 1.0))
    
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
      scale_y_continuous(limits = c(0.7, 1.0))
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
#'

fig_3_left <- function(summary_df) {
  # Define shared settings
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
  
  # Plot for Gamma
  p1 <- ggplot(summary_df) +
    geom_point(aes(x = Start_year, y = mean_gamma, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_gamma, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),   # 讓線型 legend 變長
      color    = guide_legend(keywidth = 3)    # 讓顏色 legend 也變長
    ) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    scale_x_continuous(breaks = x_breaks) +
    labs(y = "Community stability", x = "Starting year of sliding 10-year windows", color = "Number of species", shape = "Number of species") +
    base_theme
  
  # Plot for Alpha
  p2 <- ggplot(summary_df) +
    geom_point(aes(x = Start_year, y = mean_alpha, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_alpha, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),   # 讓線型 legend 變長
      color    = guide_legend(keywidth = 3)    # 讓顏色 legend 也變長
    ) +
    scale_color_manual(values = color_palette) +
    scale_shape_manual(values = shape_palette) +
    scale_x_continuous(breaks = x_breaks) +
    labs(y = "Population stability", x = "Starting year of sliding 10-year windows", color = "Number of species", shape = "Number of species") +
    base_theme
  
  # Plot for Synchrony (excluding log2_sowndiv == 1)
  p3 <- summary_df |> 
    filter(log2_sowndiv != 1) |> 
    ggplot() +
    geom_point(aes(x = Start_year, y = mean_syn, color = factor(log2_sowndiv), shape = factor(log2_sowndiv)), size = 3) +
    geom_line(aes(x = Start_year, y = mean_syn, color = factor(log2_sowndiv)), linewidth = 1) +
    facet_wrap(~ Order_q, nrow = 1) +
    guides(
      linetype = guide_legend(keywidth = 3),   # 讓線型 legend 變長
      color    = guide_legend(keywidth = 3)    # 讓顏色 legend 也變長
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
#'

slope_3 <- function(metric_name) {
  lapply(seq_along(year_windows), function(i) {
    
    # Extract 10-year data window
    data_10yr <- lapply(Data_Jena_76_metapopulations, \(df) df[, year_windows[[i]], drop = FALSE])
    names(data_10yr) <- names(Data_Jena_76_metapopulations)
    
    # Compute metric output and fit LMM
    output <- iSTAY_Multiple(data_10yr, order.q = c(0.5, 1, 2))
    if(metric_name == "Synchrony"){ 
      output <- output |> filter(Synchrony != 1)
      structure3 <- structure3 |> filter(log2_sowndiv != 0)
    }
    result <- LMM_2_to_4(output, structure3, metric_name = metric_name)
    
    # Extract slope information for Total block
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
#'
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
      linetype = guide_legend(keywidth = 3),   # 讓線型 legend 變長
      color    = guide_legend(keywidth = 3)    # 讓顏色 legend 也變長
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




# ========================================================================================================== #
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
# ========================================================================================================== #


# library(devtools)
# install_github("AnneChao/iSTAY")    # Press 'Enter' to skip update options

# source("Source R code.txt")


# ========================================================================================================== #
# Figure 1. Gamma, alpha and synchrony profiles for q between 0 and 2 within each plot.
df1 <- list(B4A14_2 = Data_Jena_76_metapopulations$B4A14_B4_2,
            B1A04_4 = Data_Jena_76_metapopulations$B1A04_B1_4)
df1 <- Map(function(x, nm) {
  rownames(x) <- paste0(nm, rownames(x))
  x
}, df1, names(df1))

# df1$B4A14_2[df1$B4A14_2 == 0] = 10^-5
# df1$B1A04_4[df1$B1A04_4 == 0] = 10^-5

output_fig_1 <- iSTAY_Multiple(df1, order.q = seq(0.1, 2, 0.01))
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


ggsave("Figure_1a.png", fig_1a, width = 8, height = 6, dpi = 1000)
ggsave("Figure_1b.png", fig_1b, width = 8, height = 6, dpi = 1000)
ggsave("Figure_1c.png", fig_1c, width = 8, height = 6, dpi = 1000)

# ========================================================================================================== #
# Figure 2. Biodiversity–stability and biodiversity–synchrony relationships based on 76 plots.

split_names2 <- str_split(names(Data_Jena_76_metapopulations), "_", simplify = TRUE)

structure2 <- data.frame(
  block = split_names2[, 2],
  log2_sowndiv = log2(as.numeric(split_names2[, 3]))
)

structure2c <- structure2 |> filter(log2_sowndiv != 0)

output2 <- iSTAY_Multiple(Data_Jena_76_metapopulations, order.q = c(0.5, 1, 2))

## Figure 2 (a)

output_fig_2a <- LMM_2_to_4(output2, structure = structure2, metric_name = "Gamma")

fig_2a = fig2_or_4(output_fig_2a, metric_name = "Gamma")

## Figure 2 (b)

output_fig_2b <- LMM_2_to_4(output2, structure = structure2, metric_name = "Alpha")

fig_2b = fig2_or_4(output_fig_2b, metric_name = "Alpha")

## Figure 2 (c)

output_fig_2c <- LMM_2_to_4(output2 |> filter(Synchrony != 1), structure = structure2c, metric_name = "Synchrony")

fig_2c = fig2_or_4(output_fig_2c, metric_name = "Synchrony")


ggsave("Figure_2a.png", fig_2a, width = 12, height = 6, dpi = 1000)
ggsave("Figure_2b.png", fig_2b, width = 12, height = 6, dpi = 1000)
ggsave("Figure_2c.png", fig_2c, width = 12, height = 6, dpi = 1000)

# ========================================================================================================== #
# Figure 3. Temporal effects of species richness on stability and synchrony based on 12 consecutive overlapping 10-year moving window.

### Left

split_names3 <- str_split(names(Data_Jena_76_metapopulations), "_", simplify = TRUE)

# Create 10-year moving windows, excluding those containing 2004
year_windows <- lapply(2003:2015, function(start) {
  yrs <- if (start == 2003) c(2003, 2005:2013) else start:(start + 9)
  if (2004 %in% yrs) return(NULL)
  as.character(yrs)
}) |> compact()

names(year_windows) <- as.character(c(2003, 2005:2015))

output_fig_3_left <- lapply(names(year_windows), function(start_year) {
  window_years <- year_windows[[start_year]]

  biomass_data <- lapply(Data_Jena_76_metapopulations, \(df) df[, window_years, drop = FALSE])

  output3_left <- iSTAY_Multiple(biomass_data, order.q = c(0.5, 1, 2))

  output3_left |>
    mutate(
      Start_year = as.numeric(start_year),
      log2_sowndiv = rep(as.numeric(split_names3[, 3]), 3)
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


### Right

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
#   geom_hline(yintercept = 0, color = "dodgerblue1", linewidth = 2)

## Figure 3 (c) right

output_fig_3c_right <- slope_3(metric_name = "Synchrony")

fig_3c_right = fig_3_right(output_fig_3c_right)


# Left
ggsave("Figure_3a_left.png", fig_3a_left, width = 12, height = 6, dpi = 1000)
ggsave("Figure_3b_left.png", fig_3b_left, width = 12, height = 6, dpi = 1000)
ggsave("Figure_3c_left.png", fig_3c_left, width = 12, height = 6, dpi = 1000)

# Right
ggsave("Figure_3a_right.png", fig_3a_right, width = 8, height = 6, dpi = 1000)
ggsave("Figure_3b_right.png", fig_3b_right, width = 8, height = 6, dpi = 1000)
ggsave("Figure_3c_right.png", fig_3c_right, width = 8, height = 6, dpi = 1000)


# ========================================================================================================== #
# Figure 4. Biodiversity–stability and biodiversity–synchrony relationships based on 20 sets

split_names4 <- str_split(names(Data_Jena_20_metacommunities), "_", simplify = TRUE)

structure4 <- data.frame(
  block = split_names4[, 1],
  log2_sowndiv = log2(as.numeric(split_names4[, 2]))
)

output4 <- iSTAY_Multiple(Data_Jena_20_metacommunities, order.q = c(0.5, 1, 2))

## Figure 4 (a)

output_fig_4a <- LMM_2_to_4(output4, structure = structure4, metric_name = "Gamma")

fig_4a = fig2_or_4(output_fig_4a, metric_name = "Gamma")

## Figure 4 (b)

output_fig_4b <- LMM_2_to_4(output4, structure = structure4, metric_name = "Alpha")

fig_4b = fig2_or_4(output_fig_4b, metric_name = "Alpha")

## Figure 4 (c)

output_fig_4c <- LMM_2_to_4(output4, structure = structure4, metric_name = "Synchrony")

fig_4c = fig2_or_4(output_fig_4c, metric_name = "Synchrony")


ggsave("Figure_4a.png", fig_4a, width = 12, height = 6, dpi = 1000)
ggsave("Figure_4b.png", fig_4b, width = 12, height = 6, dpi = 1000)
ggsave("Figure_4c.png", fig_4c, width = 12, height = 6, dpi = 1000)


# ========================================================================================================== #
# Figure S4 Relationships between biodiversity and gamma stability, alpha stability, species synchrony and asynchrony across three time periods




# library(devtools)
# install_github("chengjr1009/iStay")

library(iSTAY)
library(tidyverse)
library(lmerTest)

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
  scale_color_manual(values = c("#EA0000","#0066CC","#faa307")) +
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
  scale_color_manual(values = c("#EA0000","#0066CC","#faa307")) +
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

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
    sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
    sub_df
  })
  names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
  out <- multi_beta_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
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
             aes(x = sowndiv, y = beta, color = year), size = 2, alpha = 0.3, show.legend = FALSE) +
  geom_line(data = plotdata3, aes(x = sowndiv, y = predicted, linetype = significance, color = year), linewidth = 1.2) +
  geom_text(data = plotdata_text3,
            aes(x = 2.55, y = 0.545, label = slope_text, color = year,
                hjust = rep(-c(0,0,0), 3),
                vjust = c(rep(c(-2, 0, 2), each = 3))), size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("#EA0000","#0066CC","#faa307")) +
  scale_linetype_manual(values = c("solid","dashed"), drop = FALSE, labels = c("Significant", "Non-significant")) +
  scale_x_continuous(breaks = c(1, 2, 3, 4), labels = c(2, 4, 8, 16)) +
  facet_grid(~ version) +
  guides(linetype = guide_legend(keywidth = 3.1)) +
  labs(linetype = "", x = "Number of species (log2 scale)", y = "Beta stability") + 
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

# ggsave("Population Beta 0512.png", width = 9.5, height = 5)


# Synchroy (76) -----------------------------------------------------------

Jena_10yr_window_result <- lapply(year_windows, function(years){
  Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
    sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
    sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
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
  scale_color_manual(values = c("#EA0000","#0066CC","#faa307")) +
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
  scale_color_manual(values = c("#EA0000","#0066CC","#faa307")) +
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



# # Effect time relationship ------------------------------------------------
# 
# struct_plot <- data.frame(plot = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 1),
#                           block = map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 2),
#                           sowndiv = as.numeric(map_chr(names(Data_Jena_76_metapopulations) |> str_split("_"), 3)))
# 
# blocksowndiv_spe <- struct_plot |> select(block, sowndiv)
# colnames(blocksowndiv_spe) <- c("block", "sowndiv")
# blocksowndiv_spe$sowndiv <- log2(as.numeric(blocksowndiv_spe$sowndiv))
# 
# # Community Stability(76) ----------------------------------------------------
# 
# valid_start_years <- c(2003, 2005:2015)
# year_windows <- lapply(valid_start_years, function(start_year) {
#   years <- start_year:(start_year + 9)
#   if (start_year == 2003) {
#     years <- c(2003, 2005:2013)
#   } else if (2004 %in% years) {
#     return(NULL)
#   }
#   return(as.character(years))
# })
# year_windows <- Filter(Negate(is.null), year_windows)
# 
# Jena_10yr_window_result <- lapply(year_windows, function(years){
#   Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
#     sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
#     sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
#     sub_df
#   })
#   names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
#   out <- multi_gamma_diffk(Data_Jena_76_metapopulations_10, blocksowndiv_spe)
#   return(out)
# })
# 
# year_labels <- sapply(year_windows, function(yrs){
#   paste0(min(yrs), "-", max(yrs))
# })
# year_labels[1] <- "2003, 2005-2013"
# 
# slope <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |> 
#     mutate(year = year_labels[i],
#            starting_year = valid_start_years[i],
#            slope_num = parse_number(slope_text))
# }))
# 
# slope |> 
#   ggplot() +
#   geom_line(aes(x = starting_year, y = slope_num, color = version, linetype = version), linewidth = 1) +
#   geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
#   scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
#   scale_color_manual(values = c("#386641", "#6a994e", "#a7c957")) +
#   labs(y = "Effect of species richness", x = "Starting year of sliding 10-year windows") +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
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
# 
# # ggsave("Effect of diversity on community stability.png", width = 8, height = 6)
# 
# 
# # Alpha -------------------------------------------------------------------
# 
# Jena_10yr_window_result <- lapply(year_windows, function(years){
#   Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
#     sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
#     sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
#     sub_df
#   })
#   names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
#   out <- multi_alpha_diffk(Data_Jena_76_metapopulations_10, blocksowndiv_spe)
#   return(out)
# })
# 
# slope <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |> 
#     mutate(year = year_labels[i],
#            starting_year = valid_start_years[i],
#            slope_num = parse_number(slope_text))
# }))
# 
# slope |> 
#   ggplot() +
#   geom_line(aes(x = starting_year, y = (slope_num), color = version, linetype = version), linewidth = 1) +
#   geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
#   scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
#   scale_color_manual(values = c("#386641", "#6a994e", "#a7c957")) +
#   labs(y = "Effect of species richness", x = "Starting year of sliding 10-year windows") +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
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
# 
# # ggsave("Effect of diversity on population stability.png", width = 8, height = 6)
# 
# # Beta -------------------------------------------------------------------
# 
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
# slope <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |> 
#     mutate(year = year_labels[i],
#            starting_year = valid_start_years[i],
#            slope_num = parse_number(slope_text))
# }))
# 
# slope |> 
#   ggplot() +
#   geom_line(aes(x = starting_year, y = (slope_num), color = version, linetype = version), linewidth = 1) +
#   geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
#   scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
#   scale_color_manual(values = c("#386641", "#6a994e", "#a7c957")) +
#   labs(y = "Effect of species richness", x = "Starting year of sliding 10-year windows") +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
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
# 
# # ggsave("Effect of diversity on beta stability.png", width = 8, height = 6)
# 
# 
# # Synchrony -------------------------------------------------------------------
# 
# Jena_10yr_window_result <- lapply(year_windows, function(years){
#   Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
#     sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
#     sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
#     sub_df
#   })
#   names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
#   out <- multi_syn_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
#   return(out)
# })
# 
# slope <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |> 
#     mutate(year = year_labels[i],
#            starting_year = valid_start_years[i],
#            slope_num = parse_number(slope_text))
# }))
# 
# slope |> 
#   ggplot() +
#   geom_line(aes(x = starting_year, y = (slope_num), color = version, linetype = version), linewidth = 1) +
#   geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
#   scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
#   scale_color_manual(values = c("#386641", "#6a994e", "#a7c957")) +
#   labs(y = "Effect of species richness", x = "Starting year of sliding 10-year windows") +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
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
# 
# # ggsave("Effect of diversity on population synchrony.png", width = 8, height = 6)
# 
# 
# 
# # Asynchrony -------------------------------------------------------------------
# 
# Jena_10yr_window_result <- lapply(year_windows, function(years){
#   Data_Jena_76_metapopulations_10 <- lapply(Data_Jena_76_metapopulations, function(df){
#     sub_df <- df[, years, drop = FALSE]      # 取出這個 window 的欄位
#     sub_df[sub_df == 0] <- 1e-15             # 把 0 換成 10^(-15)
#     sub_df
#   })
#   names(Data_Jena_76_metapopulations_10) <- names(Data_Jena_76_metapopulations)
#   out <- multi_asyn_diffk(Data_Jena_76_metapopulations_10[names], blocksowndiv_spe2)
#   return(out)
# })
# 
# slope <- do.call(rbind, lapply(seq_along(Jena_10yr_window_result), function(i){
#   Jena_10yr_window_result[[i]]$slope_text |>
#     filter(block == "Total") |> 
#     mutate(year = year_labels[i],
#            starting_year = valid_start_years[i],
#            slope_num = parse_number(slope_text))
# }))
# 
# slope |> 
#   ggplot() +
#   geom_line(aes(x = starting_year, y = slope_num, color = version, linetype = version), linewidth = 1) +
#   geom_point(aes(x = starting_year, y = slope_num, color = version), size = 3) +
#   scale_x_continuous(breaks = c(2003, 2005, 2007, 2009, 2011, 2013, 2015)) +
#   scale_color_manual(values = c("#386641", "#6a994e", "#a7c957")) +
#   labs(y = "Effect of species richness", x = "Starting year of sliding 10-year windows") +
#   guides(linetype = guide_legend(keywidth = 3.1)) +
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
# 
# # ggsave("Effect of diversity on population asynchrony.png", width = 8, height = 6)
# 
# 
# 
# 

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
  
  # 這邊是在篩選weight大但是stab小的
  # 將每個plot的weight由大排到小，如果stab小於mean_alpha則選取
  # 但是更好的做法應該是，選取每個plot每個物種的stab小於alpha，老師問再改
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
    geom_point(data = temp |> filter(species == "Ono.vic"), aes(x = weight, y = Stab, color = species), alpha = 0.5, size = 5) +
    geom_line(data = fit_slope(q = q)[[2]], aes(x = weight, y = predicted, group = number_of_species), color = "#00BFC4", linewidth = 1) +
    geom_text(data = fit_slope(q = q)[[1]], 
              aes(x = 0.75, y = 0.09, group = number_of_species, label = paste0("Slope = ", sprintf("%.4f", Slope))), color = "#00BFC4",
              size = 4, key_glyph = draw_key_path) +
    geom_text(data = temp_ono, 
              aes(x = 0.61, y = 0.03, group = number_of_species, label = text), color = "#F8766D",
              size = 4, key_glyph = draw_key_path) +
    ylim(c(0, 1.1)) +
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    # labs(title = paste0("q = ", q), x = "Weight", y = "Stability", color = "Speices") +
    labs(title = paste0("q = ", q), x = "Relative biomass", y = "Stability", color = "Speices") +
    facet_wrap(~ number_of_species, nrow = 1) +
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
          axis.title.y = element_text(hjust = 0.5, size = 14))
  
}

gg_alpha_wrt_weight_group_by_species(data = Jena_tidy |> filter(number_of_species != 1), q = 1, threshold = 1 / 6)
gg_alpha_wrt_weight_group_by_species(data = Jena_tidy |> filter(number_of_species != 1), q = 2, threshold = 1 / 6)



