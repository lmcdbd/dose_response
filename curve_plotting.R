library(dplyr)
library(tidyr)
library(drc)
library(ggplot2)
library(data.table)
# import data ----
dose_response_R848_Sx <- fread("dose_response_R848_Sx.csv")
dose_response_R848_Px <- fread("dose_response_R848_Px.csv")
dose_response_TAK242_Sx <- fread("dose_response_TAK242_Sx.csv")
dose_response_TAK242_Px <- fread("dose_response_TAK242_Px.csv")

pg <- fread("data.csv")
metadata <- read.table("metadata.txt",
                       header = TRUE, sep = "\t")

metadata <- metadata[, -ncol(metadata)]
# TAK242 proteome
metadata <- metadata[-79, ]

log2norm.df <- pg[ ,7:ncol(pg)]
log2norm.df <- as.matrix(log2norm.df)
rownames(log2norm.df) <- pg$genes
# remove samples with less than 2000 IDs, delete those samples in metadata as well
metadata <- metadata[which(colSums(!is.na(log2norm.df)) >= 2000),]
log2norm.df <- log2norm.df[, colSums(!is.na(log2norm.df)) >= 2000]

# imputation
source("curve fitting/ND_imputation.R")
log2norm.imp.df <- impute_normal(log2norm.df)
hist(log2norm.imp.df, breaks = 100)

# R848 secretome
conc_levels <- c(10.000, 3.333, 1.111, 0.370, 0.123, 0.041, 0.014, 0.005, 0.002, 
                 0.001, 0.0002, 1e-8 )
concentration <- c(rep(conc_levels, 3), rep(1e8, 3))
concentration_all <- c(rep(concentration, 3))

# R848 proteome
conc_levels <- c(10.000, 3.333, 1.111, 0.370, 0.123, 0.041, 0.014, 0.005, 0.002, 
                 0.001, 0.0002, 1e-8 )
concentration <- c(rep(conc_levels, 3), rep(1e8, 12))
concentration_all <- rep(concentration, 3)
concentration_all <- concentration_all[-c(110,113,121)]

# TAK242 secretome
conc_levels <- c(10.000, 3.333, 1.111, 0.370, 0.123, 0.041, 0.014, 0.005, 0.002, 
                 0.001, 0.0002, 1e-8 )
concentration <- c(rep(conc_levels, 3), rep(1e8, 4))
concentration_all <- rep(concentration, 3)

# TAK242 proteome
conc_levels <- c(10.000, 3.333, 1.111, 0.370, 0.123, 0.041, 0.014, 0.005, 0.002, 
                 0.001, 0.0002, 1e-8 )
concentration <- c(rep(conc_levels, 3), rep(1e8, 3))
concentration_all <- c(rep(concentration, 3))

# D-R curve for individual proteins ----
# use logFC as y axis
col_index <- which(metadata$Donor == "006" )# get column index 
D006_df <- log2norm.imp.df[, col_index]
D006_logFC.df <- D006_df - rowMeans(D006_df[, 37:40], na.rm = TRUE)

col_index <- which(metadata$Donor == "S2" )# get column index 
S2_df <- log2norm.imp.df[, col_index]
S2_logFC.df <- S2_df - rowMeans(S2_df[, 37:40], na.rm = TRUE)

col_index <- which(metadata$Donor == "M3" )# get column index 
M3_df <- log2norm.imp.df[, col_index]
M3_logFC.df <- M3_df - rowMeans(M3_df[, 37:40], na.rm = TRUE)

All_logFC.df <- log2norm.imp.df - rowMeans(log2norm.imp.df[, which(concentration_all == 1e8)])


# plot D-R curve for one protein
concentration_M <- 1e-6 * concentration
concentration_all_M <- 1e-6 * concentration_all

# R848 proteome 
concentration_M <- 1e-6 * concentration
concentration_all_M <- 1e-6 * concentration_all
concentration_M_S2 <- concentration_M[-c(14, 17, 25)]

# loop through list 
# R848 secretome list
R848_protein_list_Sx <- 
  c(
    "ADM", "CCL15", "CCL2", "CCL3", "CCL3L1",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8",
    "ICAM1",
    "IL12B", "IL1B", "IL1RN", "IL6",
    "TNF", "TNFAIP6", "TNFSF9"
  )



R848_protein_list_Px <- c(
  "ARHGDIA", "C1QBP", "CCDC124", "CCL20", "CCL3", "CD83",
  "CSF1R", "CXCL1", "CXCL2", "CXCL3", "CXCL8",
  "DNAJB6", "EHD1", "ENG", "ICAM1", "IFNGR1",
  "IL1A", "IL1B", "IRAK2",
  "MERTK", "MT.CO2", "MYO1F",
  "NFKB2", "NFKBIZ", "NRP2",
  "OLR1",
  "PDCD4", "PDIA3", "PLXDC2", "PTGS2",
  "RND3", "RPS19", "RTN4",
  "SDC4", "SLAMF7", "SLC18B1", "SLCO2B1", "SQSTM1",
  "SRC", "SRXN1", "STIP1",
  "THBS1", "TNF", "TNFAIP2", "TNFAIP6",
  "TNIP1", "TNIP3",
  "TRAF1", "TRIP10", "TWF2",
  "VPS35"
)


TAK242_protein_list_Sx <- 
  c(
    "ADM", "CCL2", "CCL20", "CCL3", "CCL3L1", "CCL8",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10",
    "ICAM1",
    "IFIT2",
    "IL1B", "IL1RN", "IL6",
    "INHBA",
    "MMP1", "MMP19",
    "PLAUR",
    "PTX3",
    "SDC4", "SRGN",
    "TIMP1",
    "TNF", "TNFAIP6", "TNFSF9"
  )

  

TAK242_protein_list_Px <-  c(
  "ATP2A2", "CCL2", "CCL20", "CCL3", "CCL3L1", "CCL5",
  "CD163", "CD40", "CD83",
  "CKB", "CLIC4", "CST3",
  "CSF1R", "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10",
  "DNAJB6", "EHD1", "EIF1", "ENG", "ETV3",
  "FCGR3A", "FOSL2",
  "GTF2B",
  "HAVCR2", "HELZ2",
  "ICAM1", "IFIH1", "IFIT1", "IFIT2", "IFIT3", "IL1A", "IL1B", "IL6",
  "ISG15",
  "JUNB",
  "LCP2",
  "NAMPT", "NFKB2", "NFKBIZ", "NRP2",
  "OASL", "OTUD1",
  "PFKFB3", "PLAUR", "PLXDC2", "PMAIP1", "PTGS2",
  "REL", "RIN2", "RIPK1", "RND3", "RSAD2",
  "SDC4", "SEMA6B", "SLAMF7", "SLCO2B1", "SRC", "STAT3", "STX11",
  "TNF", "TNFAIP2", "TNFAIP3", "TNFAIP6", "TNFSF9",
  "TRAF1", "TRIP10",
  "VSIG4",
  "WARS1", "WTAP",
  "ZC3H12A", "ZC3HAV1", "ZFP36"
)
  
  
  

R848_protein_list <-   
  c(
    "ADM", "CCL2", "CCL20", "CCL3", "CCL3L1", "CCL8",
    "CXCL1", "CXCL2", "CXCL3", "CXCL8", "CXCL10",
    "ICAM1",
    "IFIT2",
    "IL1B", "IL1RN", "IL6",
    "INHBA",
    "MMP1", "MMP19",
    "PLAUR",
    "PTX3",
    "SDC4", "SRGN",
    "TIMP1",
    "TNF", "TNFAIP6", "TNFSF9"
  )


# note that for some proteins the fit to all samples converge but may not converge for some donors
# For TAK242, replace A27 with M3

plot_list <- list()


for (i in 1:length(R848_protein_list)){
  failed <- FALSE
  
  if (R848_protein_list[i] == "CD40") {    # only for TAK242 Proteome CD40 
    ylim_use <- c(-0.1, 1)
  } else {
    ylim_use <- NULL
  }
  
  
  tryCatch({
  
  df1 <- data.frame(response = as.numeric(S2_logFC.df[R848_protein_list[i], (concentration_M != 1e2 & concentration_M != 1e-14)]),
                    concentration = concentration_M[concentration_M != 1e2 & concentration_M != 1e-14])
  df2 <- data.frame(response = as.numeric(M3_logFC.df[R848_protein_list[i], (concentration_M != 1e2 & concentration_M != 1e-14)]),
                    concentration = concentration_M[concentration_M != 1e2 & concentration_M != 1e-14])
  df3 <- data.frame(response = as.numeric(D006_logFC.df[R848_protein_list[i],(concentration_M != 1e2 & concentration_M != 1e-14)]),
                    concentration = concentration_M[concentration_M != 1e2 & concentration_M != 1e-14])
  df4 <- data.frame(response = as.numeric(All_logFC.df[R848_protein_list[i], (concentration_all_M != 1e2 & concentration_all_M != 1e-14)]),
                    concentration = concentration_all_M[concentration_all_M != 1e2 & concentration_all_M != 1e-14])
  
  fit1 <- drm(response ~ concentration, data = df1, fct = LL.4(),
              start = c(
                b = 1,
                c = min(df1$response),
                d = max(df1$response),
                e = median(concentration_all_M))
  )
  fit2 <- drm(response ~ concentration, data = df2, fct = LL.4(),
              start = c(
                b = 1,
                c = min(df2$response),
                d = max(df2$response),
                e = median(concentration_all_M))
  )
  
  fit3 <- drm(response ~ concentration, data = df3, fct = LL.4(),
              start = c(
                b = 1,
                c = min(df3$response),
                d = max(df3$response),
                e = median(concentration_all_M))
  )
  fit4 <- drm(response ~ concentration, data = df4, fct = LL.4(),
              start = c(
                b = 1,
                c = min(df4$response),
                d = max(df4$response),
                e = median(concentration_all_M))
  )
  
  IC50 <- as.numeric(dose_response_TAK242_Sx[dose_response_TAK242_Sx$V1 == R848_protein_list[i], ED50]) * 1e3
  # calculate mean 
  df1_mean <-df1 %>%
    dplyr::group_by(concentration) %>%
    dplyr::summarise(response_mean = mean(response, na.rm = TRUE))
  df2_mean <-df2 %>%
    dplyr::group_by(concentration) %>%
    dplyr::summarise(response_mean = mean(response, na.rm = TRUE))
  df3_mean <-df3 %>%
    dplyr::group_by(concentration) %>%
    dplyr::summarise(response_mean = mean(response, na.rm = TRUE))
  df4_mean <-df4 %>%
    dplyr::group_by(concentration) %>%
    dplyr::summarise(response_mean = mean(response, na.rm = TRUE))
  
  # plot 95% confidence interval
  
  xseq <- exp(seq(log(min(df4$concentration) * 1e-2),  # extend the curve 
                  log(max(df4$concentration) * 10),
                  length = 200))
  
  newdata4 <- data.frame(concentration = xseq)
  
  pm4 <- predict(
    fit4,                       # LL.4 model fitted to df4_mean
    newdata = newdata4,
    interval = "confidence"
  )
  
  newdata4$p    <- pm4[,1]
  newdata4$pmin <- pm4[,2]
  newdata4$pmax <- pm4[,3]
  
  newdata123 <- data.frame(concentration = xseq)
  
  newdata123$fit1 <- predict(fit1, newdata123)
  newdata123$fit2 <- predict(fit2, newdata123)
  newdata123$fit3 <- predict(fit3, newdata123)
  
  # add curve information
  curves_bg <- newdata123 %>%
    pivot_longer(
      cols = c(fit1, fit2, fit3),
      names_to = "Curve",
      values_to = "Response"
    )
  
  curves_bg$Curve <- recode(curves_bg$Curve,
                            fit1 = "Donor 3",
                            fit2 = "Donor 2",
                            fit3 = "Donor 1")
  
  
  curve_df4 <- data.frame(
    concentration = newdata4$concentration,
    Response = newdata4$p,
    Curve = "All (mean ± 95% CI)"
  )
  
  all_curves <- rbind(curves_bg, curve_df4)
  
  p <- ggplot() +
    
    geom_line(
      data = curves_bg,
      aes(x = concentration, y = Response, color = Curve),
      linewidth = 1.5,
      alpha = 0.4
    ) +
    
    geom_ribbon(
      data = newdata4,
      aes(x = concentration, ymin = pmin, ymax = pmax),
      fill = "grey60",
      alpha = 0.35
    ) +
    
    geom_line(
      data = curve_df4,
      aes(x = concentration, y = Response, color = Curve),
      linewidth = 1.5
    ) +
    
    geom_point(
      data = df4_mean,
      aes(x = concentration, y = response_mean),
      size = 2.5,
      color = "black"
    ) +
    
    geom_point(
      data = df3_mean,
      aes(x = concentration, y = response_mean),
      size = 2.5,
      color = "#B85FB8",
      alpha = 0.4
    ) +
    geom_point(
      data = df2_mean,
      aes(x = concentration, y = response_mean),
      size = 2.5,
      color = "#00BFBF",
      alpha = 0.4
    ) +
    geom_point(
      data = df1_mean,
      aes(x = concentration, y = response_mean),
      size = 2.5,
      color = "#A67C52",
      alpha = 0.4
    ) +
    
    scale_x_log10(
      limits = c(1e-10, 1e-5),
      breaks = 10^(-10:-5),
      labels = -10:-5
    ) +
    
  
    scale_color_manual(
      name = NULL,
      values = c(
        "Donor 3" = "#A67C52",
        "Donor 2" = "#00BFBF",
        "Donor 1" = "#B85FB8",
        "All (mean ± 95% CI)" = "black"
      )
    ) +
    
    labs(
      title = R848_protein_list[i],   # ✅ dynamic title
      x = expression(log[10]~"[Inhibitor], M"),
      y = expression(log[2]~"(fold change)")
    ) +
    
    theme_classic() +
    
    theme(aspect.ratio = 1) +
    
    theme(
      # title formatting
      plot.title = element_text(
        face = "bold",
        size = 14,
        hjust = 0.5   # center title
      ),
      
      axis.line = element_line(linewidth = 1.2, color = "black"),

      axis.ticks = element_line(linewidth = 1.2),
      axis.ticks.length = unit(0.25, "cm"),
      
      
      axis.title = element_text(face = "bold", size = 13),
      axis.text = element_text(size = 11, color = "black")
    ) +
    
    annotate(
      "text",
      x = min(df4$concentration) * 2,
      y = -Inf,
      label = paste0("EC50 = ", round(IC50, 2), " nM"),
      hjust = 0.1,
      vjust = -0.5,
      size = 4,
      fontface = "bold"
    ) +
    
    coord_cartesian(
      ylim = ylim_use,
      clip = "off"
    ) +
    
  
  theme(legend.position = "none") 
  
  plot_list[[length(plot_list) + 1]] <- p
  
  
  }, error = function(e) {
    message("Protein ", R848_protein_list[i], " failed: ", e$message)
    failed <<- TRUE   # assign outside tryCatch
  })
  
  if (failed) next
  if (length(plot_list) %% 20 == 0) {
    
    idx_start <- length(plot_list) - 19
    idx_end <- length(plot_list)
    
    combined_plot <- wrap_plots(
      plot_list[idx_start:idx_end],
      ncol = 4,
      nrow = 5
    )
    
    ggsave(
      filename = paste0(
        "dose_response_plots/Supplementary/Panel_",
        idx_start, "_to_", idx_end, ".png"
      ),
      plot = combined_plot,
      width = 12,   
      height = 12
    )
  }
  
}

remaining <- length(plot_list) %% 20

if (remaining != 0) {
  
  idx_start <- length(plot_list) - remaining + 1
  idx_end <- length(plot_list)
  
  combined_plot <- wrap_plots(
    plot_list[idx_start:idx_end],
    ncol = 4,
    nrow = ceiling(remaining / 4)
  )
  
  ggsave(
    paste0("dose_response_plots/Supplementary/Panel_", idx_start, "_to_", idx_end, ".png"),
    combined_plot,
    width = 12,
    height = 12
  )
}
