library(ggplot2)
library(limma)
library(data.table)
library(tidyr)
library(plyr)
library(dplyr)
library(pscl)
#library(ggrepel)
library(vsn)
library(drc)
library(sandwich)
library(lmtest)
library(ggVennDiagram)

# data pre-analysis ----
pg <- fread("data.csv")
metadata <- read.table("metadata.txt",
                       header = TRUE, sep = "\t")
metadata <- metadata[, -ncol(metadata)]

log2norm.df <- pg[ ,7:ncol(pg)]
log2norm.df <- as.matrix(log2norm.df)
rownames(log2norm.df) <- pg$genes
# remove samples with less than 2000 IDs for proteome data, delete those samples in metadata as well
metadata <- metadata[which(colSums(!is.na(log2norm.df)) >= 2000),]
log2norm.df <- log2norm.df[, colSums(!is.na(log2norm.df)) >= 2000]


# check data distribution and missing proportion
hist(log2norm.df, breaks = 100)
sum(is.na(log2norm.df)) / (nrow(log2norm.df) * ncol(log2norm.df))
# density overlay plot of each sample
df_long <- as.data.frame(log2norm.df) %>% pivot_longer(cols = everything(),
                                                       names_to = "Sample_ID",
                                                       values_to = "Intensity")

ggplot(df_long, aes(x = Intensity, group = Sample_ID)) +
  geom_density(alpha = 0.05, color = "steelblue", size = 0.3) +
  theme_minimal() +
  labs(title = "Density overlay",
       x = "log2 intensity",
       y = "density ") +
  theme(legend.position = "none")

# batch correction (not used for dose-response) ----
log2norm.bc.df <- removeBatchEffect(log2norm.df, batch = as.factor(metadata$Donor))
# merge with pg 
pg_sub <- pg[match(rownames(log2norm.bc.df), pg$genes), 1:6]

log2norm.bc.df_merged <- cbind( pg_sub, log2norm.bc.df)
write.csv(log2norm.bc.df_merged, file = "C:/Users/ml132260/OneDrive - GSK/Documents/curve fitting/Chloe_data/TAK242/batch_corrected_Px.csv")

# data imputation ----
# ND imputation
impute_normal <- function(object, width=0.3, downshift=1.8, seed=100) {
  
  if (!is.matrix(object)) object <- as.matrix(object)
  mx <- max(object, na.rm=TRUE)
  mn <- min(object, na.rm=TRUE)
  if (mx - mn > 20) warning("Please make sure the values are log-transformed")
  
  set.seed(seed)
  object <- apply(object, 2, function(temp) {
    temp[!is.finite(temp)] <- NA
    temp_sd <- stats::sd(temp, na.rm=TRUE)
    temp_mean <- mean(temp, na.rm=TRUE)
    shrinked_sd <- width * temp_sd   # shrink sd width
    downshifted_mean <- temp_mean - downshift * temp_sd   # shift mean of imputed values
    n_missing <- sum(is.na(temp))
    temp[is.na(temp)] <- stats::rnorm(n_missing, mean=downshifted_mean, sd=shrinked_sd)
    temp
  })
  return(object)
}

log2norm.imp.df <- impute_normal(log2norm.df)
hist(log2norm.imp.df, breaks = 100)

# y is log intensity(flat vs LL4) ----
# model comparison
N <- nrow(log2norm.imp.df)
result_all <- matrix(NA, nrow = N, ncol = 6, dimnames = list(rownames(log2norm.imp.df), 
                                                             c("delta_AIC", "slop", 
                                                               "Lower_limit", "Upper_limit",
                                                               "ED50", "R2")))
result_all <- as.data.frame(result_all)

conc_levels <- c(10.000, 3.333, 1.111, 0.370, 0.123, 0.041, 0.014, 0.005, 0.002, 
                 0.001, 0.0002, 1e-8 )
concentration <- c(rep(conc_levels, 3), rep(1e8, 3))
concentration <- rep(concentration, 3)

# logFC matrix
logFC.df <- log2norm.imp.df - rowMeans(log2norm.imp.df[, which(concentration == 1e8)])


for (i in 1:N){
  
  prot_i <- logFC.df[i, (concentration != 1e8 & concentration != 1e-8)]
  prot_i.df <- data.frame(response = as.numeric(prot_i), dose = concentration[(concentration != 1e8 & concentration != 1e-8)])
  
  prot_i.fit.null <- tryCatch(
    lm(response ~ 1, data = prot_i.df),
    error = function(e) {
      message("Skipping row ", i, " due to error: ", e$message)
      return(NULL)
    }
  )
  
  prot_i.fit.LL.4 <- tryCatch(
    drm(response ~ dose, data = prot_i.df, fct = LL.4(),
        start = c(
          b = 1,
          c = min(prot_i.df$response),
          d = max(prot_i.df$response),
          e = median(concentration))),
    error = function(e) {
      message("Skipping row ", i, " due to error: ", e$message)
      return(NULL)
    }
  )
  
  if (is.null(prot_i.fit.null) | is.null(prot_i.fit.LL.4)) next
  
  rss <- sum(residuals(prot_i.fit.LL.4)^2)
  tss <- sum(residuals(prot_i.fit.null)^2)
  
 R2_i <- 1 - rss / tss
  
  
  delta_AIC_i <- AIC(prot_i.fit.null) - AIC(prot_i.fit.LL.4)
  result_all$delta_AIC[i] <- delta_AIC_i
  result_all$slop[i] <- summary(prot_i.fit.LL.4)[["coefficients"]]["b:(Intercept)", 1]
  result_all$Lower_limit[i] <- summary(prot_i.fit.LL.4)[["coefficients"]]["c:(Intercept)", 1]
  result_all$Upper_limit[i] <- summary(prot_i.fit.LL.4)[["coefficients"]]["d:(Intercept)", 1]
  result_all$ED50[i] <- summary(prot_i.fit.LL.4)[["coefficients"]]["e:(Intercept)", 1]
  result_all$R2[i] <-R2_i
}

write.csv(result_all, "result.csv", row.names = TRUE)
pg_sub <- pg[match(rownames(logFC.df), pg$genes), 1:6]
logFC.df_merged <- cbind(pg_sub, logFC.df)
write.csv(logFC.df_merged, "logFC.merged.csv", row.names = TRUE)

example_prot.df <- data.frame(response = as.numeric(log2norm.imp.df["CXCL3",]),
                              concentration = concentration)
example_prot.fit <- drm(response ~ concentration, data = example_prot.df, fct = LL.4(),
                        start = c(
                          b = 1,
                          c = min(log2norm.imp.df["CXCL3", ]),
                          d = max(log2norm.imp.df["CXCL3", ]),
                          e = median(concentration)
                        )
)
plot(example_prot.fit, broken = TRUE, type = "all", xlim = c(1e-8,1e8), log = "x")
summary(example_prot.fit)

