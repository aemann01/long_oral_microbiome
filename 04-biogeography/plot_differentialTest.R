library(ggplot2)
library(phyloseq)
library(corncob)
library(magrittr)

load("../02-diversity/master_phyloseq.RData")
ps.dat
head(sample_data(ps.dat))

# glom to species level, only doing incisors/premolars for testing purposes
ps.dat.sp <- ps.dat %>%
                phyloseq::subset_samples(tooth_type == "incisor" | tooth_type == "premolar") %>%
                tax_glom("V8")
ps.dat.sp

da_analysis <- differentialTest(formula = ~ tooth_type,
                               phi.formula = ~ 1,
                               formula_null = ~ 1,
                               phi.formula_null = ~ 1,
                               test = "Wald",
                               boot = FALSE,
                               data = ps.dat.sp,
                               fdr_cutoff = 0.01)

save.image("test.R")

level <- c("V8")
x <- da_analysis

signif_taxa <- otu_to_taxonomy(x$significant_taxa, x$data, level = level)
signif_taxa <- paste0(signif_taxa, " (", x$significant_taxa, ")") # NOTE: ASV label no longer informative after glom

var_per_mod <- length(x$restrictions_DA) + length(x$restrictions_DV)
total_var_count <- length(signif_taxa) * var_per_mod
# initialize empty dataframe 
df <- as.data.frame(matrix(NA, nrow = total_var_count, ncol = 5))
colnames(df) <- c("x", "xmin", "xmax", "taxa", "variable")

qval <- stats::qnorm(.975)
restricts_mu <- attr(x$restrictions_DA, "index")
restricts_phi <- attr(x$restrictions_DV, "index")

# loop through models to populate dataframe
count <- 1
for (i in 1:length(x$significant_models)) {
  # Below from print_summary_bbdml, just to get coefficient names
  tmp <- x$significant_models[[i]]
  coefs.mu <- tmp$coefficients[1:tmp$np.mu,, drop = FALSE]
  rownames(coefs.mu) <- paste0(substring(rownames(coefs.mu), 4), " Differential Abundance")
  coefs.mu <- coefs.mu[restricts_mu,, drop = FALSE]

  coefs.phi <- tmp$coefficients[(tmp$np.mu + 1):nrow(tmp$coefficients),, drop = FALSE]
  rownames(coefs.phi) <- paste0(substring(rownames(coefs.phi), 5), " Differential Variability")
  coefs.phi <- coefs.phi[restricts_phi - tmp$np.mu,, drop = FALSE]

  coefs <- rbind(coefs.mu, coefs.phi)
  for (j in 1:var_per_mod) {
    df[count, 1:3] <- c(coefs[j, 1], coefs[j, 1] - qval * coefs[j, 2],
                      coefs[j, 1] + qval * coefs[j, 2])
    df[count, 4:5] <- c(signif_taxa[i], rownames(coefs)[j])
    count <- count + 1
    }
}

head(df)

pdf("test.pdf")
ggplot(df, aes(x = x, y = taxa)) +
  geom_vline(xintercept = 0, color = "gray50", lty = "dashed",
                          alpha = 0.75, lwd = 1) +
  geom_point(show.lengend=FALSE) +
  geom_errorbarh(aes(xmin = xmin, xmax = xmax, colour = xmin<=0 & xmax <= 0| xmin>=0 & xmax >= 0), height = .3, show.lengend=FALSE) +
  theme_bw() +
  facet_wrap(~variable, scales = "free_x", nrow = 1) +
  labs(title = "", x = "", y = "Taxa") +
  scale_y_discrete(limits = rev(df$taxa)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






plot.differentialTest <- function(x, level = NULL, data_only = FALSE, ...) {
  signif_taxa <- x$significant_taxa
  if ("phyloseq" %in% class(x$data)) {
    if (!is.null(x$data@tax_table)) {
      signif_taxa <- otu_to_taxonomy(signif_taxa, x$data, level = level)
      if (length(unique(signif_taxa)) != length(unique(x$significant_taxa))) {
        # Make sure if repeated taxa add unique otu identifiers
        signif_taxa <- paste0(signif_taxa, " (", x$significant_taxa, ")")
      }
    }
  }
  if (length(x$significant_models) != 0) {
    var_per_mod <- length(x$restrictions_DA) + length(x$restrictions_DV)
    total_var_count <- length(signif_taxa) * var_per_mod
    df <- as.data.frame(matrix(NA, nrow = total_var_count, ncol = 5))
    colnames(df) <- c("x", "xmin", "xmax", "taxa", "variable")
    qval <- stats::qnorm(.975)
    restricts_mu <- attr(x$restrictions_DA, "index")
    restricts_phi <- attr(x$restrictions_DV, "index")

    count <- 1
    for (i in 1:length(x$significant_models)) {

      # Below from print_summary_bbdml, just to get coefficient names
      tmp <- x$significant_models[[i]]
      coefs.mu <- tmp$coefficients[1:tmp$np.mu,, drop = FALSE]
      rownames(coefs.mu) <- paste0(substring(rownames(coefs.mu), 4), "\nDifferential Abundance")
      coefs.mu <- coefs.mu[restricts_mu,, drop = FALSE]

      coefs.phi <- tmp$coefficients[(tmp$np.mu + 1):nrow(tmp$coefficients),, drop = FALSE]
      rownames(coefs.phi) <- paste0(substring(rownames(coefs.phi), 5), "\nDifferential Variability")
      coefs.phi <- coefs.phi[restricts_phi - tmp$np.mu,, drop = FALSE]

      coefs <- rbind(coefs.mu, coefs.phi)
      for (j in 1:var_per_mod) {
        df[count, 1:3] <- c(coefs[j, 1], coefs[j, 1] - qval * coefs[j, 2],
                        coefs[j, 1] + qval * coefs[j, 2])
        df[count, 4:5] <- c(signif_taxa[i], rownames(coefs)[j])
        count <- count + 1
      }
    }
    if (!(data_only)) {
    # global variables warning suppression
    taxa <- xmin <- xmax <- NULL

    ggplot2::ggplot(df, ggplot2::aes(x = x, y = taxa)) +
      ggplot2::geom_vline(xintercept = 0, color = "gray50", lty = "dashed",
                          alpha = 0.75, lwd = 1) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = xmin, xmax = xmax, colour = xmin<=0 & xmax <= 0| xmin>=0 & xmax >= 0), height = .3) +
      ggplot2::theme_bw() +
      ggplot2::facet_wrap(~variable, scales = "free_x", nrow = 1) +
      ggplot2::labs(title = "", x = "", y = "Taxa") +
      ggplot2::scale_y_discrete(limits = rev(df$taxa)) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else {
    return(df)
  }
  } else {
    message("No taxa were found to be significantly different using your model specification. \nPlease verify that your formulas are correctly specified.")
  }
}