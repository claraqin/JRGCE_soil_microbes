# R code for Qin et al., Journal of Ecology, accepted Aug. 28, 2019

# Figure 2 and related statistics on microbial diversity.
# Effects of global change treatments, soil properties, 
# and plant community properties on observed richness 
# of soil microbial communities. 

# For this script to work, you must first run 0_setup.R
# in the same R session.

library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(lemon)
theme_set(theme_bw())
standardize <- function(x) { (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

# Convenience function for plotting model coefficients, 
# either just observed richness or just Shannon div.
plot_lme_coefs_richness <- function(coefs_table, levels_display, diversity_metric) {
  cbbPalette <- c("#000000", "#FF69B4", "#009E73") #"#0072B2")
  coefs_table %>%
    filter(variable == diversity_metric) %>% # Observed richness only
    mutate(model_sig = interaction(model, sig)) %>%
    ggplot(aes(x=treatment, y=std.coef, col=model, dodge=model, shape=sig, alpha=sig)) +
    facet_rep_grid(community~., repeat.tick.labels=TRUE) +
    scale_shape_manual(values=c(1,16)) + 
    scale_alpha_manual(values=c(0.6,1)) +
    geom_abline(intercept=0, slope=0, col="grey") +
    geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), 
                  position=position_dodge(0.6),
                  width=0) +
    scale_x_discrete(labels = levels_display) +
    scale_y_continuous(breaks=c(seq(-4,4))) +
    geom_point(position=position_dodge(0.6), size=2.5) + 
    scale_colour_manual(labels=c("Main","GC interactions", "Plant-soil"),
                        values=cbbPalette) +
    guides(fill=FALSE) + xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5, size=12)) +
    theme(axis.text.y = element_text(size=14)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          strip.background = element_blank(),
          strip.text.x = element_blank()) +
    NULL
}

# Convenience function to plot lme residuals and print Shapiro non-normality test pvals
check_residuals_normality <- function(model_obj, model_name) {
  # Print results of Shapiro test for non-normality
  shapiro.pval <- shapiro.test(residuals(model_obj))$p.value
  print(paste0("Shapiro test on residuals from ", model_name, " model: ", shapiro.pval))
  # Plot residuals
  # qqnorm(residuals(model_obj), main=paste0(model_name, ",\n", "pval = ", sprintf("%.4f", shapiro.pval)))
  par(mar=c(5,4,1,2)+0.1)
  qqnorm(residuals(model_obj), main=NULL)
  qqline(residuals(model_obj))
}


# 1. GC main model: (no LR test, only plots)
par(mfcol=c(4,3))
coefs_all <- list()
for(var in c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")) {
  sampleData.long %>%
    filter(variable == var,
           complete.cases(.)) %>%
    mutate(value = standardize(value)) ->
    temp
  lme.main <- lmer(value ~ n + co2 + heat + precip + wildfire.2003 + burn.2011 + (1|plt), temp)
  coefs <- as.data.frame(summary(lme.main)$coefficients)
  colnames(coefs) <- c("est", "stderr", "df", "tval", "pval")
  coefs %>%
    mutate(std.coef = est,
           ci.upper = std.coef + qnorm(0.975) * stderr,
           ci.lower = std.coef - qnorm(0.975) * stderr,
           sig = pval < 0.05,
           treatment = rownames(coefs),
           variable = substr(var, 5, 7),
           community = substr(var, 1, 3)) %>%
    filter(treatment != "(Intercept)") %>%
    dplyr::select(community, variable, treatment, sig, std.coef, ci.lower, ci.upper) -> temp_out
  coefs_all <- c(coefs_all, list(temp_out))
  
  # Check for normality of lme residuals
  check_residuals_normality(lme.main, paste(var, "main"))
}
# Collect and save and plot coefficients
coefs_main <- do.call(rbind, coefs_all)
# write.csv(coefs_main, "lme_coefs_main.csv")

# 2. GC interactions model:
# diversity ~ N + C + T + P + NC + NT + NP + CT + CP + TP + NCT + NCP + NTP + CTP + NCTP
coefs_all <- list()
lrtest_results <- list()
for(var in c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")) {
  sampleData.long %>%
    filter(variable == var,
           complete.cases(.)) %>%
    mutate(value = standardize(value)) ->
    temp
  lme.interact <- lmer(value ~ n * co2 * heat * precip + wildfire.2003 + burn.2011 + (1|plt), temp)
  coefs <- as.data.frame(summary(lme.interact)$coefficients)
  colnames(coefs) <- c("est", "stderr", "df", "tval", "pval")
  coefs %>%
    mutate(std.coef = est,
           ci.upper = std.coef + qnorm(0.975) * stderr,
           ci.lower = std.coef - qnorm(0.975) * stderr,
           sig = pval < 0.05,
           treatment = rownames(coefs),
           variable = substr(var, 5, 7),
           community = substr(var, 1, 3)) %>%
    filter(treatment != "(Intercept)") %>%
    dplyr::select(community, variable, treatment, sig, std.coef, ci.lower, ci.upper) -> temp_out
  coefs_all <- c(coefs_all, list(temp_out))
  
  # Check for normality of lme residuals
  check_residuals_normality(lme.interact, paste(var, "interact"))
  
  # In the meantime, do some LR testing
  lme.full.lr <- lmer(value ~ n * co2 * heat * precip + wildfire.2003 + burn.2011 + (1|plt),
                      temp, REML=FALSE)
  lme.main.lr <- lmer(value ~ n + co2 + heat + precip + wildfire.2003 + burn.2011 + (1|plt),
                      temp, REML=FALSE)
  lrtest <- anova(lme.full.lr, lme.main.lr)
  lrtest_results <- c(lrtest_results, 
                      list(c(lrtest[["AIC"]][1], lrtest[["AIC"]][2], lrtest[["Chisq"]][2], lrtest[["Df"]][1], lrtest[["Df"]][2], lrtest[["Pr(>Chisq)"]][2])))
}
# Collect and save lr test results
lrtest_df <- as.data.frame(do.call(rbind, lrtest_results))
colnames(lrtest_df) <- c("AIC_a","AIC_b","Chisq","df1","df2","pval")
rownames(lrtest_df) <- c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")
# write.csv(lrtest_df, "lrtest_interact.csv")

# Collect and save and plot coefficients
coefs_interact <- do.call(rbind, coefs_all)
# write.csv(coefs_interact, "lme_coefs_interact.csv")

# 3. Plant + soil model:
# diversity ~ N + C + T + P + %n + moist. + c/n + ph + npp + p.rich + p.NMDS1 + p.NMDS2 + WF03 + B11
coefs_all <- list()
lrtest_results <- list()
for(var in c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")) {
  sampleData.long %>%
    filter(variable == var,
           complete.cases(.)) %>%
    mutate(value = standardize(value)) %>%
    mutate_at(vars(perc.n:NMDS.Plants.2), standardize) -> temp
  lme.full <- lmer(value ~ n + co2 + heat + precip + wildfire.2003 + burn.2011 + (1|plt) +
                     perc.n + water.content + c.n + ph +
                     npp.std + plant.rich + NMDS.Plants.1 + NMDS.Plants.2, temp)
  coefs <- as.data.frame(summary(lme.full)$coefficients)
  colnames(coefs) <- c("est", "stderr", "df", "tval", "pval")
  coefs %>%
    mutate(std.coef = est,
           ci.upper = std.coef + qnorm(0.975) * stderr,
           ci.lower = std.coef - qnorm(0.975) * stderr,
           sig = pval < 0.05,
           treatment = rownames(coefs),
           variable = substr(var, 5, 7),
           community = substr(var, 1, 3)) %>%
    filter(treatment != "(Intercept)") %>%
    dplyr::select(community, variable, treatment, sig, std.coef, ci.lower, ci.upper) -> temp_out
  coefs_all <- c(coefs_all, list(temp_out))
  
  # Check for normality of lme residuals
  check_residuals_normality(lme.full, paste(var, "plantsoil"))
  
  # In the meantime, do some LR testing
  lme.full.lr <- lmer(value ~ n + co2 + heat + precip + wildfire.2003 + burn.2011 + (1|plt) +
                        perc.n + water.content + c.n + ph +
                        npp.std + plant.rich + NMDS.Plants.1 + NMDS.Plants.2, temp, REML=FALSE)
  lme.main.lr <- lmer(value ~ n + co2 + heat + precip + wildfire.2003 + burn.2011 + (1|plt),
                      temp, REML=FALSE)
  lrtest <- anova(lme.full.lr, lme.main.lr)
  lrtest_results <- c(lrtest_results, 
                      list(c(lrtest[["AIC"]][1], lrtest[["AIC"]][2], lrtest[["Chisq"]][2], lrtest[["Df"]][1], lrtest[["Df"]][2], lrtest[["Pr(>Chisq)"]][2])))
}
# Collect and save lr test results
lrtest_df <- as.data.frame(do.call(rbind, lrtest_results))
colnames(lrtest_df) <- c("AIC_a","AIC_b","Chisq","df1","df2","pval")
rownames(lrtest_df) <- c("ITS.observed","ITS.shannon","BAC.observed","BAC.shannon")
# write.csv(lrtest_df, "lrtest_plantsoil.csv")

# Collect and save and plot coefficients
coefs_plantsoil <- do.call(rbind, coefs_all)
# write.csv(coefs_plantsoil, "lme_coefs_plantsoil.csv")

# 4. Combine all three model types into same plot
levels_order_threemodels <- c("n2","co22","heat2","precip2","wildfire.20032","burn.20112",
                              "n2:co22","n2:heat2","n2:precip2","co22:heat2","co22:precip2","heat2:precip2",
                              "n2:co22:heat2","n2:co22:precip2","n2:heat2:precip2","co22:heat2:precip2",
                              "n2:co22:heat2:precip2",
                              "perc.n","c.n","ph","water.content",
                              "npp.std","plant.rich","NMDS.Plants.1","NMDS.Plants.2")
levels_display_threemodels <- c("N","C","T","P","F03","B11",
                                "NC","NT","CT","NP","CP","TP",
                                "NCT","NCP","NTP","CTP",
                                "NCTP", 
                                "%n","c/n","ph","moist",
                                "npp","p.rich","p.nmds1","p.nmds2")

coefs_main$model <- "main"
coefs_interact$model <- "interact"
coefs_plantsoil$model <- "plantsoil"

# Plot observed richness
coefs_main %>%
  full_join(coefs_interact) %>%
  full_join(coefs_plantsoil) %>%
  mutate(model = factor(model, levels=c("main","interact","plantsoil"))) %>%
  mutate(treatment = factor(treatment, levels=levels_order_threemodels),
         variable = factor(variable, levels=c("obs","sha")),
         community = factor(community, levels=c("ITS","BAC"))) %>%
  plot_lme_coefs_richness(levels_display_threemodels, "obs")
# Save as 7"x11" PDF

# Plot Shannon diversity
coefs_main %>%
  full_join(coefs_interact) %>%
  full_join(coefs_plantsoil) %>%
  mutate(model = factor(model, levels=c("main","interact","plantsoil"))) %>%
  mutate(treatment = factor(treatment, levels=levels_order_threemodels),
         variable = factor(variable, levels=c("obs","sha")),
         community = factor(community, levels=c("ITS","BAC"))) %>%
  plot_lme_coefs_richness(levels_display_threemodels, "sha")
# Save as 7"x11" PDF

par(mfrow=c(1,1))
