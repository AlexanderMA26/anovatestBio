# ============================================================
# One-Way ANOVA: RGB Distance to Closest Reference Color
# Treatment Groups: pH 2, pH 4, pH 7, pH 10, pH 12
# Uses base R only (stats, graphics packages)
# ============================================================

# --- 1. Data entry ---
# Euclidean distance from each trial RGB to its closest reference color

data <- data.frame(
  treatment = factor(rep(c("pH 2", "pH 4", "pH 7", "pH 10", "pH 12"), each = 5),
                     levels = c("pH 2", "pH 4", "pH 7", "pH 10", "pH 12")),
  distance = c(
    # pH 2 (closest ref: 0 mg/dL = (235,222,75))
    as.numeric(dist(rbind(c(159,168,90),  c(235,222,75)))),
    as.numeric(dist(rbind(c(162,171,91),  c(235,222,75)))),
    as.numeric(dist(rbind(c(197,192,85),  c(235,222,75)))),
    as.numeric(dist(rbind(c(196,189,75),  c(235,222,75)))),
    as.numeric(dist(rbind(c(185,181,55),  c(235,222,75)))),

    # pH 4 (mixed closest refs)
    as.numeric(dist(rbind(c(185,184,91),  c(191,219,127)))),  # 100 mg/dL
    as.numeric(dist(rbind(c(164,169,70),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(207,182,80),  c(235,222,75)))),   # 0 mg/dL
    as.numeric(dist(rbind(c(170,179,80),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(191,185,77),  c(235,222,75)))),   # 0 mg/dL

    # pH 7 (closest ref: 300 mg/dL = (139,183,122), except Trial 2)
    as.numeric(dist(rbind(c(136,158,76),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(111,145,70),  c(80,159,122)))),   # 1000 mg/dL
    as.numeric(dist(rbind(c(136,154,60),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(165,172,71),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(134,155,67),  c(139,183,122)))),  # 300 mg/dL

    # pH 10 (mixed closest refs)
    as.numeric(dist(rbind(c(213,196,86),  c(235,222,75)))),   # 0 mg/dL
    as.numeric(dist(rbind(c(188,182,87),  c(191,219,127)))),  # 100 mg/dL
    as.numeric(dist(rbind(c(191,192,82),  c(191,219,127)))),  # 100 mg/dL
    as.numeric(dist(rbind(c(183,179,84),  c(139,183,122)))),  # 300 mg/dL
    as.numeric(dist(rbind(c(218,199,90),  c(235,222,75)))),   # 0 mg/dL

    # pH 12 (mixed closest refs)
    as.numeric(dist(rbind(c(212,203,109), c(191,219,127)))),  # 100 mg/dL
    as.numeric(dist(rbind(c(195,195,110), c(191,219,127)))),  # 100 mg/dL
    as.numeric(dist(rbind(c(216,208,102), c(235,222,75)))),   # 0 mg/dL
    as.numeric(dist(rbind(c(219,208,102), c(235,222,75)))),   # 0 mg/dL
    as.numeric(dist(rbind(c(217,200,90),  c(235,222,75))))    # 0 mg/dL
  )
)

# --- 2. Descriptive summary ---
cat("=== DATA SUMMARY ===\n")
print(aggregate(distance ~ treatment, data,
                function(x) round(c(mean=mean(x), sd=sd(x), min=min(x), max=max(x)), 2)))

# --- 3. Assumption checks ---

cat("\n=== ASSUMPTION 1: Bartlett's Test for Equality of Variance ===\n")
bart <- bartlett.test(distance ~ treatment, data = data)
print(bart)
if (bart$p.value < 0.05) {
  cat("  WARNING: Variances are NOT equal (p < 0.05) - use Welch ANOVA as primary test\n")
} else {
  cat("  OK: Variances are equal (p >= 0.05) - standard ANOVA assumption met\n")
}

cat("\n=== ASSUMPTION 2: Shapiro-Wilk Normality Test (per group) ===\n")
for (grp in levels(data$treatment)) {
  vals <- data$distance[data$treatment == grp]
  s <- shapiro.test(vals)
  flag <- if (s$p.value < 0.05) "NOT normal" else "Normal"
  cat(sprintf("  %-6s W = %.4f, p = %.4f  [%s]\n", grp, s$statistic, s$p.value, flag))
}

# --- 4. One-Way ANOVA ---

cat("\n=== ONE-WAY ANOVA ===\n")
model <- aov(distance ~ treatment, data = data)
anova_table <- summary(model)
print(anova_table)

# Manual eta-squared effect size
ss   <- anova_table[[1]]$"Sum Sq"
eta2 <- ss[1] / sum(ss)
cat(sprintf("\nEffect size (eta-squared): eta2 = %.4f\n", eta2))
cat("  Interpretation: 0.01 = small, 0.06 = medium, 0.14 = large\n")

# --- 5. Welch ANOVA (robust to unequal variances) ---

cat("\n=== WELCH ANOVA (recommended when variances are unequal) ===\n")
welch <- oneway.test(distance ~ treatment, data = data, var.equal = FALSE)
print(welch)

# --- 6. Post-hoc tests ---

cat("\n=== POST-HOC: Tukey HSD (follows standard ANOVA) ===\n")
print(TukeyHSD(model))

cat("\n=== KRUSKAL-WALLIS TEST (non-parametric alternative to ANOVA) ===\n")
kw <- kruskal.test(distance ~ treatment, data = data)
print(kw)

cat("\n=== POST-HOC: Pairwise Wilcoxon with Bonferroni correction (follows Kruskal-Wallis) ===\n")
print(pairwise.wilcox.test(data$distance, data$treatment, p.adjust.method = "bonferroni"))

# --- 7. Visualization ---

cols <- c("pH 2"="#E07B8A", "pH 4"="#E8A95C", "pH 7"="#7BAF6E",
          "pH 10"="#5B8DBE", "pH 12"="#9B72B0")

# Plot 1: Boxplot with jittered individual points
png("boxplot_distances.png", width = 900, height = 550, res = 120)
par(mar = c(5, 5, 4, 2), bg = "white")
boxplot(distance ~ treatment, data = data,
        col    = cols[levels(data$treatment)],
        border = "gray30",
        main   = "Distance to Closest Reference Color by pH Treatment",
        xlab   = "Treatment Group",
        ylab   = "Euclidean Distance to Closest Reference",
        cex.main = 1.1, cex.lab = 1.0, las = 1)
set.seed(42)
for (i in seq_along(levels(data$treatment))) {
  grp  <- levels(data$treatment)[i]
  vals <- data$distance[data$treatment == grp]
  jx   <- jitter(rep(i, length(vals)), amount = 0.12)
  points(jx, vals, pch = 21, bg = "white", col = "gray30", cex = 1.4, lwd = 1.5)
}
means <- tapply(data$distance, data$treatment, mean)
points(seq_along(means), means, pch = 18, col = "red", cex = 2.2)
legend("topright", legend = "Group mean", pch = 18, col = "red", cex = 0.9, bty = "n")
dev.off()

# Plot 2: Mean +/- SD bar chart
png("barplot_mean_sd.png", width = 900, height = 550, res = 120)
par(mar = c(5, 5, 4, 2), bg = "white")
grps  <- levels(data$treatment)
means <- tapply(data$distance, data$treatment, mean)[grps]
sds   <- tapply(data$distance, data$treatment, sd)[grps]
x <- barplot(means,
             col    = cols[grps],
             border = "gray30",
             ylim   = c(0, max(means + sds) * 1.25),
             main   = "Mean Distance to Closest Reference Color +/- SD",
             xlab   = "Treatment Group",
             ylab   = "Mean Euclidean Distance (+/- SD)",
             cex.main = 1.1, cex.lab = 1.0, las = 1)
arrows(x, means - sds, x, means + sds,
       angle = 90, code = 3, length = 0.08, lwd = 2, col = "gray30")
dev.off()

cat("\nPlots saved: boxplot_distances.png, barplot_mean_sd.png\n")
cat("=== ANALYSIS COMPLETE ===\n")