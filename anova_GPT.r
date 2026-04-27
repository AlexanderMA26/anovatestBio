# -----------------------------
# Initialize Data
# -----------------------------

# Treatment groups (pH levels)
Treatment <- factor(rep(c("pH2", "pH4", "pH7", "pH10", "pH12"), each = 5))

# Trial distances (from your dataset)
Response <- c(
  # pH 2
  94.43, 90.48, 49.44, 51.09, 67.68,
  
  # pH 4
  379.34, 363.01, 402.47, 366.24, 386.76,
  
  # pH 7
  279.94, 258.16, 283.39, 308.44, 279.89,
  
  # pH 10
  297.51, 271.25, 275.92, 266.48, 302.37,
  
  # pH 12
  213.36, 194.74, 218.64, 221.46, 238.13
)

# Combine into dataframe
data <- data.frame(Treatment, Response)

# -----------------------------
# Run ANOVA
# -----------------------------
anova_model <- aov(Response ~ Treatment, data = data)

# Show ANOVA table
summary(anova_model)

# -----------------------------
# Post-hoc Test (Tukey)
# -----------------------------
TukeyHSD(anova_model)

# -----------------------------
# Assumption Checks
# -----------------------------

# Normality of residuals
shapiro.test(residuals(anova_model))

# Homogeneity of variance
bartlett.test(Response ~ Treatment, data = data)

# Optional: Levene's test (more robust)
if (!require(car)) install.packages("car")
library(car)
leveneTest(Response ~ Treatment, data = data)

# -----------------------------
# Optional Visualization
# -----------------------------
boxplot(Response ~ Treatment,
        data = data,
        main = "Distance by Treatment Group",
        xlab = "Treatment (pH)",
        ylab = "Distance",
        border = "black")