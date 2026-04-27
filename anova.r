# 1. Create the dataset
# We will use the 'Closest Reference' (Glucose mg/dL) as our dependent variable
# based on the matching results from the previous step.

ph_levels <- factor(c(rep("pH_2", 5), rep("pH_4", 5), rep("pH_7", 5), rep("pH_10", 5), rep("pH_12", 5)),
                    levels = c("pH_2", "pH_4", "pH_7", "pH_10", "pH_12"))

# These are the glucose values (mg/dL) identified by our distance calculations
glucose_values <- c(300, 300, 0, 0, 0,        # pH 2
                    100, 300, 0, 300, 0,      # pH 4
                    300, 1000, 300, 300, 300, # pH 7
                    0, 100, 100, 300, 0,      # pH 10
                    100, 100, 0, 0, 0)        # pH 12

lab_data <- data.frame(pH = ph_levels, Glucose = glucose_values)

# 2. Perform the One-Way ANOVA
anova_results <- aov(Glucose ~ pH, data = lab_data)

# 3. View the summary table
# Look for the 'Pr(>F)' value. If it is < 0.05, your results are significant.
summary(anova_results)

# 4. Perform Tukey's HSD Post-Hoc Test
# This tells you exactly WHICH groups are significantly different from each other.
tukey_test <- TukeyHSD(anova_results)
print(tukey_test)

# 5. Visualize the data with a Boxplot
boxplot(Glucose ~ pH, data = lab_data,
        main = "Lactase Activity Across pH Levels",
        xlab = "pH Treatment Group",
        ylab = "Glucose Concentration (mg/dL)",
        col = "lightblue",
        border = "darkblue")