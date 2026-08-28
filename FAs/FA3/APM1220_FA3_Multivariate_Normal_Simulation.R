# APM1220 - Formative Assessment 3
# Custom Multivariate Normal Simulation Lab
# RStudio / R Script
# Scenario: Delivery Operations

# ============================================================
# SETUP
# ============================================================

# MASS provides mvrnorm(), which generates multivariate normal data.
library(MASS)

set.seed(1220)

# Create a folder for graphs
if (!dir.exists("figures")) {
  dir.create("figures")
}

# ============================================================
# PART A - DEFINE THE MULTIVARIATE NORMAL DISTRIBUTION
# ============================================================

# Scenario:
# We model three delivery-operation measurements:
# X1 = delivery distance (km)
# X2 = delivery time (minutes)
# X3 = fuel used (liters)

# Population model:
# X ~ N_3(mu, Sigma)

mu <- c(
  Distance_km = 12,
  Delivery_Time_min = 35,
  Fuel_Liters = 4.5
)

# Covariance matrix
Sigma <- matrix(
  c(
     9, 12, 2.7,
    12, 64, 7.2,
     2.7, 7.2, 2.25
  ),
  nrow = 3,
  byrow = TRUE
)

rownames(Sigma) <- names(mu)
colnames(Sigma) <- names(mu)

# Corresponding theoretical correlation matrix
R_theoretical <- cov2cor(Sigma)

cat("\n================ PART A ================\n")
cat("Mean vector mu:\n")
print(mu)

cat("\nCovariance matrix Sigma:\n")
print(Sigma)

cat("\nCorrelation matrix R:\n")
print(round(R_theoretical, 3))

# Check that Sigma is symmetric and positive definite.
cat("\nIs Sigma symmetric? ", isTRUE(all.equal(Sigma, t(Sigma))), "\n")

eigenvalues <- eigen(Sigma, symmetric = TRUE)$values

cat("Eigenvalues of Sigma:\n")
print(round(eigenvalues, 6))

cat("Determinant of Sigma:", det(Sigma), "\n")

# Interpretation:
# A valid covariance matrix must be symmetric and positive definite.
# The positive eigenvalues and positive determinant support validity.
#
# The diagonal entries are variances:
# Var(X1)=9, Var(X2)=64, Var(X3)=2.25.
#
# The off-diagonal entries are covariances:
# Cov(X1,X2)=12, Cov(X1,X3)=2.7, Cov(X2,X3)=7.2.
#
# The corresponding correlations are:
# Corr(X1,X2)=0.50
# Corr(X1,X3)=0.60
# Corr(X2,X3)=0.60

# ============================================================
# PART B - GENERATE AND EXPLORE THE SIMULATED DATA
# ============================================================

n <- 5000

X <- mvrnorm(
  n = n,
  mu = mu,
  Sigma = Sigma
)

X <- as.data.frame(X)

cat("\n================ PART B ================\n")
cat("Number of observations:", n, "\n")

# Simulated means, variances, covariances, correlations
sample_means <- colMeans(X)
sample_variances <- sapply(X, var)
sample_covariances <- cov(X)
sample_correlations <- cor(X)

cat("\nSimulated means:\n")
print(round(sample_means, 3))

cat("\nSimulated variances:\n")
print(round(sample_variances, 3))

cat("\nSimulated covariance matrix:\n")
print(round(sample_covariances, 3))

cat("\nSimulated correlation matrix:\n")
print(round(sample_correlations, 3))

# Direct theoretical-vs-simulated comparison
comparison <- data.frame(
  Variable = names(mu),
  Theoretical_Mean = as.numeric(mu),
  Simulated_Mean = as.numeric(sample_means),
  Theoretical_Variance = diag(Sigma),
  Simulated_Variance = as.numeric(sample_variances)
)

comparison$Mean_Difference <- comparison$Simulated_Mean -
  comparison$Theoretical_Mean

comparison$Variance_Difference <- comparison$Simulated_Variance -
  comparison$Theoretical_Variance

cat("\nTheoretical vs simulated quantities:\n")
print(round(comparison, 3))

# ---- B1: Scatterplot matrix ----
png("figures/Part_B_pairs_plot.png", width = 1000, height = 900)

pairs(
  X,
  main = "Part B: Delivery Operations Scatterplot Matrix",
  pch = 16,
  cex = 0.45
)

dev.off()

# ---- B2: Histograms for each variable ----
png("figures/Part_B_histograms.png", width = 1200, height = 400)

par(mfrow = c(1, 3))

hist(
  X$Distance_km,
  breaks = 30,
  main = "Delivery Distance",
  xlab = "Distance (km)"
)

hist(
  X$Delivery_Time_min,
  breaks = 30,
  main = "Delivery Time",
  xlab = "Time (minutes)"
)

hist(
  X$Fuel_Liters,
  breaks = 30,
  main = "Fuel Used",
  xlab = "Fuel (liters)"
)

dev.off()

# ---- B3: Correlation visualization ----
png("figures/Part_B_correlation_matrix.png", width = 800, height = 700)

par(mar = c(6, 6, 4, 2))

image(
  1:3,
  1:3,
  sample_correlations[3:1, ],
  axes = FALSE,
  xlab = "",
  ylab = "",
  main = "Part B: Simulated Correlation Matrix"
)

axis(1, at = 1:3, labels = names(mu), las = 2)
axis(2, at = 1:3, labels = rev(names(mu)), las = 2)

for (i in 1:3) {
  for (j in 1:3) {
    text(
      i,
      4 - j,
      round(sample_correlations[j, i], 2)
    )
  }
}

box()

dev.off()

# ============================================================
# PART C - EFFECT OF CORRELATION
# ============================================================

# We keep the means and standard deviations fixed.
# Only the common pairwise correlation is changed.

correlation_levels <- c(0, 0.50, 0.85)

C_data <- list()

cat("\n================ PART C ================\n")

png(
  "figures/Part_C_effect_of_correlation.png",
  width = 1500,
  height = 450
)

par(mfrow = c(1, 3))

for (rho in correlation_levels) {

  R <- matrix(rho, nrow = 3, ncol = 3)
  diag(R) <- 1

  Sigma_rho <- diag(sqrt(diag(Sigma))) %*%
    R %*%
    diag(sqrt(diag(Sigma)))

  X_rho <- mvrnorm(
    n = 1500,
    mu = mu,
    Sigma = Sigma_rho
  )

  C_data[[as.character(rho)]] <- X_rho

  plot(
    X_rho[, 1],
    X_rho[, 2],
    pch = 16,
    cex = 0.45,
    xlab = "Distance (km)",
    ylab = "Delivery Time (minutes)",
    main = paste("Common correlation =", rho)
  )

  abline(
    lm(X_rho[, 2] ~ X_rho[, 1]),
    lwd = 2
  )

  cat("\nCorrelation setting:", rho, "\n")
  cat("Sample Corr(Distance, Time):",
      round(cor(X_rho[, 1], X_rho[, 2]), 3), "\n")
}

dev.off()

# The plot should become more elongated and concentrated along
# a line as the positive correlation increases.

# ============================================================
# PART D - LINEAR COMBINATION
# ============================================================

# Define a weighted operational burden:
#
# Y = 0.5(X1) + 1(X2) + 2(X3)
#
# The weights create a simple composite score from the
# three operational measurements.

a <- c(
  0.5,
  1,
  2
)

Y <- as.numeric(as.matrix(X) %*% a)

# Theoretical result:
# a'X ~ N(a'mu, a'Sigma a)

theoretical_mean_Y <- as.numeric(t(a) %*% mu)
theoretical_variance_Y <- as.numeric(t(a) %*% Sigma %*% a)

simulated_mean_Y <- mean(Y)
simulated_variance_Y <- var(Y)

cat("\n================ PART D ================\n")
cat("Linear combination:\n")
cat("Y = 0.5(X1) + 1(X2) + 2(X3)\n")

cat("\nTheoretical mean:", theoretical_mean_Y, "\n")
cat("Simulated mean:", round(simulated_mean_Y, 3), "\n")

cat("\nTheoretical variance:", theoretical_variance_Y, "\n")
cat("Simulated variance:", round(simulated_variance_Y, 3), "\n")

# Graph of the simulated linear combination
png("figures/Part_D_linear_combination.png", width = 900, height = 650)

hist(
  Y,
  breaks = 40,
  probability = TRUE,
  main = "Part D: Weighted Operational Burden",
  xlab = "Y = 0.5(Distance) + Time + 2(Fuel)"
)

curve(
  dnorm(
    x,
    mean = theoretical_mean_Y,
    sd = sqrt(theoretical_variance_Y)
  ),
  add = TRUE,
  lwd = 2
)

abline(
  v = theoretical_mean_Y,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c("Theoretical normal density", "Theoretical mean"),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dev.off()

# ============================================================
# PART E - MAHALANOBIS DISTANCE AND CHI-SQUARE
# ============================================================

# For each observation x:
#
# D^2 = (x-mu)' Sigma^(-1) (x-mu)
#
# For a multivariate normal distribution with p = 3,
# D^2 theoretically follows chi-square with 3 degrees of freedom.

Sigma_inverse <- solve(Sigma)

X_centered <- sweep(
  as.matrix(X),
  2,
  mu,
  "-"
)

D2 <- rowSums(
  (X_centered %*% Sigma_inverse) * X_centered
)

chi_square_cutoff_90 <- qchisq(
  0.90,
  df = 3
)

inside_90 <- mean(D2 <= chi_square_cutoff_90)

cat("\n================ PART E ================\n")
cat("Degrees of freedom:", length(mu), "\n")
cat("Theoretical 90% chi-square cutoff:",
    round(chi_square_cutoff_90, 4), "\n")

cat("Simulated proportion inside the 90% region:",
    round(inside_90, 4), "\n")

cat("Theoretical proportion:", 0.90, "\n")

cat("\nMean of simulated D^2:",
    round(mean(D2), 3), "\n")

cat("Theoretical mean of Chi-square(3):", 3, "\n")

# Histogram with theoretical chi-square density
png("figures/Part_E_mahalanobis_chisquare.png", width = 950, height = 650)

hist(
  D2,
  breaks = 40,
  probability = TRUE,
  main = expression(
    "Part E: Squared Mahalanobis Distance vs " * chi[3]^2
  ),
  xlab = expression(D^2)
)

curve(
  dchisq(x, df = 3),
  add = TRUE,
  lwd = 2
)

abline(
  v = chi_square_cutoff_90,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    expression(chi[3]^2 ~ " theoretical density"),
    "90% cutoff"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dev.off()

# ============================================================
# PART F - CUSTOM INVESTIGATION: CONDITIONAL DISTRIBUTION
# ============================================================

# We investigate X2 | X1 = 15.
#
# The lecture gives:
#
# E(X2 | X1=x1)
# = mu2 + Sigma21 Sigma11^(-1) (x1-mu1)
#
# Var(X2 | X1)
# = Sigma22 - Sigma21 Sigma11^(-1) Sigma12
#
# For this 3-variable setup, the relevant values are:
# mu1 = 12
# mu2 = 35
# Sigma11 = 9
# Sigma21 = 12
# Sigma22 = 64

x1_fixed <- 15

conditional_mean <- mu["Delivery_Time_min"] +
  Sigma["Delivery_Time_min", "Distance_km"] /
  Sigma["Distance_km", "Distance_km"] *
  (x1_fixed - mu["Distance_km"])

conditional_variance <- Sigma["Delivery_Time_min", "Delivery_Time_min"] -
  Sigma["Delivery_Time_min", "Distance_km"]^2 /
  Sigma["Distance_km", "Distance_km"]

conditional_sd <- sqrt(conditional_variance)

# Approximate X1 = 15 by selecting observations in a narrow window.
window <- 0.5

conditional_sample <- X[
  abs(X$Distance_km - x1_fixed) <= window,
]

cat("\n================ PART F ================\n")
cat("Conditional distribution investigated: X2 | X1 = 15\n")

cat("\nTheoretical conditional mean:",
    round(conditional_mean, 3), "\n")

cat("Theoretical conditional variance:",
    round(conditional_variance, 3), "\n")

cat("Theoretical conditional SD:",
    round(conditional_sd, 3), "\n")

cat("\nNumber of observations used to approximate X1 = 15:",
    nrow(conditional_sample), "\n")

cat("Approximate simulated conditional mean:",
    round(mean(conditional_sample$Delivery_Time_min), 3), "\n")

cat("Approximate simulated conditional variance:",
    round(var(conditional_sample$Delivery_Time_min), 3), "\n")

# Conditional distribution graph
png("figures/Part_F_conditional_distribution.png", width = 900, height = 650)

hist(
  conditional_sample$Delivery_Time_min,
  breaks = 25,
  probability = TRUE,
  main = "Part F: Delivery Time Given Distance Near 15 km",
  xlab = "Delivery Time (minutes)"
)

curve(
  dnorm(
    x,
    mean = conditional_mean,
    sd = conditional_sd
  ),
  add = TRUE,
  lwd = 2
)

abline(
  v = conditional_mean,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = c(
    "Theoretical conditional normal density",
    "Theoretical conditional mean"
  ),
  lty = c(1, 2),
  lwd = 2,
  bty = "n"
)

dev.off()

# ============================================================
# FINAL SUMMARY
# ============================================================

cat("\n================ FINAL SUMMARY ================\n")
cat("Simulation scenario: Delivery Operations\n")
cat("Variables: Distance, Delivery Time, Fuel Used\n")
cat("Sample size:", n, "\n")

cat("\nPart D theoretical Y distribution:\n")
cat(
  "Y ~ N(",
  theoretical_mean_Y,
  ", ",
  theoretical_variance_Y,
  ")\n",
  sep = ""
)

cat("\nPart E theoretical result:\n")
cat("D^2 ~ Chi-square(3)\n")

cat("\nPart F theoretical conditional distribution:\n")
cat(
  "X2 | X1 = 15 ~ N(",
  round(conditional_mean, 3),
  ", ",
  round(conditional_variance, 3),
  ")\n",
  sep = ""
)

cat("\nAll graphs were saved in the 'figures' folder.\n")
