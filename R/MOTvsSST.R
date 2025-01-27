# segmented regression of MOT on SST 
library(segmented)

csvpath <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/"
csvfile <- "MOTvsSST_v05.csv"

ocean_tmp <- read.csv(paste(csvpath, csvfile, sep = ""))

# sort data by SST_norm values
ocean_tmp <- ocean_tmp[order(ocean_tmp$SST_norm), ]
head(ocean_tmp)

n <- dim(ocean_tmp)[1]
n
ocean_tmp$obs_num <- seq(1, n, by = 1)

attach(ocean_tmp)

plot(MOT_norm ~ SST_norm, data = ocean_tmp)
ylim=range(ocean_tmp$MOT_norm, ocean_tmp$SST_norm, na.rm=T, lwd=2)
ylim
ylim = c(-6.0, 20.0)
ylim_resid = c(-5, 5)

# new data for predictions
newdata <- data.frame(ocean_tmp$SST_norm)
names(newdata) <- "SST_norm"

# Model 1: basic linear trend =================================================

lm1 <- lm(MOT_norm ~ SST_norm, data = ocean_tmp)
AIC_lm1 <- round(AIC(lm1), 2)
summary(lm1)
AIC_lm1

rss <- sum(lm1$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
np <- 2
rsq_lm1 <- 1 - (rss/tss)
adj_rsq_lm1 <- 1 - (rss / (n - np)) / (tss / (n - 1))
rsq_lm1 <- round(rsq_lm1, digits = 3)
adj_rsq_lm1 <- round(adj_rsq_lm1, digits = 3)
rsq_lm1; adj_rsq_lm1

hist(lm1$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm1$residuals); qqline(lm1$residuals)
shapiro.test(lm1$residuals)

lm1_text <- paste("Linear regression: AIC = ", as.character(AIC_lm1),  ", R² = ", as.character(rsq_lm1),
                  ", adj R² = ", as.character(adj_rsq_lm1), sep = "")

# get and plot prediction intervals and confidence intervals
pred_int <- predict(lm1, newdata, int="p")
conf_int <- predict(lm1, newdata, int="c")

# linear_v05
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm1_text)
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")

# linear_residual_plot_v05
plot(lm1$residuals ~ lm1$fitted, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Linear Regression", xlab =  "Linear Regression Fitted Values", sub = lm1_text)
lo_lm1_resid <- loess(lm1$residuals ~ lm1$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lm1_resid$fitted ~ lm1$fitted, col = "red", lwd = 2)

linear_residuals_out <- data.frame(cbind(lm1$fitted, lm1$residuals, lm1$fitted, lo_lm1_resid$fitted))
names(linear_residuals_out) <- c("linear_fitted", "linear_residuals", "smoothed_x", "smoothed_y")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "linear_residual_plot_v05.csv"
write.csv(linear_residuals_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# residual diagnostic plots
# Cook's distance
lm1_cook <- cooks.distance(lm1)
plot(lm1, which = c(4))

# linear_residual_diagnostics_v05
linear_residual_diagnostics_out <- data.frame(cbind(SST_norm, MOT_norm, lm1$fitted, lm1$residuals, lm1_cook, obs_num, Model))
names(linear_residual_diagnostics_out) <- c("SST_norm", "MOT_norm, linear_fitted", "linear_fitted", 
                                            "linear_residuals", "cooks_distance", "obs_num", "Model")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "linear_residual_diagnostics_v05.csv"
write.csv(linear_residual_diagnostics_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# Model 2: lowess =============================================================

span = 0.8; deg = 2
lo1 <- loess(MOT_norm ~ SST_norm, span = span, degree = deg, iterations = 3)
lo1

enp <- ceiling(lo1$enp)
sigma2 <- lo1$s
max_log_lik = -1.0 * (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
max_log_lik = -1.0 * (n / 2) * log(2 * pi * sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
AIC_lo1 <- - 2.0 * max_log_lik + 2.0 * enp
AIC_lo1 <- round(AIC_lo1, 2)
AIC_lo1

rss <- sum(lo1$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_lo1 <- 1 - (rss/tss)
adj_rsq_lo1 <- 1 - (rss / (n - enp)) / (tss / (n - 1))
rsq_lo1 <- round(rsq_lo1, digits = 3)
adj_rsq_lo1 <- round(adj_rsq_lo1, digits = 3)
rsq_lo1; adj_rsq_lo1

hist(lo1$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lo1$residuals); qqline(lo1$residuals)
shapiro.test(lo1$residuals)

lo1_text <- paste("Loess: AIC = ", as.character(AIC_lo1), ", R² = ", as.character(rsq_lo1),
                  ", adj R² = ", as.character(adj_rsq_lo1), sep = "")

# get and plot prediction intervals and confidence intervals
lo1_pred <- predict(lo1, se = TRUE)

# loess_v05
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lo1_text)
lines(lo1$fitted ~ SST_norm, lwd = 2, col = "red")
lines((lo1_pred$fit + qt(0.975, lo1_pred$df) * lo1_pred$se) ~ SST_norm, lty=2, col = "red")
lines((lo1_pred$fit - qt(0.975, lo1_pred$df) * lo1_pred$se) ~ SST_norm, lty=2, col = "red")
lines((lo1_pred$fit + 2.0 * sd(lo1$residuals)) ~ SST_norm, lty=2, col = "black")
lines((lo1_pred$fit - 2.0 * sd(lo1$residuals)) ~ SST_norm, lty=2, col = "black")

# # residual plot
# plot(lo1$residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lo1_text)
# lo_lo1_resid <- loess(lo1$residuals ~ SST_norm, span=0.90, degree=1, iterations = 3)
# lines(lo_lo1_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# loess_residual_plot_v05
plot(lo1$residuals ~ lo1$fitted, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Loess Local Regression", xlab =  "Loess Local Regression Fitted Values", sub = lo1_text)
lo_lo1_resid <- loess(lo1$residuals ~ lo1$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lo1_resid$fitted ~ lo1$fitted, col = "red", lwd = 2)

loess_residuals_out <- data.frame(cbind(lo1$fitted, lo1$residuals, lo1$fitted, lo_lo1_resid$fitted))
names(loess_residuals_out) <- c("loess_fitted", "loess_residuals", "smoothed_x", "smoothed_y")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "loess_residual_plot_v05.csv"
write.csv(loess_residuals_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# Model 3: Third-order polynomial trend =======================================

lm_poly1 <- lm(MOT_norm ~ SST_norm + I(SST_norm^2) + I(SST_norm^3), data = ocean_tmp)
AIC_lm_poly1 <- round(AIC(lm_poly1), 2)
summary(lm_poly1)
AIC_lm_poly1

rss <- sum(lm_poly1$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_lm_poly1 <- 1 - (rss/tss)
adj_rsq_lm_poly1 <- 1 - (rss / (n - enp)) / (tss / (n - 1))
rsq_lm_poly1 <- round(rsq_lo1, digits = 3)
adj_rsq_lm_poly1 <- round(adj_rsq_lm_poly1, digits = 3)
rsq_lm_poly1; adj_rsq_lm_poly1

hist(lm_poly1$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm_poly1$residuals); qqline(lm_poly1$residuals)
shapiro.test(lm_poly1$residuals)

lm_poly1_text <- paste("Third-Order Polynomial Trend: AIC = ", as.character(AIC_lm_poly1), ", R² = ", as.character(rsq_lo1),
                       ", adj R² = ", as.character(adj_rsq_lo1), sep = "")

# get and plot prediction intervals and confidence intervals
pred_int <- predict(lm_poly1, newdata, int="p")
conf_int <- predict(lm_poly1, newdata, int="c")
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim=ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm_poly1_text)
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")

# residual plot
plot(lm_poly1$residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lm_poly1_text)
lo_lm_poly1_resid <- loess(lo1$residuals ~ SST_norm, span=0.90, degree=1, iterations = 3)
lines(lo_lm_poly1_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# Model 4: Segmented linear trend =============================================

lm1_seg <- segmented(lm1, seg.Z = ~ SST_norm, psi=list(SST_norm = c(1.0, 5.0)), 
                     control=seg.control(display=T, h=0.2) )
summary(lm1_seg)
AIC_lm1_seg <- AIC(lm1_seg)
AIC_lm1_seg <- round(AIC_lm1_seg, 2)
slope(lm1_seg)
AIC_lm1_seg

rss <- sum(lm1_seg$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_lm1_seg <- 1 - (rss/tss)
adj_rsq_lm1_seg <- 1 - (rss / (n - 6)) / (tss / (n - 1))
rsq_lm1_seg <- round(rsq_lm1_seg, digits = 3)
adj_rsq_lm1_seg <- round(adj_rsq_lm1_seg, digits = 3)
rsq_lm1_seg; adj_rsq_lm1_seg

hist(lm1_seg$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm1_seg$residuals); qqline(lm1_seg$residuals)
shapiro.test(lm1_seg$residuals)

lm1_seg_text <- paste("Segmented regression: AIC = ", as.character(AIC_lm1_seg), ", R² = ", as.character(rsq_lm1_seg),
                      ", adj R² = ", as.character(adj_rsq_lm1_seg), sep = "")
# segmented_v05
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim=ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm1_seg_text)
lines(lm1_seg$fitted.values ~ SST_norm, col="purple", lwd=2)
lines(lm1_seg, col="purple", lwd=2)

pred_int <- predict(lm1_seg, newdata, int="p")
conf_int <- predict(lm1_seg, newdata, int="c")
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")

# # residual plot
# plot(lm1_seg$residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lm1_seg_text)
# lo_lm1_seg_resid <- loess(lm1_seg$residuals ~ SST_norm, span=0.90, degree = 1, iterations = 5)
# lines(lo_lm1_seg_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# segmented_residual_plot
plot(lm1_seg$residuals ~ lm1_seg$fitted, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Segmented Regression", xlab =  "Segmented Regression Fitted Values", sub = lm1_seg_text)
lo_lm1_seg_resid <- loess(lm1_seg$residuals ~ lm1_seg$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lm1_seg_resid$fitted ~ lm1_seg$fitted, col = "red", lwd = 2)

segmented_residuals_out <- data.frame(cbind(lm1_seg$fitted, lm1_seg$residuals, lm1_seg$fitted, lo_lm1_seg_resid$fitted))
names(segmented_residuals_out) <- c("segmented_fitted", "segmented_residuals", "smoothed_x", "smoothed_y")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "segmented_residual_plot_v05.csv"
write.csv(segmented_residuals_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "segmented_residual_plot_v05.csv"
write.csv(seqmaented_residuals_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# linear_residual_diagnostics_v05
# Cook's distance
lm1_seg_cook <- cooks.distance(lm1_seg)
plot(lm1_seg, which = c(4))
segmented_residual_diagnostics_out <- data.frame(cbind(SST_norm, MOT_norm, lm1_seg$fitted, lm1_seg$residuals, lm1_seg_cook, obs_num, Model))
names(segmented_residual_diagnostics_out) <- c("SST_norm", "MOT_norm, segmented_fitted", "segmented_fitted", 
                                            "segmented_residuals", "cooks_distance", "obs_num", "Model")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "segmented_residual_diagnostics_v05.csv"
write.csv(segmented_residual_diagnostics_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# Model 5: Segmented linear trend -- one breakpoint ===========================

lm2_seg <- segmented(lm1, seg.Z = ~ SST_norm, psi=list(SST_norm = c(5.0)), 
                     control=seg.control(display=T, h=0.2) )
summary(lm2_seg)
AIC_lm2_seg <- AIC(lm2_seg)
AIC_lm2_seg <- round(AIC_lm2_seg, 2)
slope(lm2_seg)
AIC_lm2_seg

rss <- sum(lm2_seg$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_lm2_seg <- 1 - (rss/tss)
adj_rsq_lm2_seg <- 1 - (rss / (n - 6)) / (tss / (n - 1))
rsq_lm2_seg <- round(rsq_lm2_seg, digits = 3)
adj_rsq_lm2_seg <- round(adj_rsq_lm2_seg, digits = 3)
rsq_lm2_seg; adj_rsq_lm2_seg

hist(lm2_seg$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm2_seg$residuals); qqline(lm2_seg$residuals)
shapiro.test(lm2_seg$residuals)

lm2_seg1_text <- paste("Segmented regression: AIC = ", as.character(AIC_lm2_seg), ", R² = ", as.character(rsq_lm2_seg),
                      ", adj R² = ", as.character(adj_rsq_lm2_seg), sep = "")

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim=ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm2_seg1_text)
lines(lm2_seg$fitted.values ~ SST_norm, col="purple", lwd=2)
lines(lm2_seg, col="purple", lwd=2)

pred_int <- predict(lm2_seg, newdata, int="p")
conf_int <- predict(lm2_seg, newdata, int="c")
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")

# residual plot
plot(lm2_seg$residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lm2_seg1_text)
lo_lm2_seg_resid <- loess(lm2_seg$residuals ~ SST_norm, span=0.90, degree = 1, iterations = 5)
lines(lo_lm2_seg_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# Model 6: Segmented linear trend -- fixed breakpoints ========================

# generate "design matrix" variables
bp1 <- 1.0; bp2 <- 5.290
n
x1 <- SST_norm
x2 <- rep(0, n)
x2[SST_norm >= bp1 & SST_norm <= bp2 + 0.001] <- 1.0
x3 <- rep(0, n)
x3[SST_norm >= bp1 & SST_norm <= bp2 + 0.001] <- SST_norm[SST_norm >= bp1 & SST_norm < bp2 + 0.001]
x4 <- rep(0, n)
x4[SST_norm > bp2 + 0.001] <- 1.0
x5 <- rep(0, n)
x5[ SST_norm > bp2 + 0.001] <- SST_norm[SST_norm > bp2 + 0.001]

lm3_seg <- lm(MOT_norm ~ x1 + x2 + x3 + x4 + x5)
summary(lm3_seg)
AIC_lm3_seg <- AIC(lm3_seg)
AIC_lm3_seg <- round(AIC_lm3_seg, 2)
# slope(lm3_seg)
AIC_lm3_seg

rss <- sum(lm3_seg$residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_lm3_seg <- 1 - (rss/tss)
adj_rsq_lm3_seg <- 1 - (rss / (n - 6)) / (tss / (n - 1))
rsq_lm3_seg <- round(rsq_lm3_seg, digits = 3)
adj_rsq_lm3_seg <- round(adj_rsq_lm3_seg, digits = 3)
rsq_lm3_seg; adj_rsq_lm3_seg

hist(lm3_seg$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm3_seg$residuals); qqline(lm3_seg$residuals)
shapiro.test(lm3_seg$residuals)

# new data for fitted values
pred_range <- range(ocean_tmp$SST_norm)
pred_range <- round(pred_range, 2)
x1 <- seq(pred_range[1], pred_range[2], by =  0.01) # 0.01) # 
nd <- length(x1)
x2 <- rep(0, nd)
x2[x1 >= bp1 & x1 <= bp2 + 0.001] <- 1.0
x3 <- rep(0, nd)
x3[x1 >= bp1 & x1 <= bp2 + 0.001] <- x1[x1 >= bp1 & x1 <= bp2 + 0.001]
x4 <- rep(0, nd)
x4[x1 >= bp2 + 0.001] <- 1.0
x5 <- rep(0, nd)
x5[x1 >= bp2 + 0.001] <- x1[x1 >= bp2 + 0.001]
newdata <- data.frame(cbind(x1, x2, x3, x4, x5))
names(newdata) <- c("x1", "x2", "x3", "x4", "x5")

pred_int <- predict(lm3_seg, newdata, int="p")
conf_int <- predict(lm3_seg, newdata, int="c")

lm3_seg_text <- paste("Segmented regression: AIC = ", as.character(AIC_lm3_seg), ", R² = ", as.character(rsq_lm3_seg),
                      ", adj R² = ", as.character(adj_rsq_lm3_seg), sep = "")

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm3_seg_text)
abline(a = 0.0, b = 1.0, col = "gray")
lines(conf_int[, 1] ~ x1, lwd = 2, col = "red")
matlines(x1, pred_int, lty=c(1,2,2), col="black")
matlines(x1, conf_int, lty=c(1,2,2), col="red")

# residual plot
plot(lm3_seg$residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lm3_seg_text)
lo_lm3_seg_resid <- loess(lm3_seg$residuals ~ SST_norm, span=0.90, degree = 1, iterations = 5)
lines(lo_lm3_seg_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# Model 7: Smoothing spline ===================================================

spar = 0.95
ss1 <- smooth.spline(MOT_norm ~ SST_norm, spar = spar, keep.data = TRUE)
ss1
ss1_fit <- predict(ss1, SST_norm, deriv = 0)$y
ss1_residuals <- MOT_norm - ss1_fit
ss1_slope <- predict(ss1, SST_norm, deriv = 1)$y

df <- ss1$df # ceiling(ss1$df)
sigma2 <- sum(ss1_residuals^2)/(n - 1.0)
rss <- sum(ss1_residuals^2 )
# http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
max_log_lik <- (-1.0 / (2.0 * sigma2)) * rss - 1.0 * (n / 2.0) * log(sigma2) - 1.0 * (n / 2.0) * log(2 * pi)
AIC_ss1 <- - 2.0 * max_log_lik + 2.0 * df
AIC_ss1 <- round(AIC_ss1, 2)
AIC_ss1

rss <- sum(ss1_residuals^2 )
tss <- sum((MOT_norm - mean(MOT_norm))^2)
rsq_ss1 <- 1 - (rss/tss)
adj_rsq_ss1 <- 1 - (rss / (n - ss1$df)) / (tss / (n - 1))
rsq_ss1 <- round(rsq_ss1, digits = 3)
adj_rsq_ss1 <- round(adj_rsq_ss1, digits = 3)
rsq_ss1; adj_rsq_ss1

hist(ss1_residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(ss1_residuals); qqline(lo1$residuals)
shapiro.test(ss1_residuals)

ss1_text <- paste("Smoothing spline: AIC = ", as.character(AIC_ss1), ", R² = ", as.character(rsq_ss1),
                  ", adj R² = ", as.character(adj_rsq_ss1), sep = "")

# get and plot prediction intervals and confidence intervals
res <- (ss1$y - ss1$yin)/(1 - ss1$lev)
sig <- sqrt(var(res))
ci_upper <- ss1$y + 2.0 * sig * sqrt(ss1$lev)
ci_lower <- ss1$y - 2.0 * sig * sqrt(ss1$lev)

# smooth_spline_v05
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = ss1_text)
lines(ss1_fit ~ SST_norm, lwd = 2, col = "red")
# lines(ss1$y ~ ss1$x, lwd = 2, col = "blue")
lines(ci_upper ~ ss1$x, lty=2, col = "red")
lines(ci_lower ~ ss1$x, lty=2, col = "red")
lines(ss1_fit + 2.0 * sig ~ SST_norm, lty=2, col = "black")
lines(ss1_fit - 2.0 * sig ~ SST_norm, lty=2, col = "black")

# # residual plot
# plot(ss1_residuals ~ SST_norm, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, sub = lo1_text)
# lo_ss1_resid <- loess(ss1_residuals ~ SST_norm, span=0.90, degree=1, iterations = 3)
# lines(lo_ss1_resid$fitted ~ SST_norm, col = "red", lwd = 2)

# smooth_spline_residual_plot_v05
plot(ss1_residuals ~ ss1_fit, data = ocean_tmp, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values", sub = ss1_text)
lo_ss1_resid <- loess(ss1_residuals ~ ss1_fit, span=0.90, degree = 1, iterations = 3)
lines(lo_ss1_resid$fitted ~ ss1_fit, col = "red", lwd = 2)

smooth_spline_residuals_out <- data.frame(cbind(ss1_fit, ss1_residuals, ss1_fit, lo_ss1_resid$fitted))
names(smooth_spline_residuals_out) <- c("smoothing_spline_fitted", "smoothing_spline_residuals", "smoothed_x", "smoothed_y")

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "smooth_spline_residual_plot_v05.csv"
write.csv(smooth_spline_residuals_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# slope plot
plot(ss1_slope ~ SST_norm, ylim = c(0.0, 2.0), xlim = ylim, pch = 16, cex = 0.8, type = "o", col = "red", sub = ss1_text)

# linear_loess_segmented_smooth_spline_v05
oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8,
     sub = "Linear, loess, segmented-regression and smoothing spline fits")
# plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = c(-4, 6), xlim = c(-4, 6), pch = 16, cex = 0.8,
#      sub = "Linear, loess, segmented-regression and smoothing spline fits")
abline(a = 0.0, b = 1.0, col = "gray")
lines(lm1$fitted ~ SST_norm, lwd = 2, col = "red")
lines(lo1$fitted ~ SST_norm, lwd = 2, col = "blue")
lines(lm1_seg$fitted.values ~ SST_norm, col="purple", lwd=2)
# lines(conf_int[, 1] ~ x1, lwd = 2, col = "purple")
lines(ss1_fit ~ SST_norm, lwd = 2, col = "magenta")

legend("bottomright", c(lm1_text, lo1_text, lm1_seg_text, ss1_text), col = c("red", "blue", "purple", "magenta"), lwd = 2, cex = 0.85)

# predictions =================================================================

pred_range <- range(ocean_tmp$SST_norm)
pred_range <- round(pred_range, 2)
x <- seq(pred_range[1], pred_range[2], by =  0.01) # 0.01) # 
newdata <- data.frame(x)
names(newdata) <- "SST_norm"
npred <- length(x)

# linear
# get and plot prediction intervals and confidence intervals
pred_int <- predict(lm1, newdata, int="p")
conf_int <- predict(lm1, newdata, int="c")

linear_out <- data.frame(cbind(newdata, pred_int, conf_int[, 2:3]))
names(linear_out) <- c("SST_norm", "linear_fit", "linear_pred_lwr", "linear_pred_upr", "linear_conf_lwr", "linear_conf_upr")
head(linear_out)

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm1_text)
abline(a = 0.0, b = 1.0, col = "gray")
abline(lm1, col = "black", lwd = 2)
points(linear_out$linear_fit ~ linear_out$SST_norm, pch = 16, cex = 0.6, col = "red")

# loess
# get and plot prediction intervals and confidence intervals
lo1_pred <- predict(lo1, se = TRUE, newdata = newdata)
conf_lwr <- lo1_pred$fit + qt(0.975, lo1_pred$df) * lo1_pred$se
conf_upr <- lo1_pred$fit - qt(0.975, lo1_pred$df) * lo1_pred$se
pred_lwr <- lo1_pred$fit + 2.0 * sd(lo1$residuals)
pred_upr <- lo1_pred$fit - 2.0 * sd(lo1$residuals)

loess_out <- data.frame(cbind(newdata, lo1_pred$fit, pred_lwr, pred_upr, conf_lwr, conf_upr))
names(loess_out) <- c("SST_norm", "loess_fit", "loess_pred_lwr", "loess_pred_upr", "loess_conf_lwr", "loess_conf_upr")
head(loess_out)

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim = ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lo1_text)
abline(a = 0.0, b = 1.0, col = "gray")
lines(lo1$fitted ~ SST_norm, lwd = 2, col = "black")
points(loess_out$loess_fit ~ loess_out$SST_norm, pch = 16, cex = 0.6, col = "red")

# slope plot
lo1_slope <- diff(lo1_pred$fit)/diff(x)
plot(lo1_slope ~ x[2:npred], ylim = c(0.0, 2.0), xlim = ylim, xlab = "SST_norm", ylab = "Loess Slope", pch = 16, cex = 0.8, type = "o", col = "red", sub = lo1_text)

# segmented
# get and plot prediction intervals and confidence intervals
pred_int <- predict(lm1_seg, newdata, int="p")
conf_int <- predict(lm1_seg, newdata, int="c")

segmented_out <- data.frame(cbind(newdata, pred_int, conf_int[, 2:3]))
names(segmented_out) <- c("SST_norm", "segmented_fit", "segmented_pred_lwr", "segmented_pred_upr", "segmented_conf_lwr", "segmented_conf_upr")
head(segmented_out)

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim=ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm1_seg_text)
abline(a = 0.0, b = 1.0, col = "gray")
lines(lm1_seg$fitted.values ~ SST_norm, col="black", lwd=2)
points(segmented_out$segmented_fit ~ segmented_out$SST_norm, pch = 16, cex = 0.6, col = "red")

# segmented -- fixed breakpoints
# get and plot prediction intervals and confidence intervals
# new data for fitted values
pred_range <- range(ocean_tmp$SST_norm)
pred_range <- round(pred_range, 2)
x1 <- seq(pred_range[1], pred_range[2], by =  0.01) # 0.01) # 
nd <- length(x1)
x2 <- rep(0, nd)
x2[x1 >= bp1 & x1 <= bp2 + 0.001] <- 1.0
x3 <- rep(0, nd)
x3[x1 >= bp1 & x1 <= bp2 + 0.001] <- x1[x1 >= bp1 & x1 <= bp2 + 0.001]
x4 <- rep(0, nd)
x4[x1 >= bp2 + 0.001] <- 1.0
x5 <- rep(0, nd)
x5[x1 >= bp2 + 0.001] <- x1[x1 >= bp2 + 0.001]
newdata <- data.frame(cbind(x1, x2, x3, x4, x5))
names(newdata) <- c("x1", "x2", "x3", "x4", "x5")

pred_int <- predict(lm3_seg, newdata, int="p")
conf_int <- predict(lm3_seg, newdata, int="c")

segmented_fixed_bp_out <- data.frame(cbind(newdata, pred_int, conf_int[, 2:3]))
names(segmented_fixed_bp_out) <- c("x1", "x2", "x3", "x4", "x5", "segmented_fit", "segmented_pred_lwr", "segmented_pred_upr", "segmented_conf_lwr", "segmented_conf_upr")
head(segmented_fixed_bp_out)

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, data = ocean_tmp, ylim=ylim, xlim = ylim, pch = 16, cex = 0.8, sub = lm3_seg_text)
abline(a = 0.0, b = 1.0, col = "gray")
lines(lm3_seg$fitted.values ~ SST_norm, col="black", lwd=2)
points(conf_int[, 1] ~ x1, pch = 16, cex = 0.6, col = "red")

# smoothing spline
# get and plot prediction intervals and confidence intervals

# get smooth leverages (hat-matrix diagonal values) (using (what else?) a smoothing spline)
lev_smooth <- smooth.spline(ss1$lev ~ ss1$x, spar = 0.675, keep.data = TRUE)
plot(lev_smooth$y ~ lev_smooth$x, type = "l", col = "black")
points(lev_smooth$yin ~ lev_smooth$x, pch = 16, cex = 0.6, col = "blue")
lev_predict <- predict(lev_smooth, x, deriv = 0)$y
points(lev_predict ~ x, pch = 16, cex = 0.6, col = "red")

# other CI components
ss1_pred <- predict(ss1, x, deriv = 0)$y
res <- (ss1$y - ss1$yin)/(1 - ss1$lev)
sig <- sqrt(var(res))

conf_upr <- ss1_pred + 2.0 * sig * sqrt(lev_predict)
conf_lwr <- ss1_pred - 2.0 * sig * sqrt(lev_predict)
pred_lwr <- ss1_pred + 2.0 * sig
pred_upr <- ss1_pred - 2.0 * sig

oldpar <- par(pty="s")
plot(ss1_pred ~ x, ylim = ylim, xlim = ylim, type = "l", lwd = 2, sub = ss1_text)
abline(a = 0.0, b = 1.0, col = "gray")
lines(conf_upr ~ x, lty=2, col = "red")
lines(conf_lwr ~ x, lty=2, col = "red")
lines(pred_upr ~ x, lty=2, col = "black")
lines(pred_lwr ~ x, lty=2, col = "black")

# slope plot
# ss1_pred_slope <- diff(ss1_pred)/diff(x)
ss1_pred_slope <- predict(ss1, x, deriv = 1)$y
max(ss1_pred_slope); min(ss1_pred_slope)
# plot(ss1_pred_slope ~ x[2:npred], ylim = c(0.0, 2.0), type = "l", xlim = ylim, xlab = "SST_norm", 
#      ylab = "Smoothing Spline Slope", pch = 16, sub = ss1_text)
plot(ss1_pred_slope ~ x, ylim = c(0.0, 2.0), type = "l", xlim = ylim, xlab = "SST_norm", 
     ylab = "Smoothing Spline Slope", pch = 16, sub = ss1_text)
points(ss1_slope ~ SST_norm, pch = 16, cex = 0.6, col = "red")

ss_out <- data.frame(cbind(x, ss1_pred, pred_lwr, pred_upr, conf_lwr, conf_upr, ss1_pred_slope))
names(ss_out) <- c("SST_norm", "ss_fit", "ss_pred_lwr", "ss_pred_upr", "ss_conf_lwr", "ss_conf_upr", "ss_slope")
head(ss_out)

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "linear_predictions_v05.csv"
write.csv(linear_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)
out_csv <- "loess_predictions_v05.csv"
write.csv(loess_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)
out_csv <- "segmented_predictions_v05.csv"
write.csv(segmented_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)
out_csv <- "segmented_predictions_fixed_bp_v05.csv"
write.csv(segmented_fixed_bp_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)
out_csv <- "smooth_spine_predictions_v05.csv"
write.csv(ss_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

# estimate MOT_norm and get slopes

csvpath <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/"
csvfile <- "dSST_dGMST.csv"
tmps <- read.csv(paste(csvpath, csvfile, sep = ""))
head(tmps)
summary(tmps)

tmps$dMOT_est <- predict(ss1, tmps$dSST, deriv = 0)$y
head(tmps$dMOT_est); tail(tmps$dMOT_est)
plot(tmps$dMOT_est ~ tmps$dSST, pch = 16, cex = 0.6, col = "blue")

tmps$HSE <- predict(ss1, tmps$dSST, deriv = 1)$y
summary(tmps)
oldpar <- par(pty="s")
plot(tmps$dSST ~ tmps$Age_Ma, type = "o", pch = 16, cex = 0.5, col = "lightblue", ylim = c(-4, 4), xlim = c(4.5, 0),
     xlab = "Age Ma", ylab = "dSST and dMOT (estimated)")
lines(tmps$dMOT_est ~ tmps$Age_Ma, type = "o", pch = 16, cex = 0.5, col = "blue")
legend("bottomleft", c("dSST", "dMOT (estimated)"), col = c("lightblue", "blue"), lwd = 2, cex = 0.85)

plot(tmps$HSE ~ tmps$Age_Ma, type = "o", pch = 16, cex = 0.5, col = "magenta", ylim = c(0, 2), xlim = c(4.5, 0),
     xlab = "Age Ma", ylab = "HSE")

max(tmps$HSE[tmps$dSST < 4]); min(tmps$HSE[tmps$dSST < 4])
max(tmps$HSE); min(tmps$HSE)
max(ss1_pred_slope[x < 4]); min(ss1_pred_slope[x < 4])

dMOT_out <- data.frame(cbind(tmps$Age_Ma, tmps$dSST, tmps$dMOT_est, tmps$HSE))
dMOT_out <- tmps[, c(1,2,8,9)]
head(dMOT_out)

out_path <- "/Users/bartlein/Dropbox/WorkCurrent/globalT/paper2/MOTvsSST/v05/"
out_csv <- "dMOT_est_v05.csv"
write.csv(dMOT_out, paste(out_path, out_csv, sep = ""), row.names = FALSE)

oldpar <- par(pty="s")
plot(MOT_norm ~ SST_norm, ylim = ylim, xlim = ylim, type = "n")
abline(a = 0.0, b = 1.0, col = "gray")
text(MOT_norm ~ SST_norm, labels = ocean_tmp$obs, cex = 0.5)

plot(MOT_norm ~ SST_norm, ylim = c(-4, 7), xlim = c(-4, 7), type = "n")
abline(a = 0.0, b = 1.0, col = "gray")
text(MOT_norm ~ SST_norm, labels = ocean_tmp$obs, cex = 0.5)

cor(data.frame(lm1$residuals, lm1_seg$residuals, lo1$residuals, ss1_residuals))

# # loess AIC test
# lm_AIC_test <- lm(MOT_norm ~ lo1$fitted)
# summary(lm_AIC_test)
# AIC_lm_AIC_test <- round(AIC(lm_AIC_test), 2)
# AIC_lm_AIC_test
# 
# enp <- ceiling(lo1$enp) # enp <- lo1$enp
# sigma2 <- lo1$s
# max_log_lik = -1.0 * (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
# max_log_lik = -1.0 * (n / 2) * log(2 * pi * sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
# AIC_lo1 <- - 2.0 * max_log_lik + 2.0 * enp
# AIC_lo1 <- round(AIC_lo1, 2)
# AIC_lo1

# library(interp)
# lev_interp <- approx(ss1$x, ss1$lev, x)
# length(lev_interp$y)

# conf_upr <- ss1_pred + 2.0 * sig * sqrt(lev_interp$y)
# conf_lwr <- ss1_pred - 2.0 * sig * sqrt(lev_interp$y)
