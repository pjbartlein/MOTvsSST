# regression analyses of dMOT vs dSST
library(segmented)

# read the delta SST and delta MOT climate-model output
csv_path <- "~/Projects/MOTvsSST/csv/"
csv_file <- "MOTvsSST.csv"

# path to saved plots
plotpath <- "~/Projects/MOTvsSST/pdf/"

ocean_temps <- read.csv(paste(csv_path, csv_file, sep = ""))
head(ocean_temps)

# sort data by dSST values
ocean_temps <- ocean_temps[order(ocean_temps$dSST), ]
head(ocean_temps)

n <- dim(ocean_temps)[1]
n
ocean_temps$obs_num <- seq(1, n, by = 1)

attach(ocean_temps)
summary(ocean_temps)

# basic scatter diagram
ylim <- c(-5.0, 20.0); xlim <- ylim
ylim_resid <- c(-5, 5)
SST_label <- expression(paste(Delta, "SST", sep = ""))
MOT_label <- expression(paste(Delta, "MOT", sep = ""))

oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = xlim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8)

# text plot with (sorted) observation labels (for identifying individual points)
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label,type = "n")
abline(a = 0.0, b = 1.0, col = "gray")
text(dMOT ~ dSST, labels = ocean_temps$obs_num, cex = 0.5)
dev.print(pdf, paste(plotpath, "obs_num.pdf", sep = ""))

plot(dMOT ~ dSST, ylim = c(-4, 7), xlim = c(-4, 7), ylab = MOT_label, xlab = SST_label,type = "n")
abline(a = 0.0, b = 1.0, col = "gray")
text(dMOT ~ dSST, labels = ocean_temps$obs, cex = 0.5)
dev.print(pdf, paste(plotpath, "obs_num_zoom.pdf", sep = ""))

# Model 1: basic linear trend =========================================================================================

lm1 <- lm(dMOT ~ dSST)
AIC_lm1 <- round(AIC(lm1), 2)
summary(lm1)
AIC_lm1

# R-sq values
rss <- sum(lm1$residuals^2 )
tss <- sum((dMOT - mean(dMOT))^2)
np <- 2
rsq_lm1 <- 1 - (rss/tss)
adj_rsq_lm1 <- 1 - (rss / (n - np)) / (tss / (n - 1))
rsq_lm1 <- round(rsq_lm1, digits = 3)
adj_rsq_lm1 <- round(adj_rsq_lm1, digits = 3)
rsq_lm1; adj_rsq_lm1

hist(lm1$residuals, breaks = seq(-4.25, 4.25, by = 0.5))
qqnorm(lm1$residuals); qqline(lm1$residuals)
shapiro.test(lm1$residuals)
plot(lm1, which = 4)

lm1_text <- paste("Linear regression: AIC = ", as.character(AIC_lm1),  ", R² = ", as.character(rsq_lm1),
                  ", adj R² = ", as.character(adj_rsq_lm1), sep = "")

# new data for prediction intervals and confidence intervals
newdata <- data.frame(ocean_temps$dSST)
names(newdata) <- "dSST"

# get prediction intervals and confidence intervals
pred_int <- predict(lm1, newdata, int="p")
conf_int <- predict(lm1, newdata, int="c")

# linear regression scatter diagram
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, sub = lm1_text)
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")
dev.print(pdf, paste(plotpath, "linear.pdf", sep = ""))

# linear residual plot
plot(lm1$residuals ~ lm1$fitted, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Linear Regression", xlab =  "Linear Regression Fitted Values", sub = lm1_text)
lo_lm1_resid <- loess(lm1$residuals ~ lm1$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lm1_resid$fitted ~ lm1$fitted, col = "red", lwd = 2)
dev.print(pdf, paste(plotpath, "linear_residual_plot.pdf", sep = ""))

# residual plot data
linear_residuals_out <- data.frame(cbind(lm1$fitted, lm1$residuals, lm1$fitted, lo_lm1_resid$fitted))
names(linear_residuals_out) <- c("linear_fitted", "linear_residuals", "smoothed_x", "smoothed_y")
out_csv <- "linear_residual_plot.csv"
write.csv(linear_residuals_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# residual diagnostic plots
# Cook's distance
lm1_cook <- cooks.distance(lm1)
plot(lm1, which = c(4))

# linear_residual_diagnostics
linear_residual_diagnostics_out <- data.frame(cbind(dSST, dMOT, lm1$fitted, lm1$residuals, lm1_cook, obs_num, Model))
names(linear_residual_diagnostics_out) <- c("dSST", "dMOT, linear_fitted", "linear_fitted", 
                                            "linear_residuals", "cooks_distance", "obs_num", "Model")
out_csv <- "linear_residual_diagnostics.csv"
write.csv(linear_residual_diagnostics_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# generate predictions over range of dSST

pred_range <- range(ocean_temps$dSST)
pred_range <- round(pred_range, 2)
x <- seq(pred_range[1], pred_range[2], by =  0.01) 
newdata <- data.frame(x)
names(newdata) <- "dSST"
npred <- length(newdata)

pred_int <- predict(lm1, newdata, int="p")
conf_int <- predict(lm1, newdata, int="c")

linear_out <- data.frame(cbind(newdata, pred_int, conf_int[, 2:3]))
names(linear_out) <- c("dSST", "linear_fit", "linear_pred_lwr", "linear_pred_upr", "linear_conf_lwr", "linear_conf_upr")
head(linear_out)
out_csv <- "linear_predictions.csv"
write.csv(linear_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# Model 2: Segmented linear trend =====================================================================================

lm1_seg <- segmented(lm1, seg.Z = ~ dSST, psi=list(dSST = c(1.0, 5.0)), 
                     control=seg.control(display=T, h=0.2) )
summary(lm1_seg)
AIC_lm1_seg <- AIC(lm1_seg)
AIC_lm1_seg <- round(AIC_lm1_seg, 2)
slope(lm1_seg)
AIC_lm1_seg

# R-sq values
rss <- sum(lm1_seg$residuals^2 )
tss <- sum((dMOT - mean(dMOT))^2)
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

# get prediction intervals and confidence intervals
pred_int <- predict(lm1_seg, newdata, int="p")
conf_int <- predict(lm1_seg, newdata, int="c")

# new data for prediction intervals and confidence intervals
newdata <- data.frame(ocean_temps$dSST)
names(newdata) <- "dSST"

# segmented regression scatter diagram
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim=ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, sub = lm1_seg_text)
lines(lm1_seg$fitted.values ~ dSST, col="purple", lwd=2)
lines(lm1_seg, col="purple", lwd=2)
abline(a = 0.0, b = 1.0, col = "gray")
matlines(newdata, pred_int, lty=c(1,2,2), col="black")
matlines(newdata, conf_int, lty=c(1,2,2), col="red")
dev.print(pdf, paste(plotpath, "segmented.pdf", sep = ""))

# segmented residual plot
plot(lm1_seg$residuals ~ lm1_seg$fitted, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Segmented Regression", xlab =  "Segmented Regression Fitted Values", sub = lm1_seg_text)
lo_lm1_seg_resid <- loess(lm1_seg$residuals ~ lm1_seg$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lm1_seg_resid$fitted ~ lm1_seg$fitted, col = "red", lwd = 2)
dev.print(pdf, paste(plotpath, "segmented_residual_plot.pdf", sep = ""))

# residual plot data
segmented_residuals_out <- data.frame(cbind(lm1_seg$fitted, lm1_seg$residuals, lm1_seg$fitted, lo_lm1_seg_resid$fitted))
names(segmented_residuals_out) <- c("segmented_fitted", "segmented_residuals", "smoothed_x", "smoothed_y")
out_csv <- "segmented_residual_plot.csv"
write.csv(segmented_residuals_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# segmented residual diagnostics
# Cook's distance
lm1_seg_cook <- cooks.distance(lm1_seg)
plot(lm1_seg_cook ~ obs_num, ylab = "Cook's distance", xlab = "Obs. number", sub = "segmented(dMOT ~ dSST)", type = "h")

# segmenteed residual diagnostics
segmented_residual_diagnostics_out <- data.frame(cbind(dSST, dMOT, lm1_seg$fitted, lm1_seg$residuals, lm1_seg_cook, obs_num, Model))
names(segmented_residual_diagnostics_out) <- c("dSST", "dMOT, segmented_fitted", "segmented_fitted", 
                                               "segmented_residuals", "cooks_distance", "obs_num", "Model")
out_csv <- "segmented_residual_diagnostics.csv"
write.csv(segmented_residual_diagnostics_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# generate predictions over range of dSST

pred_range <- range(ocean_temps$dSST)
pred_range <- round(pred_range, 2)
x <- seq(pred_range[1], pred_range[2], by =  0.01) 
newdata <- data.frame(x)
names(newdata) <- "dSST"
npred <- length(newdata)

pred_int <- predict(lm1_seg, newdata, int="p")
conf_int <- predict(lm1_seg, newdata, int="c")

segmented_out <- data.frame(cbind(newdata, pred_int, conf_int[, 2:3]))
names(segmented_out) <- c("dSST", "segmented_fit", "segmented_pred_lwr", "segmented_pred_upr", "segmented_conf_lwr", "segmented_conf_upr")
head(segmented_out)
out_csv <- "segmented_predictions.csv"
write.csv(segmented_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# Model 3: lowess =====================================================================================================

span = 0.8; deg = 2
lo1 <- loess(dMOT ~ dSST, span = span, degree = deg, iterations = 3)
lo1

enp <- ceiling(lo1$enp)
sigma2 <- lo1$s
max_log_lik = -1.0 * (n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
max_log_lik = -1.0 * (n / 2) * log(2 * pi * sigma2) - (1/(2 * sigma2)) * sum(lo1$residuals^2) 
AIC_lo1 <- - 2.0 * max_log_lik + 2.0 * enp
AIC_lo1 <- round(AIC_lo1, 2)
AIC_lo1

# R-sq values
rss <- sum(lo1$residuals^2 )
tss <- sum((dMOT - mean(dMOT))^2)
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

# loess scatter diagram
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, sub = lo1_text)
lines(lo1$fitted ~ dSST, lwd = 2, col = "red")
lines((lo1_pred$fit + qt(0.975, lo1_pred$df) * lo1_pred$se) ~ dSST, lty=2, col = "red")
lines((lo1_pred$fit - qt(0.975, lo1_pred$df) * lo1_pred$se) ~ dSST, lty=2, col = "red")
lines((lo1_pred$fit + 2.0 * sd(lo1$residuals)) ~ dSST, lty=2, col = "black")
lines((lo1_pred$fit - 2.0 * sd(lo1$residuals)) ~ dSST, lty=2, col = "black")
dev.print(pdf, paste(plotpath, "loess.pdf", sep = ""))

# loess residual plot
plot(lo1$residuals ~ lo1$fitted, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Loess Local Regression", xlab =  "Loess Local Regression Fitted Values", sub = lo1_text)
lo_lo1_resid <- loess(lo1$residuals ~ lo1$fitted, span=0.90, degree = 1, iterations = 3)
lines(lo_lo1_resid$fitted ~ lo1$fitted, col = "red", lwd = 2)
dev.print(pdf, paste(plotpath, "loess_residual_plot.pdf", sep = ""))

# residual plot data
loess_residuals_out <- data.frame(cbind(lo1$fitted, lo1$residuals, lo1$fitted, lo_lo1_resid$fitted))
names(loess_residuals_out) <- c("loess_fitted", "loess_residuals", "smoothed_x", "smoothed_y")
out_csv <- "loess_residual_plot.csv"
write.csv(loess_residuals_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# Cook's distances are not available for loess() models
# slope plot
lo1_slope <- diff(lo1_pred$fit)/diff(x)
plot(lo1_slope ~ x[2:npred], ylim = c(0.0, 2.0), xlim = ylim, xlab = "dSST", ylab = "Loess Slope", pch = 16, cex = 0.8, type = "o", col = "red", sub = lo1_text)

# generate predictions over range of dSST
pred_range <- range(ocean_temps$dSST)
pred_range <- round(pred_range, 2)
x <- seq(pred_range[1], pred_range[2], by =  0.01) 
newdata <- data.frame(x)
names(newdata) <- "dSST"
dim(newdata)[1]

# generate confidence intervals and predication intervals
lo1_pred <- predict(lo1, se = TRUE, newdata = newdata)
conf_lwr <- lo1_pred$fit + qt(0.975, lo1_pred$df) * lo1_pred$se
conf_upr <- lo1_pred$fit - qt(0.975, lo1_pred$df) * lo1_pred$se
pred_lwr <- lo1_pred$fit + 2.0 * sd(lo1$residuals)
pred_upr <- lo1_pred$fit - 2.0 * sd(lo1$residuals)

loess_out <- data.frame(cbind(newdata, lo1_pred$fit, pred_lwr, pred_upr, conf_lwr, conf_upr))
names(loess_out) <- c("dSST", "loess_fit", "loess_pred_lwr", "loess_pred_upr", "loess_conf_lwr", "loess_conf_upr")
head(loess_out)
out_csv <- "loess_predictions.csv"
write.csv(loess_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# slope plot
lo1_slope <- diff(lo1_pred$fit)/diff(x)
plot(lo1_slope ~ x[2:npred], ylim = c(0.0, 2.0), xlim = ylim, xlab = "dSST", ylab = "Loess Slope", pch = 16, cex = 0.8, type = "o", col = "red", sub = lo1_text)

# Model 7: Smoothing spline ===========================================================================================

spar = 0.95
ss1 <- smooth.spline(dMOT ~ dSST, spar = spar, keep.data = TRUE)
ss1
ss1_fit <- predict(ss1, dSST, deriv = 0)$y
ss1_residuals <- dMOT - ss1_fit
ss1_slope <- predict(ss1, dSST, deriv = 1)$y

# extract components of the AIC
df <- ss1$df 
sigma2 <- sum(ss1_residuals^2)/(n - 1.0)
rss <- sum(ss1_residuals^2 )
# http://users.stat.umn.edu/~helwig/notes/smooth-spline-notes.html#aic-and-bic
max_log_lik <- (-1.0 / (2.0 * sigma2)) * rss - 1.0 * (n / 2.0) * log(sigma2) - 1.0 * (n / 2.0) * log(2 * pi)
AIC_ss1 <- - 2.0 * max_log_lik + 2.0 * df
AIC_ss1 <- round(AIC_ss1, 2)
AIC_ss1

# R-sq values
rss <- sum(ss1_residuals^2 )
tss <- sum((dMOT - mean(dMOT))^2)
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

# get prediction intervals and confidence intervals
res <- (ss1$y - ss1$yin)/(1 - ss1$lev)
sig <- sqrt(var(res))
ci_upper <- ss1$y + 2.0 * sig * sqrt(ss1$lev)
ci_lower <- ss1$y - 2.0 * sig * sqrt(ss1$lev)

# smoothing spline scatter diagram
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, sub = ss1_text)
lines(ss1_fit ~ dSST, lwd = 2, col = "red")
lines(ci_upper ~ ss1$x, lty=2, col = "red")
lines(ci_lower ~ ss1$x, lty=2, col = "red")
lines(ss1_fit + 2.0 * sig ~ dSST, lty=2, col = "black")
lines(ss1_fit - 2.0 * sig ~ dSST, lty=2, col = "black")
dev.print(pdf, paste(plotpath, "smooth_spline.pdf", sep = ""))

# smoothin spline residual plot
oldpar <- par(pty="s")
plot(ss1_residuals ~ ss1_fit, ylim = ylim_resid, xlim = ylim, pch = 16, cex = 0.8, 
     ylab = "Residuals from Smoothing Spline", xlab =  "Smoothing Spline Fitted Values", sub = ss1_text)
lo_ss1_resid <- loess(ss1_residuals ~ ss1_fit, span=0.90, degree = 1, iterations = 3)
lines(lo_ss1_resid$fitted ~ ss1_fit, col = "red", lwd = 2)
dev.print(pdf, paste(plotpath, "smooth_spline_residual_plot.pdf", sep = ""))

# residual plot data
smooth_spline_residuals_out <- data.frame(cbind(ss1_fit, ss1_residuals, ss1_fit, lo_ss1_resid$fitted))
names(smooth_spline_residuals_out) <- c("smoothing_spline_fitted", "smoothing_spline_residuals", "smoothed_x", "smoothed_y")
out_csv <- "smooth_spline_residual_plot.csv"
write.csv(smooth_spline_residuals_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)

# slope plot from input data
oldpar <- par(pty="s")
plot(ss1_slope ~ dSST, ylim = c(0.0, 2.0), xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, type = "o", col = "red", sub = ss1_text)

# generate predictions over range of dSST
pred_range <- range(ocean_temps$dSST)
pred_range <- round(pred_range, 2)
x <- seq(pred_range[1], pred_range[2], by =  0.01) 

# to get predictions and C.I.'s over the range of the data, the ss1 leverages are used
oldpar <- par(pty="s")
plot(ss1$lev ~  ss1$x, ylim = c(0.0, 1.0), xlim = ylim, ylab = "Leverage", xlab = SST_label, pch = 16, cex = 0.8, sub = ss1_text)

# the leverage values have occasional outliers, so smooth them
# get smooth leverages (hat-matrix diagonal values) (using (what else?) a smoothing spline)
lev_smooth <- smooth.spline(ss1$lev ~ ss1$x, spar = 0.675, keep.data = TRUE)
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

# also get the local slope (first derivative)

# slope plot
ss1_pred_slope <- predict(ss1, x, deriv = 1)$y

oldpar <- par(pty="s")
plot(ss1_pred_slope ~ x, ylim = c(0.0, 2.0), type = "l", xlim = ylim, xlab = SST_label,
     ylab = "Smoothing Spline Slope", pch = 16, sub = ss1_text)
points(ss1_slope ~ dSST, pch = 16, cex = 0.6, col = "red")

ss_out <- data.frame(cbind(x, ss1_pred, pred_lwr, pred_upr, conf_lwr, conf_upr, ss1_pred_slope))
names(ss_out) <- c("dSST", "ss_fit", "ss_pred_lwr", "ss_pred_upr", "ss_conf_lwr", "ss_conf_upr", "ss_slope")
head(ss_out)
out_csv <- "smooth_spine_predictions.csv"
write.csv(ss_out, paste(csv_path, out_csv, sep = ""), row.names = FALSE)


# plot all four models
# linear_segmented_loess_smooth_spline
oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = ylim, xlim = ylim, ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8,
     sub = "Linear, segmented-regression, loess, and smoothing spline fits")
abline(a = 0.0, b = 1.0, col = "gray")
lines(lm1$fitted ~ dSST, lwd = 2, col = "red")
lines(lo1$fitted ~ dSST, lwd = 2, col = "blue")
lines(lm1_seg$fitted.values ~ dSST, col="purple", lwd=2)
lines(ss1_fit ~ dSST, lwd = 2, col = "magenta")
legend("bottomright", c(lm1_text, lm1_seg_text, lo1_text, ss1_text), 
       col = c("red", "purple", "blue", "magenta"), lwd = 2, cex = 0.85)
dev.print(pdf, paste(plotpath, "linear_segmented_loess_smooth_spline.pdf", sep = ""))

oldpar <- par(pty="s")
plot(dMOT ~ dSST, ylim = c(-4, 6), xlim = c(-4, 6), ylab = MOT_label, xlab = SST_label, pch = 16, cex = 0.8, 
     sub = "Linear, segmented-regression, loess, and smoothing spline fits")
abline(a = 0.0, b = 1.0, col = "gray")
lines(lm1$fitted ~ dSST, lwd = 2, col = "red")
lines(lo1$fitted ~ dSST, lwd = 2, col = "blue")
lines(lm1_seg$fitted.values ~ dSST, col="purple", lwd=2)
lines(ss1_fit ~ dSST, lwd = 2, col = "magenta")
legend("bottomright", c(lm1_text, lm1_seg_text, lo1_text, ss1_text), 
       col = c("red", "purple", "blue", "magenta"), lwd = 2, cex = 0.85)
dev.print(pdf, paste(plotpath, "linear_segmented_loess_smooth_spline_zoom.pdf", sep = ""))

