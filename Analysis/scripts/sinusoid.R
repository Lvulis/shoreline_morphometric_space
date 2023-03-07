## Generate synthetic sinusoids for curvature example in supplementary material
## Author: lvulis
## Date: 2022 - 06 - 22

# Change to appropriate wording directory
# dirname <- "~/deltas/shorelines/"
dirname <- "G:/My Drive/NSFDELTAS/Analysis/"
setwd(dirname)

source("./scripts/shoreline_paper/utils.R")
## Construct timeseries

t <- seq(0, 100*pi, by = 1)
A1 = 5
A2 = 20
x1 <- A1*sin(t) + rnorm(length(t), 0, A1/10)
x2 <- A2*sin(t/2/pi) + rnorm(length(t), 0, A2/10)

# Curvature mapping
ksine <- compute_curvature(cbind(t, x1+x2))
ksine_fast <- compute_curvature(cbind(t, x1))
ksine_slow <- compute_curvature(cbind(t, x2))

# plot(t, x1)
# lines(t, x2)
# filename = paste0(dirname, "/plots/sine_ex/sinusoid.png")
# png(file = filename, width = 1400, height = 700, pointsize = 14)

filename = paste0(dirname, "/plots/sine_ex/sinusoid.pdf")
pdf(file = filename, width = 10, height = 7.5, pointsize = 14)
par(mar = c(5.1, 6.1, 4.1, 2.1), mfrow = c(3, 2))

plot(t, x1, type = 'l', main = bquote(z[1]),
     ylab = bquote(z[1]),
     xlab = bquote(italic(s)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2)

psd1= spectrum(x1, plot = F)
plot(psd1$freq, psd1$spec/1e3, type = 'l', log = '', 
     main = bquote(widehat(z)[1]),
     ylab = bquote("PSD("*italic(k)*")"~10^3),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2 )

plot(t, x2, type = 'l', main = bquote(z[2]),
     ylab = bquote(z[2]),
     xlab = bquote(italic(s)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2)

psd2 = spectrum(x2, plot = F)
plot(psd2$freq, psd2$spec/1e3, type = 'l', log = '', 
     main = bquote(widehat(z)[2]),
     ylab = bquote("PSD("*italic(k)*")"~10^3),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2 )


plot(t, x1+x2, type = 'l', main = bquote(z[3] == x[1] + x[2]),
     ylab = bquote(z[3]),
     xlab = bquote(italic(s)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2)

psd3 = spectrum(x1 + x2, plot = F)
plot(psd3$freq, psd3$spec/1e3, type = 'l', log = '', 
     main = bquote(widehat(z)[3]),
     ylab = bquote("PSD("*italic(k)*")"~10^3),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3, las = 1,
     lwd = 2 )

dev.off()


filename = paste0(dirname, "/plots/sine_ex/sinusoid_curvature.pdf")
# png(file = filename, width = 1400, height = 700, pointsize = 14)
# par(mar = c(5.1, 5.1, 4.1, 2.1), mfrow = c(3, 2))
pdf(file = filename, width = 10, height = 7.5, pointsize = 14)
par(mar = c(5.1, 6.1, 4.1, 2.1), mfrow = c(3, 2))

plot(t[-c(1, length(t))], ksine_fast, type = 'l',
     xlab = bquote(italic(l)),  las = 1,
     ylab = "",
     main = bquote(kappa[1]*"("*italic(l)*")"),
     cex.axis = 2, cex.lab = 2, cex.main = 3)
title(ylab = bquote(kappa), line = 4, cex.lab = 2)

psd1= spectrum(ksine_fast, plot = F)
plot(psd1$freq, psd1$spec, type = 'l', log = '',  las = 1,
     main = bquote(widehat(kappa)[1]),
     ylab = bquote("PSD("*italic(k)*")"),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3,
     lwd = 2 )


plot(t[-c(1, length(t))], ksine_slow, type = 'l', xlab = bquote(italic(l)), 
     ylab = "",
     las = 1, main = bquote(kappa[2]*"("*italic(l)*")"),
     cex.axis = 2, cex.lab = 2, cex.main = 3)
title(ylab = bquote(kappa), line = 4, cex.lab = 2)

psd2= spectrum(ksine_slow, plot = F)
plot(psd2$freq, psd2$spec, type = 'l', log = '',  las = 1,
     main = bquote(widehat(kappa)[2]),
     ylab = bquote("PSD("*italic(k)*")"),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3,
     lwd = 2 )

plot(t[-c(1, length(t))], ksine, type = 'l', 
     xlab = bquote(italic(l)), 
     ylab = "", 
     main = bquote(kappa[3]*"("*italic(l)*")"),
     cex.axis = 2, cex.lab = 2,  las = 1, cex.main = 3)
title(ylab = bquote(kappa), line = 4, cex.lab = 2)

psd3= spectrum(ksine, plot = F)
plot(psd3$freq, psd3$spec, type = 'l', log = '', 
     main = bquote(widehat(kappa)[3]), las = 1,
     ylab = bquote("PSD("*italic(k)*")"),
     xlab = bquote(italic(k)),
     cex.axis = 2, cex.lab = 2, cex.main = 3,
     lwd = 2 )

dev.off()

