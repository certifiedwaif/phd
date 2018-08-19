#################################################################################################### 

set.seed(1)

source("functions.R")

library(latex2exp)

n = 120

ln = vector("list", 6)
ln[[1]] = n
ln[[2]] = n
ln[[3]] = n
ln[[4]] = n
ln[[5]] = n
ln[[6]] = n

lp = vector("list", 6)
lp[[1]] = 3
lp[[2]] = 20
lp[[3]] = 100
lp[[4]] = 3
lp[[5]] = 20
lp[[6]] = 100

lR2 = vector("list", 6)
lR2[[1]] = 0.05
lR2[[2]] = 0.05
lR2[[3]] = 0.05
lR2[[4]] = 0.95
lR2[[5]] = 0.95
lR2[[6]] = 0.95

PLOT.G = TRUE
PLOT.U = TRUE
PLOT.SIGMA = TRUE
PLOT.ALPHA = TRUE

PDF = TRUE
CHECK = FALSE

#################################################################################################### 

if (PLOT.G) {
    if (PDF) {
        pdf("gGivenY.pdf", width = 12, height = 8)
    }
    
    par(mfrow = c(2, 3))
    
    for (i in 1:6) {
        
        n = ln[[i]]
        p = lp[[i]]
        R2 = lR2[[i]]
        
        den1 = posterior.gGivenY.MonteCarlo(n, p, R2, N = 1000)
        
        g = den1$x  # den1$samples
        g = g[g > 0]
        
        # eps = 0.001 r = quantile( g, 1-epsR )
        g = seq(min(c(1, g)), max(g), , 501)
        
        den2 = posterior.gGivenY.exact(g, n, p, R2)
        
        xlim = range(g)
        ylim = range(c(den1$y, den2$y))
        
        plot(NA, type = "n", xlim = xlim, ylim = ylim, ylab = "density", xlab = "g", 
            cex.lab = 1.5, cex.main = 1.5, main = paste("g|y   (n=", n, ", p=", p, 
                ", R2=", R2, ")", sep = ""))
        lines(den1, col = "red", lwd = 2, lty = 2)
        lines(den2, col = "black", lwd = 2, lty = 1)
        
        
        if (i == 1) {
            legend("topright", legend = c("Exact", "Samples"), col = c("black", "red"), 
                lwd = c(2, 2), lty = c(1, 2), cex = 1.5, bty = "n")
        }
        
        if (CHECK) {
            print(trapint(den1$x, den1$y))
            print(trapint(den2$x, den2$y))
        }
    }
    
    if (PDF) {
        dev.off()
    }
}

#################################################################################################### 

if (PLOT.U) {
    if (PDF) {
        pdf("uGivenY.pdf", width = 12, height = 8)
    }
    
    par(mfrow = c(2, 3))
    
    for (i in 1:6) {
        
        n = ln[[i]]
        p = lp[[i]]
        R2 = lR2[[i]]
        
        den1 = posterior.uGivenY.MonteCarlo(n, p, R2, N = 1000)
        
        u = den1$x  # den1$samples
        u = u[(u > 0) & (u < 1)]
        
        # eps = 0.001 r = quantile( u, 1-epsR )
        g = seq(min(u), max(u), , 501)
        
        den2 = posterior.uGivenY.exact(u, n, p, R2)
        
        xlim = range(u)
        ylim = range(c(den1$y, den2$y))
        
        plot(NA, type = "n", xlim = xlim, ylim = ylim, ylab = "density", xlab = "u", 
            cex.lab = 1.5, cex.main = 1.5, main = paste("u|y   (n=", n, ", p=", p, 
                ", R2=", R2, ")", sep = ""))
        lines(den1, col = "red", lwd = 2, lty = 2)
        lines(den2, col = "black", lwd = 2, lty = 1)
        
        if (i == 1) {
            legend("topleft", legend = c("Exact", "Samples"), col = c("black", "red"), 
                lwd = c(2, 2), lty = c(1, 2), cex = 1.5, bty = "n")
        }
        
        if (CHECK) {
            print(trapint(den1$x, den1$y))
            print(trapint(den2$x, den2$y))
        }
    }
    
    if (PDF) {
        dev.off()
    }
}

#################################################################################################### 

if (PLOT.SIGMA) {
    if (PDF) {
        pdf("sigma2GivenY.pdf", width = 12, height = 8)
    }
    
    par(mfrow = c(2, 3))
    
    for (i in 1:6) {
        
        n = ln[[i]]
        p = lp[[i]]
        R2 = lR2[[i]]
        
        den1 = posterior.sigma2GivenY.MonteCarlo(n, p, R2, N = 1000)
        
        sigma2 = den1$x  # den1$samples
        sigma2 = sigma2[sigma2 > 0]
        
        # eps = 0.001 r = quantile( sigma2, 1-epsR )
        sigma2 = seq(min(sigma2), max(sigma2), , 501)
        
        den2 = posterior.sigma2GivenY.exact(sigma2, n, p, R2)
        den3 = posterior.sigma2GivenY.RaoBlackwell(sigma2, n, p, R2, N = 100)
        den4 = posterior.sigma2GivenY.approx(sigma2, n, p, R2)
        
        xlim = range(sigma2)
        ylim = range(c(den1$y, den3$y, den4$y))
        
        plot(NA, type = "n", xlim = xlim, ylim = ylim, ylab = "density", xlab = expression(sigma^2), 
            cex.lab = 1.5, cex.main = 1.5, main = TeX(paste("$\\sigma^2|y$   (n=", 
                n, ", p=", p, ", R2=", R2, ")", sep = "")))
        lines(den1, col = "red", lwd = 1, lty = 1)
        lines(den2, col = "black", lwd = 3, lty = 3)
        lines(den3, col = "blue", lwd = 3, lty = 3)
        lines(den4, col = "green", lwd = 3, lty = 3)
        
        if (i == 1) {
            legend("topright", legend = c("Exact", "Samples", "RB", "Delta"), col = c("black", 
                "red", "blue", "green"), lwd = c(3, 1, 3, 3), lty = c(3, 1, 3, 3), 
                cex = 1.5, bty = "n")
        }
        
        if (CHECK) {
            print(trapint(den1$x, den1$y))
            print(trapint(den2$x, den2$y))
            print(trapint(den3$x, den3$y))
            print(trapint(den4$x, den4$y))
        }
    }
    
    if (PDF) {
        dev.off()
    }
}

#################################################################################################### 

if (PLOT.ALPHA) {
    
    if (PDF) {
        pdf("alphaGivenY.pdf", width = 12, height = 8)
    }
    
    par(mfrow = c(2, 3))
    
    for (i in 1:6) {
        
        n = ln[[i]]
        p = lp[[i]]
        R2 = lR2[[i]]
        
        den1 = posterior.alphaGivenY.MonteCarlo(n, p, R2, N = 1000)
        
        alpha = den1$x  # den1$samples
        
        # eps = 0.001 r = quantile( sigma2, 1-epsR )
        alpha = seq(min(alpha), max(alpha), , 501)
        
        den2 = posterior.alphaGivenY.exact(alpha, n, p, R2)
        den3 = posterior.alphaGivenY.RaoBlackwell(alpha, n, p, R2, N = 100)
        den4 = posterior.alphaGivenY.approx(alpha, n, p, R2)
        
        xlim = range(alpha)
        ylim = range(c(den1$y, den3$y, den4$y))
        
        plot(NA, type = "n", xlim = xlim, ylim = ylim, ylab = "density", xlab = expression(alpha), 
            cex.lab = 1.5, cex.main = 1.5, main = TeX(paste("$\\alpha |y$   (n=", 
                n, ", p=", p, ", R2=", R2, ")", sep = "")))
        lines(den1, col = "red", lwd = 1, lty = 1)
        lines(den2, col = "black", lwd = 3, lty = 3)
        lines(den3, col = "blue", lwd = 3, lty = 3)
        lines(den4, col = "green", lwd = 3, lty = 3)
        
        if (i == 1) {
            legend("topright", legend = c("Exact", "Samples", "RB", "Delta"), col = c("black", 
                "red", "blue", "green"), lwd = c(3, 1, 3, 3), lty = c(3, 1, 3, 3), 
                cex = 1.5, bty = "n")
        }
        
        if (CHECK) {
            print(trapint(den1$x, den1$y))
            print(trapint(den2$x, den2$y))
            print(trapint(den3$x, den3$y))
            print(trapint(den4$x, den4$y))
        }
    }
    
    if (PDF) {
        dev.off()
    }
}

