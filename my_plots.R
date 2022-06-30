adjacency.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
  JJ <- dim(mat)[1]
  colorscale <- c("white", rev(heat.colors(100)))
  if(missing(labs))     labs <- 1:JJ
  if(missing(col.axis)) col.axis <- rep("black", JJ)
  if(missing(cex.axis)) cex.axis <- 0.5
  if(missing(tick))     tick <- TRUE
  ## adjacency matrix
  image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
        col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
  for(j in 1:JJ){
    axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
    axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
  }
  box()
  if(show.grid) grid(nx=JJ, ny=JJ)
}
 
heat.plot0 <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis,...)
{ 
  JJ <- dim(mat)[1]
  colorscale <- c("white", rev(heat.colors(100)))
  if(missing(labs))     labs <- 1:JJ
  if(missing(col.axis)) col.axis <- rep("black", JJ)
  if(missing(cex.axis)) cex.axis <- 0.5
  if(missing(tick))     tick <- TRUE
  ## adjacency matrix
  image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
        col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
  for(j in 1:JJ){
    axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
    axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
  }
  box()
  if(show.grid) grid(nx = JJ, ny = JJ)
}
 
heat.plot <- function (mat, show.grid = FALSE, cex.axis, tick, labs, col.axis, ...)
{ 
  JJ <- dim(mat)[1]
  colorscale <- c("white", rev(heat.colors(100)))
  if(missing(labs))     labs <- 1:JJ
  if(missing(col.axis)) col.axis <- rep("black", JJ)
  if(missing(cex.axis)) cex.axis <- 0.5
  if(missing(tick))     tick <- TRUE
  ## adjacency matrix
  par(layout(mat = matrix(1:2, 1, 2), widths = c(9, 1)), mar = c(4, 2.5, 2, 1))
  image(seq(1, JJ), seq(1, JJ), mat, axes = FALSE, xlab = "", ylab = "",
        col = colorscale[seq(floor(100*min(mat)), floor(100*max(mat)))], ...)
  for(j in 1:JJ){
    axis(1, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
    axis(2, at = j, labels = labs[j], las = 2, cex.axis = cex.axis, tick, col.axis = col.axis[j])
  }
  box()
  if(show.grid) grid(nx = JJ, ny = JJ)
  ## scale
  par(mar = c(3.1, 0.1, 1.3, 0.1))
  plot(1:100, 1:100, xlim = c(0,2), ylim = c(0, 100), axes = FALSE, type = "n", 
       xlab = "", ylab = "")
  yposr <- 1:100
  rect(0, yposr-.5, 0.5, yposr+.5, col = rev(heat.colors(100)), border = FALSE)
  rect(0, 0.5, 0.5, 100.5, col = "transparent")
  text(0.42, c(yposr[1], yposr[25], yposr[50], yposr[75], yposr[100]), 
       c("0", "0.25", "0.5", "0.75", "1"), pos = 4, cex = 0.6)
}
 
chain.plot <- function (x, plot.mean = TRUE, plot.ci = TRUE, xtruth,
                        alpha = 0.05, ... )
{
  N <- length(x)
  plot(1:N, x, type = "l", xlab = "Iteration", col = "lightgray", ... )
  if (!missing(xtruth))
    abline(h = xtruth, col = "blue", lty = 1, lwd = 1)
  if (missing(xtruth)) {
    pmcol <- "blue"
    pmlwd <- 1
  } else {
    CI <- quantile(x, probs = c(alpha/2, 1 - alpha/2))
    if ((CI[1] < xtruth) & (xtruth < CI[2])) {
      pmcol <- "green3"
      pmlwd <- 2
    } else {
      pmcol <- "red3"
      pmlwd <- 1
    }
    
  }
  if (plot.mean) 
    abline(h = mean(x), lty = 2, col = pmcol, lwd = pmlwd)
  if (plot.ci)   
    abline(h = quantile(x, probs = c(alpha/2, 1 - alpha/2)), lty = 2, col = pmcol, lwd = pmlwd)
}
 
post.plot <- function (x, main, xtruth, ylim, alpha = 0.05, adj = 1.15, 
                       lwd.density = 1, tck = 0.015, col.border = "gray", 
                       col.shade = "gray95", info = FALSE, 
                       summaries = TRUE, plot.data = FALSE, ...)
{
  den  <- density(x, adj = adj, n = 10000) 
  xx   <- den$x
  yy   <- den$y
  if(missing(ylim))
    ylim <- 1.05*c(0, max(range(yy)[2], hist(x, plot = FALSE)$density)) 
  
  info <- ifelse(info, paste("n =",length(x),"   ","bandwith =", round(den$bw, 2)),"")
  if(missing(main)) main = ""
  
  hist(x, freq = FALSE, xaxt = "n", ann = FALSE, border = col.border, ylim = ylim, ... )
  polygon(c(xx,rev(xx)), c(rep(0,length(xx)), rev(yy)), lty = 1, col = col.shade, xaxt = "n", ann = FALSE)
  hist(x, freq = FALSE, xaxt = "n", ann = FALSE, border = col.border, add = TRUE)
  if (!missing(xtruth)) {
    abline(v = xtruth, col = "blue", lty = 1)
  }
  title(main)
  lines(den, lwd = lwd.density)
  
  if(plot.data) axis(side = 1, at = x, labels = FALSE, tck = tck)
  m <- abs(max(x))
  r <- 1
  while (m < 1) {
    m <- m*10^r
    r <- r + 1
  }
  xticks <- round(seq(min(x), max(x), len = 5), r)
  axis(side = 1, at = xticks, labels = xticks)
  
  if(summaries){
    closest <- function(x, x0) which(abs(x - x0) == min(abs(x - x0)))
    idx1    <- closest(xx, quantile(x, alpha/2))
    idx2    <- closest(xx, quantile(x, 1 - alpha/2))
    idx3    <- closest(xx, mean(x))
    idx     <- c(idx1, idx2, idx3)
    if(missing(xtruth)) {
      pmcol <- "blue"
      pmlwd <- 1
    } else {
      if ((xx[idx1] < xtruth) & (xtruth < xx[idx2])) {
        pmcol <- "green3"
        pmlwd <- 2
      } else {
        pmcol <- "red3"
        pmlwd <- 1
      }
    }
    segments(x0 = xx[idx], y0 = 0, x1 = xx[idx], y1 = yy[idx], col = pmcol, 
             lty = 2, lwd = pmlwd)
    axis(side = 1, at = xx[idx], col.ticks = pmcol, lwd.ticks = pmlwd, labels = FALSE)
  }
}
 
post.CIs.plot <- function(x, alpha = 0.05, xtruth, ylim, plot.xaxis = TRUE, 
                          cex.axis = 1, ntop = 25, ...)
{
  n  <- dim(x)[2]
  li <- apply(x, 2, quantile, probs = c(alpha/2))
  ul <- apply(x, 2, quantile, probs = c(1 - alpha/2))
  pm <- colMeans(x)
  
  if (!missing(xtruth)) {
    idx <- ((li <= xtruth) & (xtruth <= ul))
    ok  <- round(sum(idx)/n, 3) * 100
    notok <- round(100 - ok, 1)
    if(missing(ylim)) ylim <- 1.05 * range(li, ul, xtruth)
  } else {
    if(missing(ylim)) ylim <- 1.05 * range(li, ul)
  }
  
  if (n > ntop) {
    id.lab <- sort(sample(1:n, ntop))
    labs   <- NULL
    li <- li[id.lab]
    ul <- ul[id.lab]
    pm <- pm[id.lab]
    if (!missing(xtruth)) { 
      xtruth <- xtruth[id.lab] 
      idx <- idx[id.lab]
      if(missing(ylim)) ylim = 1.05 * range(li, ul, xtruth)
    } else {
      if(missing(ylim)) ylim = 1.05 * range(li, ul)
    }
    n <- ntop
  } else {
    labs <- 1:n
  }
  
  plot(1:n, ylim = ylim, type = "n", xaxt = "n", ...)
  abline(v = 1:n, lty = 3, col = "lightgray")
  if(plot.xaxis) axis(side = 1, at = 1:n, labels = labs, cex.axis = cex.axis)
  for (i in 1:n) {
    if(missing(xtruth)) {
      pmlwd <- 1
      pmcol <- "blue"
    } else {
      if ((li[i] < xtruth[i]) & (xtruth[i] < ul[i])) {
        pmlwd <- 2
        pmcol <- "green3"
      } else {
        pmlwd <- 1
        pmcol <- "red3"
      }
    }
    segments(i, li[i], i, ul[i], lwd = pmlwd, col = pmcol)
    points(i, pm[i], pch = 20, col = pmcol, cex = 1)
  }
  if (!missing(xtruth)) {
    lines(xtruth, type = "p", pch = 18, col = "blue")
    legend("top", legend = c(paste0("Right ", ok, "%"), paste0("Wrong ", notok, "%")), 
           text.col = c(3, 2), cex = 0.75, bty = "n", horiz = TRUE)
  }
}

IC.plot <- function (criteria, Kmin, Kmax, model, prefix, path.sams, path.puts) 
{
  IC.vals <- rep(NA, Kmax - Kmin + 1)
  for (K in Kmin:Kmax) IC.vals[K] <- as.numeric(read.table(paste0(path.sams, model, "_", prefix,  "/", criteria,"_K_", K, ".txt")))
  pdf(paste0(path.outs, model, "_", prefix, "_ic_", criteria, ".pdf"), height = 10, width = 10, pointsize = 12)
  par(mfrow = c(1, 1), mar = c(4, 4, 3, 2) - 0.0, mgp = c(2, 1, 0), oma = c(0, 0, 0, 0))
  plot(IC.vals, type = "b", lty = 1, lwd = 2, pch = 1, xaxt = "n", xlab = "Latent dimension, K", ylab = paste0(criteria, "(K)"), main = criteria)
  axis(side = 1, at = 1:(Kmax - Kmin + 1), labels = Kmin:Kmax, cex.axis = 0.8)
  grid()
  dev.off()
}