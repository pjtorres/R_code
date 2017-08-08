# Great way to combine the pairs(),cor(), and hist() functions into one great graph. THIS IS NOT MY CODE THIS IS TAKEN FROM ane handles man 
#AS SEEN IN THIS BLOG https://www.blogger.com/profile/07863739093681446569 and http://handlesman.blogspot.com/2011/03/matrix-plot-with-confidence-intervals.html

## put histograms on the diagonal
panel.hist <- function(x, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
rect(breaks[-nB], 0, breaks[-1], y, col="lavender", ...)
}
## put correlations & 95% CIs on the upper panels,
panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...)
{
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- cor(x, y,use="complete.obs")
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    prefix <- "r = "
    rc <- cor.test(x,y)
    rci <- rc$conf.int
    txt2 <- format(c(rci, 0.123456789), digits=digits)[1]
    txt3 <- format(c(rci, 0.123456789), digits=digits)[2]
    prefix2 <- "\nCI = "
    txt <- paste(prefix, txt, prefix2, txt2, ", ", txt3, sep="")
    if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
    text(0.5, 0.5, txt, cex = 1)
}
pairs(iris[1:4], lower.panel=panel.smooth, cex = .8, pch = 21, bg="steelblue",
       diag.panel=panel.hist, cex.labels = 1.2, font.labels=2, upper.panel=panel.cor)
