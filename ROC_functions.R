#####################################
###ROC curves####
#####################################


#code of Lennart Raman

exec.ROC <- function(truth, vals1, t1){
  o <- order(vals1)
  truth <- truth[o]
  vals1 <- vals1[o]
  
  c1 <- which(truth == t1)
  c2 <- which(truth != t1)
  
  TPs <- c() ; FPs <- c() ; TNs <- c() ; FNs <- c()
  for (s in vals1){
    these.c1 <- which(vals1 > s)
    these.c2 <- which(vals1 <= s)
    TPs <- c(TPs, length(which(these.c1 %in% c1)))
    FPs <- c(FPs, length(which(these.c1 %in% c2)))
    TNs <- c(TNs, length(which(these.c2 %in% c2)))
    FNs <- c(FNs, length(which(these.c2 %in% c1)))
  }
  return(list(FPRs = c(1, FPs / (FPs + TNs), 0), TPRs = c(1, TPs / (TPs + FNs), 0), o, truth, vals1))
  
}
exec.combineROC <- function(r1, r2, r3){
  FPRss = c()
  for (i in 1:46){
    avg= signif((r1$FPRs[i]+r2$FPRs[i]+r3$FPRs[i])/3, digits=7)
    FPRss = c(FPRss, avg)
  }
  TPRss = c()
  for (i in 1:46){
    avg= (r1$TPRs[i]+r2$TPRs[i]+r3$TPRs[i])/3
    TPRss = c(TPRss, avg)
  }
  
  avgresult = list(FPRs = FPRss, TPRs = TPRss)
  return(avgresult)
}
simple.auc <- function(roc){ # https://www.r-bloggers.com/calculating-auc-the-area-under-a-roc-curve/
  dFPR <- c(diff(roc$FPRs), 0)
  dTPR <- c(diff(roc$TPRs), 0)
  auc = -(sum(roc$TPRs * dFPR) + sum(dTPR * dFPR)/2)
  return(round(auc, 3))
}

plot.ROC <- function(roc, col='black', lty=1, add=F, cex.lab = 1,
                     line = 2, lwd = 1, plt.mtext1 = T, plt.mtext2 = T){
  if (add) par(new=T)
  plot(roc$FPRs, roc$TPRs, axes = F, ylab='', xlab='', type='l', col=col, lwd=lwd, lty=lty, xlim=c(0,1), ylim=c(0,1))
  polygon(c(rev(roc$FPRs), 0), c(rev(roc$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  if (!add){
    if (plt.mtext1){
      axis(1, las=1, tcl=0.5)
      mtext('false positive rate', 1, line, cex = cex.lab)
    }
    if (plt.mtext2){
      axis(2, las=1, tcl=0.5)
      mtext('true positive rate', 2, line + .3, cex = cex.lab)
    }
    segments(0,0,1,1, lwd=.05)
  }
}


#own code

plot.ROC.combined <- function(roc1, roc2, roc3, roc4="x", roc5=F, title="", col='black', add=F, cex.lab = 1.5,
                              line = 2, lwd = 1.5, plt.mtext1 = T, plt.mtext2 = T, plt.mtext3 = T){
  if (add) par(new=T)
  plot(roc1$FPRs, roc1$TPRs, axes = F, ylab='', xlab='', type='l', col="#66c2a5", lwd=lwd, lty=1, xlim=c(0,1), ylim=c(0,1))
  polygon(c(rev(roc1$FPRs), 0), c(rev(roc1$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  par(new=T)
  plot(roc2$FPRs, roc2$TPRs, axes = F, ylab='', xlab='', type='l', col="#fc8d62", lwd=lwd, lty=1, xlim=c(0,1), ylim=c(0,1))
  polygon(c(rev(roc2$FPRs), 0), c(rev(roc2$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  par(new=T)
  plot(roc3$FPRs, roc3$TPRs, axes = F, ylab='', xlab='', type='l', col="#8da0cb", lwd=lwd, lty=4, xlim=c(0,1), ylim=c(0,1))
  polygon(c(rev(roc3$FPRs), 0), c(rev(roc3$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  par(new=T)
  if (roc4!="x"){
    plot(roc4$FPRs, roc4$TPRs, axes = F, ylab='', xlab='', type='l', col="#e78ac3", lwd=lwd, lty=3, xlim=c(0,1), ylim=c(0,1))
    polygon(c(rev(roc4$FPRs), 0), c(rev(roc4$TPRs), 0), border = NA, col = adjustcolor(col, .04))
    par(new=T)
    plot(roc5$FPRs, roc5$TPRs, axes = F, ylab='', xlab='', type='l', col="#a6d854", lwd=lwd, lty=3, xlim=c(0,1), ylim=c(0,1))
    polygon(c(rev(roc5$FPRs), 0), c(rev(roc5$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  }
  if (!add){
    if (plt.mtext1){
      axis(1, las=1, tcl=0.5)
      mtext('false positive rate', 1, line, cex = cex.lab)
    }
    if (plt.mtext2){
      axis(2, las=1, tcl=0.5)
      mtext('true positive rate', 2, line + .3, cex = cex.lab)
    }
    if (plt.mtext3){
      mtext(title, outer = F, cex = cex.lab+0.5)
    }
    segments(0,0,1,1, lwd=.05)
  }
}

plot.ROC.final <- function(roc, col='black', lty=1, add=F, cex.lab = 1,
                     line = 2, lwd = 1, plt.mtext1 = T, plt.mtext2 = T){
  if (add) par(new=T)
  plot(roc$FPRs, roc$TPRs, axes = F, ylab='', xlab='', type='l', col="#a6d854", lwd=lwd, lty=lty, xlim=c(0,1), ylim=c(0,1))
  polygon(c(rev(roc$FPRs), 0), c(rev(roc$TPRs), 0), border = NA, col = adjustcolor(col, .04))
  if (!add){
    if (plt.mtext1){
      axis(1, las=1, tcl=0.5)
      mtext('false positive rate', 1, line, cex = cex.lab)
    }
    if (plt.mtext2){
      axis(2, las=1, tcl=0.5)
      mtext('true positive rate', 2, line + .3, cex = cex.lab)
    }
    segments(0,0,1,1, lwd=.05)
  }
}


