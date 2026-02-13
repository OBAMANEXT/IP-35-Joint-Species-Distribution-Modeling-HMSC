#lisab graafikutele toorandmete keskmised (20 võrdtõenäoses vahemikus)
#teisendab logaritmitud tunnused tagasi algskaalale

library(abind)
plotGradientlog=function (hM, Gradient, predY, measure, xlabel = NULL, ylabel = NULL, 
    index = 1, q = c(0.025, 0.5, 0.975), cicol = rgb(0, 0, 1, 
        alpha = 0.5), pointcol = "lightgrey", pointsize = 1, 
    showData = FALSE, jigger = 0, yshow = NA, showPosteriorSupport = TRUE, 
    main, ...) 
{
    Pr = NA
    if (is.null(xlabel)) {
        switch(class(hM$X)[1L], matrix = {
            xlabel = colnames(Gradient$XDataNew)[[1]]
        }, list = {
            xlabel = colnames(Gradient$XDataNew[[1]])[[1]]
        })
    }
    switch(class(hM$X)[1L], matrix = {
        xx = Gradient$XDataNew[, 1]
    }, list = {
        if (measure == "Y") {
            xx = Gradient$XDataNew[[index]][, 1]
        } else {
            xx = Gradient$XDataNew[[1]][, 1]
        }
    })
    ngrid = length(xx)
    if (measure == "S") {
        predS = abind(lapply(predY, rowSums), along = 2)
        Pr = mean(predS[ngrid, ] > predS[1, ])
        qpred = apply(predS, c(1), quantile, probs = q, na.rm = TRUE)
        if (is.null(ylabel)) {
            ylabel = "Summed response"
            if (all(hM$distr[, 1] == 2)) {
                ylabel = "Species richness"
            }
            if (all(hM$distr[, 1] == 3)) {
                ylabel = "Total count"
            }
        }
    }
    if (measure == "Y") {
        tmp = abind(predY, along = 3)
        Pr = mean(tmp[ngrid, index, ] > tmp[1, index, ])
        qpred = apply(tmp, c(1, 2), quantile, probs = q, na.rm = TRUE)
        qpred = qpred[, , index]
        if (is.null(ylabel)) {
            ylabel = hM$spNames[[index]]
        }
    }
    if (measure == "T") {
        if (all(hM$distr[, 1] == 1)) {
            predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)), 
                hM$nt), ncol = hM$nt))
        }
        else {
            predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a), 
                hM$nt), ncol = hM$nt))
        }
        predT = abind(predT, along = 3)
        Pr = mean(predT[ngrid, index, ] > predT[1, index, ])
        qpred = apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)
        qpred = qpred[, , index]
        if (is.null(ylabel)) {
            ylabel = hM$trNames[[index]]
        }
    }
    lo = qpred[1, ]
    hi = qpred[3, ]
    me = qpred[2, ]
    lo1 = min(lo, yshow, na.rm = TRUE)
    hi1 = max(hi, yshow, na.rm = TRUE)
    if (showData) {
        switch(class(hM$X)[1L], matrix = {
            XDatacol = which(colnames(Gradient$XDataNew)[[1]] == 
                colnames(hM$XData))
        }, list = {
            XDatacol = which(colnames(Gradient$XDataNew[[1]])[[1]] == 
                colnames(hM$XData[[1]]))
        })
        if (measure == "S") {
            pY = rowSums(hM$Y, na.rm = TRUE)
        }
        if (measure == "Y") {
            pY = hM$Y[, index]
        }
        if (measure == "T") {
            if (all(hM$distr[, 1] == 1)) {
                tmp = (exp(hM$Y) %*% hM$Tr)/matrix(rep(rowSums(exp(hM$Y)), 
                  hM$nt), ncol = hM$nt)
            }
            else {
                tmp = (hM$Y %*% hM$Tr)/matrix(rep(rowSums(hM$Y), 
                  hM$nt), ncol = hM$nt)
            }
            pY = tmp[, index]
        }
        if (!is.numeric(pY)) 
            pY <- as.numeric(pY)
        switch(class(hM$X)[1L], matrix = {
            pX = hM$XData[, XDatacol]
        }, list = {
            pX = hM$XData[[1]][, XDatacol]
        })
        hi1 <- max(hi1, max(pY, na.rm = TRUE))
        lo1 <- min(lo1, min(pY, na.rm = TRUE))
    }
    if (is.factor(xx)) {
        toPlot = data.frame(xx, me, lo, hi, stringsAsFactors = TRUE)
        pl = ggplot(toPlot, aes_string(x = xx, y = me)) + geom_bar(position = position_dodge(), 
            stat = "identity") + xlab(xlabel) + ylab(ylabel) + 
            geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, 
                position = position_dodge(0.9))
        if (showData) {
            if (jigger > 0) {
                pX = as.numeric(pX)
                pX = pX + runif(n = length(pY), min = -jigger, 
                  max = jigger)
            }
            dataToPlot = data.frame(pX = pX, pY = pY, stringsAsFactors = TRUE)
            pl = pl + geom_point(data = dataToPlot, aes_string(x = pX, 
                y = pY), size = pointsize)
        }
        if (!missing(main)) 
            pl = pl + labs(title = main)
    }
    else {
        if (inherits(hM$X, "list") && !measure == "Y") {
            plot(exp(xx), qpred[2, ], ylim = c(lo1, hi1), type = "l", 
                xaxt = "n", xlab = xlabel, ylab = ylabel, ...)
            axis(1, c(min(xx), (min(xx) + max(xx))/2, max(xx)), 
                c("min", "mean", "max"))
        }
        else {
            plot(exp(xx), qpred[2, ], ylim = c(lo1, hi1), type = "l", 
                xlab = xlabel, ylab = ylabel, ...)
        }
        if (showData) {
            if (jigger > 0) {
                de = (hi1 - lo1) * jigger
                pY = lo1 + de + (hi1 - lo1 - 2 * de) * (pY - 
                  lo1)/(hi1 - lo1) + runif(n = length(pY), min = -jigger, 
                  max = jigger)
            }
            if(esinemine_v_biomass!="logbiomass"){dataToPlot = cbind(exp(pX), pY)}else{dataToPlot = cbind(exp(pX), exp(pY)-1)}
            points(dataToPlot, pch = 16, col = pointcol)
        }
        polygon(c(exp(xx), rev(exp(xx))), c(qpred[1, ], rev(qpred[3, ])), 
            col = cicol, border = FALSE)
        lines(exp(xx), qpred[2, ], lwd = 2)
        if (!missing(main)) 
            title(main = main)
    }
    if (showPosteriorSupport && !is.factor(xx)) {
        mtext(gettextf("Pr[pred(%s=min) %s pred(%s=max)] = %.2f", 
            xlabel, ifelse(Pr < 0.5, ">", "<"), xlabel, ifelse(Pr < 
                0.5, 1 - Pr, Pr)))
    }
klassid=as.numeric(.bincode(pX,quantile(pX,c(0,1:20/20)),include.lowest=T,right=T)) #tegelikud vaatlused tükkideks (20 võrdse jaotusega klassi)
table=data.table(klassid,pX,pY)
mediaanid=table[,median(pX),keyby=klassid]
if(esinemine_v_biomass!="logbiomass"){osakaalud=table[,mean(pY),keyby=klassid]}else{osakaalud=table[,mean(exp(pY)-1),keyby=klassid]}
points(exp(mediaanid$V1),osakaalud$V1,col="red",pch=".",cex=4)
    if (is.factor(xx)) 
        pl
    else Pr
}



plotGradientmod=function (hM, Gradient, predY, measure, xlabel = NULL, ylabel = NULL, 
    index = 1, q = c(0.025, 0.5, 0.975), cicol = rgb(0, 0, 1, 
        alpha = 0.5), pointcol = "lightgrey", pointsize = 1, 
    showData = FALSE, jigger = 0, yshow = NA, showPosteriorSupport = TRUE, 
    main, ...) 
{
    Pr = NA
    if (is.null(xlabel)) {
        switch(class(hM$X)[1L], matrix = {
            xlabel = colnames(Gradient$XDataNew)[[1]]
        }, list = {
            xlabel = colnames(Gradient$XDataNew[[1]])[[1]]
        })
    }
    switch(class(hM$X)[1L], matrix = {
        xx = Gradient$XDataNew[, 1]
    }, list = {
        if (measure == "Y") {
            xx = Gradient$XDataNew[[index]][, 1]
        } else {
            xx = Gradient$XDataNew[[1]][, 1]
        }
    })
    ngrid = length(xx)
    if (measure == "S") {
        predS = abind(lapply(predY, rowSums), along = 2)
        Pr = mean(predS[ngrid, ] > predS[1, ])
        qpred = apply(predS, c(1), quantile, probs = q, na.rm = TRUE)
        if (is.null(ylabel)) {
            ylabel = "Summed response"
            if (all(hM$distr[, 1] == 2)) {
                ylabel = "Species richness"
            }
            if (all(hM$distr[, 1] == 3)) {
                ylabel = "Total count"
            }
        }
    }
    if (measure == "Y") {
        tmp = abind(predY, along = 3)
        Pr = mean(tmp[ngrid, index, ] > tmp[1, index, ])
        qpred = apply(tmp, c(1, 2), quantile, probs = q, na.rm = TRUE)
        qpred = qpred[, , index]
        if (is.null(ylabel)) {
            ylabel = hM$spNames[[index]]
        }
    }
    if (measure == "T") {
        if (all(hM$distr[, 1] == 1)) {
            predT = lapply(predY, function(a) (exp(a) %*% hM$Tr)/matrix(rep(rowSums(exp(a)), 
                hM$nt), ncol = hM$nt))
        }
        else {
            predT = lapply(predY, function(a) (a %*% hM$Tr)/matrix(rep(rowSums(a), 
                hM$nt), ncol = hM$nt))
        }
        predT = abind(predT, along = 3)
        Pr = mean(predT[ngrid, index, ] > predT[1, index, ])
        qpred = apply(predT, c(1, 2), quantile, probs = q, na.rm = TRUE)
        qpred = qpred[, , index]
        if (is.null(ylabel)) {
            ylabel = hM$trNames[[index]]
        }
    }
    lo = qpred[1, ]
    hi = qpred[3, ]
    me = qpred[2, ]
    lo1 = min(lo, yshow, na.rm = TRUE)
    hi1 = max(hi, yshow, na.rm = TRUE)
    if (showData) {
        switch(class(hM$X)[1L], matrix = {
            XDatacol = which(colnames(Gradient$XDataNew)[[1]] == 
                colnames(hM$XData))
        }, list = {
            XDatacol = which(colnames(Gradient$XDataNew[[1]])[[1]] == 
                colnames(hM$XData[[1]]))
        })
        if (measure == "S") {
            pY = rowSums(hM$Y, na.rm = TRUE)
        }
        if (measure == "Y") {
            pY = hM$Y[, index]
        }
        if (measure == "T") {
            if (all(hM$distr[, 1] == 1)) {
                tmp = (exp(hM$Y) %*% hM$Tr)/matrix(rep(rowSums(exp(hM$Y)), 
                  hM$nt), ncol = hM$nt)
            }
            else {
                tmp = (hM$Y %*% hM$Tr)/matrix(rep(rowSums(hM$Y), 
                  hM$nt), ncol = hM$nt)
            }
            pY = tmp[, index]
        }
        if (!is.numeric(pY)) 
            pY <- as.numeric(pY)
        switch(class(hM$X)[1L], matrix = {
            pX = hM$XData[, XDatacol]
        }, list = {
            pX = hM$XData[[1]][, XDatacol]
        })
        hi1 <- max(hi1, max(pY, na.rm = TRUE))
        lo1 <- min(lo1, min(pY, na.rm = TRUE))
    }
    if (is.factor(xx)) {
        toPlot = data.frame(xx, me, lo, hi, stringsAsFactors = TRUE)
        pl = ggplot(toPlot, aes_string(x = xx, y = me)) + geom_bar(position = position_dodge(), 
            stat = "identity") + xlab(xlabel) + ylab(ylabel) + 
            geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2, 
                position = position_dodge(0.9))
        if (showData) {
            if (jigger > 0) {
                pX = as.numeric(pX)
                pX = pX + runif(n = length(pY), min = -jigger, 
                  max = jigger)
            }
            dataToPlot = data.frame(pX = pX, pY = pY, stringsAsFactors = TRUE)
            pl = pl + geom_point(data = dataToPlot, aes_string(x = pX, 
                y = pY), size = pointsize)
        }
        if (!missing(main)) 
            pl = pl + labs(title = main)
    }
    else {
        if (inherits(hM$X, "list") && !measure == "Y") {
            plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", 
                xaxt = "n", xlab = xlabel, ylab = ylabel, ...)
            axis(1, c(min(xx), (min(xx) + max(xx))/2, max(xx)), 
                c("min", "mean", "max"))
        }
        else {
            plot(xx, qpred[2, ], ylim = c(lo1, hi1), type = "l", 
                xlab = xlabel, ylab = ylabel, ...)
        }
        if (showData) {
            if (jigger > 0) {
                de = (hi1 - lo1) * jigger
                pY = lo1 + de + (hi1 - lo1 - 2 * de) * (pY - 
                  lo1)/(hi1 - lo1) + runif(n = length(pY), min = -jigger, 
                  max = jigger)
            }
            if(esinemine_v_biomass!="logbiomass"){dataToPlot = cbind(pX, pY)}else{dataToPlot = cbind(pX, exp(pY)-1)}
            points(dataToPlot, pch = 16, col = pointcol)
        }
        polygon(c(xx, rev(xx)), c(qpred[1, ], rev(qpred[3, ])), 
            col = cicol, border = FALSE)
        lines(xx, qpred[2, ], lwd = 2)
        if (!missing(main)) 
            title(main = main)
    }
    if (showPosteriorSupport && !is.factor(xx)) {
        mtext(gettextf("Pr[pred(%s=min) %s pred(%s=max)] = %.2f", 
            xlabel, ifelse(Pr < 0.5, ">", "<"), xlabel, ifelse(Pr < 
                0.5, 1 - Pr, Pr)))
    }
klassid=as.numeric(.bincode(pX,quantile(pX,c(0,1:20/20)),include.lowest=T,right=T)) #tegelikud vaatlused tükkideks (20 võrdse jaotusega klassi)
table=data.table(klassid,pX,pY)
mediaanid=table[,median(pX),keyby=klassid]
if(esinemine_v_biomass!="logbiomass"){osakaalud=table[,mean(pY),keyby=klassid]}else{osakaalud=table[,mean(exp(pY)-1),keyby=klassid]}
points(mediaanid$V1,osakaalud$V1,col="red",pch=".",cex=4)
    if (is.factor(xx)) 
        pl
    else Pr
}



