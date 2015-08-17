comp_kids <- function(simk, imputek){
    #
    out <- data.frame()
    for(i in 1:length(simk)){
        
        imputeki <- imputek[[i]][[3]]
        simki <- simk[[i]][[2]][imputeki$idx, ]
        
        names(simki) <- c("shap1", "shap2", "obs")
        comb <- cbind(simki, imputeki)
        
        err = 0
        for(j in unique(comb$chunk)){
            sub <- subset(comb, chunk == j)
            idx1 <- which.max(c(cor(sub$k1, sub$shap1), cor(sub$k1, sub$shap2)) )
            err1 <- sum(sub$k1 != sub[, idx1])
            idx2 <- which.max(c(cor(sub$k2, sub$shap1), cor(sub$k2, sub$shap2)) )
            err2 <- sum(sub$k2 != sub[, idx2])
            err <- err + err1 + err2
        }
        tem <- data.frame(kid=i, err=err, tot=nrow(simki))
        out <- rbind(out, tem)
    }
    out$rate <- out$err/out$tot
    message(sprintf("###>>> Average error rate [ %s ]", mean(out$rate)))
    return(out)
}

