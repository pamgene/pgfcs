#' @import corpcor
#' @import combinat
#' @import dplyr
#' @import reshape2

DFT.PERMS = 500

#' @export
fcs = function(dataMatrix, classMatrix, phenoGrp, statFun = stat.snr, phenoPerms = TRUE, featurePerms = TRUE, nPerms = DFT.PERMS){
  upStats = setStats(dataMatrix, classMatrix, phenoGrp, statFun)
  nClass = apply(classMatrix > 0 ,2, sum)
  aFcsResult = data.frame(ClassName = colnames(classMatrix), nFeatures = nClass, SetStat = upStats, NormalizedSetStat = upStats/nClass)
  if(phenoPerms){
    aPhenoScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = fcsPhenoPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aPhenoScore = NaN
  }
  if(featurePerms){
    aFeatureScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = fcsSpotPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aFeatureScore = NaN
  }
  aFcsResult = data.frame(aFcsResult,
                          phenoScore = aPhenoScore,
                          featureScore = aFeatureScore,
                          pPhenoScore = -log10(aPhenoScore),
                          pFeatureScore = -log10(aFeatureScore) )
  combinedScore = apply(as.matrix(aFcsResult[, c("pPhenoScore", "pFeatureScore")]),1,sum, na.rm = TRUE)
  delta = setStats(dataMatrix, classMatrix, phenoGrp, statfun = stat.delta)
  aFcsResult = data.frame(aFcsResult, combinedScore = combinedScore, delta = delta)
  return(aFcsResult)
}
#' @export
psea = function(dataMatrix, classMatrix, phenoGrp, statFun = stat.snr, phenoPerms = TRUE, featurePerms = TRUE, nPerms = DFT.PERMS){
  upStats = pseStats(dataMatrix, classMatrix, phenoGrp, statFun)
  nClass = apply(classMatrix > 0 ,2, sum)
  aPseaResult = data.frame(ClassName = colnames(classMatrix), nFeatures = nClass, SetStat = upStats)
  if(phenoPerms){
    aPhenoScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = pseaPhenoPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aPhenoScore = NaN
  }
  #print(featureScoreFun)
  if(featurePerms){
    aFeatureScore = pScore(dataMatrix, classMatrix, phenoGrp, upStats, permFun = pseaSpotPerms, statFun = statFun, nPerms = nPerms)
  } else {
    aFeatureScore = NaN
  }
  aPseaResult = data.frame(aPseaResult, phenoScore = aPhenoScore,
                           featureScore = aFeatureScore,
                           pPhenoScore = -log10(aPhenoScore),
                           pFeatureScore = -log10(aFeatureScore))

  combinedScore = apply(as.matrix(aPseaResult[, c("pPhenoScore", "pFeatureScore")]),1,sum, na.rm = TRUE)
  aPseaResult = data.frame(aPseaResult, combinedScore = combinedScore)
  return(aPseaResult)
}
#' @export
stat.identity = function(X, grp = NULL, pair = NULL){
  # return M as a columns vector, ignore grp
  # this is intended for the case that the dataMatrix is in fact a column vector with stats.
  M = matrix(nrow = length(X), ncol = 1, data = X)
  return(M)
}

#' @export
signOfStat = function(grp){
  # Note: nor for stat.identity
  s = c(-1,1)
  names(s) = levels(grp)
  return(s)
}

#' @export
stat.snr = function(M, grp, pair = NULL){
  # snr statistic per peptide
  bGrp1 = grp == levels(grp)[1]
  bGrp2 = grp == levels(grp)[2]
  m1 = apply(M[bGrp1,], 2, mean)
  s1 = apply(M[bGrp1,], 2, var)
  m2 = apply(M[bGrp2,], 2, mean)
  s2 = apply(M[bGrp2,], 2, var)
  aStat = (m2 - m1) / sqrt(s1+s2)
  # mask any constant columns by 0
  aStat[is.nan(aStat)] = 0
  return(aStat)
}
#' @export
stat.delta = function(M, grp, pair = NULL){
  # difference statistic per peptide
  Mgrp1 = as.matrix(M[grp == levels(grp)[1],])
  Mgrp2 = as.matrix(M[grp == levels(grp)[2],])
  if(dim(Mgrp1)[2] > 1){
    m1 = apply(Mgrp1, 2, mean)
  } else {
    m1 = Mgrp1
  }
  if (dim(Mgrp2)[2] > 1){
    m2 = apply(Mgrp2, 2, mean)
  } else {
    m2 = Mgrp2
  }
  aStat = m2-m1
  return(aStat)
}

#' @export
stat.onesample = function(M, grp = NULL, pair = NULL){
  m0 = apply(M, 2, mean)
  s0 = apply(M, 2, sd)
  aStat = m0/s0
  aStat[s0 == 0]  = 0
  return(aStat)
}

#' @export
stat.paired = function(M, grp, pair){
  d = data.frame(X, grp, pair)
  dm = melt(d, id.vars = c("grp", "pair"))
  Mp = acast(dm,pair ~ variable ~ grp)
  if (dim(Mp)[3]!= 2) stop("Improper pairing")
  return(stat.onesample(Mp[,,2]-Mp[,,1]))

}

mvCorrelationScore = function(dataMatrix, classMatrix, phenoGrp, centerData = TRUE, scaleData = FALSE){
  X0 = scale(dataMatrix, center = centerData, scale = scaleData)
  Xp = X0 %*% classMatrix
  y = matrix(nrow = length(phenoGrp), ncol = 1, data = as.numeric(phenoGrp))
  r = cor(y, Xp)
  mr = r %*% corpcor::invcor.shrink(Xp) %*% t(r)
  return(mr)
}

setStats = function(M, classMat, grp, statfun){
  # statistic per set
  setStat = t(classMat) %*% statfun(M, grp)
  return(setStat)
}

pseStats = function(M, classMat, grp, statfun){
  # enrichment stats
  pepStats = statfun(M, grp)
  idx = order(-pepStats)
  pepStats = pepStats[idx]
  classMat = classMat[idx,]
  pse = apply(classMat,2, pseSingleSet, pepStats)
  return(pse)
}

pseSingleSet = function(class, stats, p = 1){
  Nhit = abs(stats)^p %*% class
  Nmiss = length(class)-sum(class > 0)
  hit  = cumsum( (class* abs(stats)^p)/Nhit)
  miss = cumsum( (!(class > 0))/Nmiss)
  ES = max(abs(hit-miss))
  return(ES)
}

pScore = function(dataMatrix, classMatrix, phenoGrp, upStat, permFun, statFun = stat.snr, nPerms = DFT.PERMS){
  pSetStat = permFun(dataMatrix, classMatrix, grp = phenoGrp, statFun = statFun, nPerms = nPerms)
  nPerms = dim(pSetStat)[2]
  myScore = vector(length = length(upStat))
  for (i in 1:length(upStat)){
    myScore[i] = max(sum(abs(pSetStat[i,])>=abs(upStat[i]))/nPerms, 1/nPerms)
  }
  return(myScore)
}

fcsPhenoPerms = function(dataMatrix, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  grpPerms = mkPerms(grp, nPerms)
  nPerms = dim(grpPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pGrp = grp[grpPerms[,i]]
    pStat[,i] = setStats(dataMatrix, classMatrix, pGrp, statFun)
  }
  return(pStat)
}

pseaPhenoPerms = function(dataMatrix, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  grpPerms = mkPerms(grp, nPerms)
  nPerms = dim(grpPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pGrp = grp[grpPerms[,i]]
    pStat[,i] = pseStats(dataMatrix, classMatrix, pGrp, statFun)
  }
  return(pStat)
}

fcsSpotPerms = function(X, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  if(is.matrix(X)){
    dataMatrix = X
  } else {
    dataMatrix = matrix(nrow = 1, ncol = length(X), data = X)
  }

  spotPerms = mkPerms(1:dim(dataMatrix)[2], nPerms)
  nPerms = dim(spotPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pMatrix  = dataMatrix[, spotPerms[,i]]
    pStat[,i] = setStats(pMatrix, classMatrix, grp, statFun)
  }
  return(pStat)
}

pseaSpotPerms = function(dataMatrix, classMatrix, grp, statFun, nPerms = DFT.PERMS){
  spotPerms = mkPerms(1:dim(dataMatrix)[2], nPerms)
  nPerms = dim(spotPerms)[2]
  pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
  for (i in 1:nPerms){
    pMatrix  = dataMatrix[, spotPerms[,i]]
    pStat[,i] = pseStats(pMatrix, classMatrix, grp, statFun)
  }
  return(pStat)
}

mkPerms = function(grp, maxPerms = DFT.PERMS){
  if ( length(grp) < 7){
    pIdx = permn(length(grp))
    p = matrix(nrow = length(grp), ncol = length(pIdx), data = unlist(pIdx))
    if (length(pIdx) > maxPerms){
      p = p[, sample(1:maxPerms)]
    }
  } else {
    mIdx = matrix(nrow = length(grp), ncol = maxPerms, data = 1:length(grp))
    p = apply(mIdx, 2, sample)
  }
  return(p)
}
