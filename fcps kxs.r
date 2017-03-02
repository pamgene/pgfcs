#source("C:\\Users\\rdwijn.PAMGENE\\Documents\\140-700 Bioinformatics\\pamapps\\components\\FPCS KXS\\operator\\impl\\fpcs kxs.r")
library(Biobase)
checkpack <- function(packagename){
  if (! (packagename %in% installed.packages() ) ){
    install.packages(packagename, repos = "http://cran.rstudio.com/", dependencies = TRUE)
  }
  library(packagename,character.only = TRUE)
}

try(detach(package:reshape), silent = TRUE)
checkpack("plyr")
checkpack("reshape")
checkpack("combinat")
checkpack("ggplot2")
checkpack("corpcor")
checkpack("shiny")
checkpack("scales")

DFT.PERMS = 500

dataFrameOperator <- function(data,properties=properties,folder=folder) {
	# main entry point, handling I/O
	metaData <- varMetadata(data)

	spotLabel 		=  	colnames(pData(data))[metaData$groupingType=='xAxis']
	grpLabel 		=  	colnames(pData(data))[metaData$groupingType == 'Color']
	sampleLabel  	=  	colnames(pData(data))[metaData$groupingType == 'Array']
	classLabel 		= 	getPropertyValue(properties, 'Kinase factor')
	rankLabel 		=   getPropertyValue(properties, 'Kinase rank factor')
	scoreLabel      =   getPropertyValue(properties, 'Score factor')
	
	if( length(grpLabel) != 1){
		stop("Need exactly one factor to define the grouping")
	}
	
	if (metaData$labelDescription[metaData$groupingType == 'Color'] %in% (metaData$labelDescription[metaData$groupingType == 'Array'])){
		stop("The factor used to define the grouping cannot be used to also define the samples")
	}
	

	
	kName = as.factor(pData(data)[[classLabel]])
	kRank = as.numeric(pData(data)[[rankLabel]])
	score = as.numeric(pData(data)[[scoreLabel]])
	
	samples = droplevels(interaction(pData(data)[sampleLabel]))
	aD = data.frame(rowSeq = pData(data)$rowSeq, colSeq = pData(data)$colSeq
			, samples = samples
			, spot = data[[spotLabel]]
			, pClass = kName
			, kRank = kRank 
			, grp = as.factor(data[[grpLabel]])
			, score = score
			, val = data[["value"]]) 

	lvGrp = levels(as.factor(aD$grp))
	if (length(lvGrp)!= 2){
		stop("Need exactly two groups to perform the analysis")
	}
	grp = cast(aD, colSeq~spot, fun.aggregate = function(x) return(x[1]), value = "grp")[,2]
	print(summary(grp))
	
	# parse properties
	MaxSS = getPropertyValue(properties, "Max kinase set size", TRUE)
	
	if (getPropertyValue(properties, "Weight") == " yes"){
		use.weigth = TRUE
	}else{
		use.weight = FALSE
	}
	
	kStatType = getPropertyValue(properties, "Kinase statistic")
	if (kStatType == "Auto"){
		if (min(summary(grp)) <= 2){
			kStatType = "Difference"
		} else {
			kStatType = "SNR"
		}
	}
	statfun = switch(kStatType, SNR = snr, Difference = delta)
	if (is.null(statfun) ){
		stop("Invalid value for property kStatType")
	}
	kScoreMethod = getPropertyValue(properties, "Rank kinases on")
	if (kScoreMethod == "Auto"){
		if ( min(summary(grp)) <= 2 | length(grp) < 8){
			kScoreMethod = "Specificity score only"
		}else{
			kScoreMethod = "Specificity score + Significance score"
		}
	}
	if (kScoreMethod == "Specificity score + Significance score"){
		phenoScoreFun = phenoPerms
	} else {
		phenoScoreFun = NULL
	}

	# scan (run the analysis for selected kinase ranks)
	aScan = getPropertyValue(properties, "Scan rank from", TRUE) : getPropertyValue(properties, "Scan rank to", TRUE) 
	aScanResult = alply(aScan, .margins = 1, .fun = dfScore, aFrame = aD, 
								grp = grp, 
								use.weight = use.weight, 
								statfun = statfun,
								scorefun = phenoScoreFun,
								maxss = MaxSS,
								.progress=progress_win(title="Scanning ...")) 	
	
	# collect data for graph output
	aRunFile = file.path(folder, "runData.RData")
	save(file = aRunFile, aScanResult, lvGrp, grpLabel, kStatType, kScoreMethod)

	#Return optim point data to BN.
	anOptim = ldply(aScanResult, .fun = function(aList)return(data.frame(mx.rank = aList$maxRank, m.cor = aList$mr)))
	optIdx = which.max(anOptim$m.cor)
	scoreFrame = aScanResult[[optIdx]]$aResult
	# map back "kinaseName" to "classLabel" in pdata
	aResult = ddply(scoreFrame, ~kinaseName, .fun = mapback, backFrame = aD)[,-1]
	aResult[is.na(aResult)] = NaN
	print("Upstream kinase analysis completed")
	print("The folowing settings were used:")
	print(paste("Kinase statistic:",kStatType))
	print(paste("Scoring method:", kScoreMethod))
	return (aResult)
}

rankspots = function(aFrame){
  # do the ranking for 1 column of the ss, then expand to the other columns
  aFrame = aFrame[order(aFrame$colSeq, aFrame$rowSeq),]
  aSub = subset(aFrame, colSeq == 1)
  aRank = rank(-aSub$score, ties.method = "min")
  aFrame = data.frame(aFrame, spot.rank = aRank)
  return(aFrame)
}

dfScore = function(maxRank, aFrame, grp, use.weight = FALSE, statfun = snrStatistic, scorefun = phenoPerms, maxss = 100){
	# run the analysis for selected kinase rank

	aFrame = subset(aFrame, kRank <= maxRank)
	aFrame = ddply(aFrame, ~pClass, .fun = rankspots)
	aFrame = subset(aFrame, spot.rank <= maxss)
	aClassMatrix = cast(aFrame, spot ~ pClass, fun.aggregate = mean, value = "score")
	rownames(aClassMatrix) = aClassMatrix[,1]
	if(!use.weight){
		aClassMatrix = !is.na(aClassMatrix[,2:ncol(aClassMatrix)])
	} else {
		aClassMatrix = aClassMatrix[,2:ncol(aClassMatrix)]
		aClassMatrix[is.na(aClassMatrix)] = 0
	}
	
	aDataMatrix  = cast(aFrame, colSeq~spot, fun.aggregate = mean, value = "val")  
	
	aDataMatrix = as.matrix(aDataMatrix[, 2:dim(aDataMatrix)[2]])
	if(any(is.na(aDataMatrix))){
		stop("Missing values are not allowed.")
	}
	
	# get the multiple correlation using the kinase projection of the data
	X0 = scale(aDataMatrix, center = TRUE, scale = FALSE)
	Xp = X0 %*% aClassMatrix 
	y = matrix(nrow = length(grp), ncol = 1, data = as.numeric(grp))
	r = cor(y, Xp)
	mr = r %*% invcor.shrink(Xp) %*% t(r)
	# and do the main analysis
	aSensResult= phenoScore(aDataMatrix, aClassMatrix, grp, statfun = statfun, scorefun = scorefun)
	spcScore = specScore(dataMatrix = aDataMatrix, classMatrix = aClassMatrix, grp = grp, setStat = aSensResult$setStat, statfun = statfun, nPerms = DFT.PERMS)
	nClass = apply(aClassMatrix > 0 ,2, sum)
	aResult = data.frame(	kRank = maxRank, 
							kinaseName = colnames(aClassMatrix), 
							sensitivityScore = aSensResult$setScore, 
							specificityScore = spcScore, 
							setStat = aSensResult$setStat,
							Peptides.in.Set = nClass,
							Normalized.Kinase.Stat = aSensResult$setStat/nClass,
							pSpecificityScore = -log10(spcScore),
							pSensitivityScore = -log10(aSensResult$setScore) )
							
							
	if (!is.null(scorefun)){
		aResult = data.frame(aResult, finalScore = aResult$pSpecificityScore + aResult$pSensitivityScore)
	} else {
		aResult = data.frame(aResult, finalScore = aResult$pSpecificityScore)
	}
	aResultList = list(maxRank = maxRank, aResult = aResult, classMatrix = aClassMatrix, dataMatrix = aDataMatrix, grp = grp, mr = mr)
	return(aResultList)
}

mapback = function(kResult, backFrame){
	# I/O helper function
	bMap = backFrame$pClass %in% kResult$kinaseName
	aResult = unique(data.frame(rowSeq = backFrame$rowSeq[bMap], colSeq = backFrame$colSeq[bMap]))
	aResult = data.frame(aResult,	pSensitivityScore = kResult$pSensitivityScore,
									pSpecificityScore = kResult$pSpecificityScore,
									sensitivityScore = kResult$sensitivityScore, 
									specificityScore = kResult$specificityScore,
									finalScore = kResult$finalScore,
									setStat = kResult$setStat,
									Normalized.Kinase.Stat = kResult$Normalized.Kinase.Stat)
	return(aResult)
}

   
 operatorProperties <- function() {
		propList = list(  	
						
							list('Kinase statistic', list('Auto', 'SNR', 'Difference')),
							list('Rank kinases on', list('Auto', 'Specificity score + Significance score', 'Specificity score only')),
							list('Scan rank from', 5), 
							list('Scan rank to', 15),
							list('Max kinase set size', 50),
							list('Kinase factor', 'Kinase.name'),
							list('Kinase rank factor', 'Kinase.Rank'),
							list('Weight', list('none','score') ),
							list('Score factor', 'Kinexus.Predictor.Version.2.Score') 
						)    
    return (propList)
}

getPropertyValue <- function(properties=properties, name=name, prop.is.numeric = FALSE) {
	for ( i in 1:length(properties) ) {
		if (properties[[i]][1] == name){
			val = properties[[i]][2]
			if (prop.is.numeric){
				return(as.numeric(val))
			} else {
				return(val)
			}
		}
	}
	stop(paste("property not found: ", name))
	return (NULL)
}

snr = function(M, grp){
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

delta = function(M, grp){
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

setStats = function(M, classMat, grp, statfun){
	# snr statistic per set
	setStat = t(classMat) %*% statfun(M, grp)
	return(setStat)
}

phenoScore = function(dataMatrix, classMatrix, grp, nPerms = DFT.PERMS, statfun = snr, scorefun = phenoPerms){	
	mySetStat = setStats(dataMatrix, classMatrix, grp, statfun = statfun)
	if (!is.null(scorefun)){
		mySetScore = scorefun(dataMatrix, classMatrix, grp = grp, setStat = mySetStat, statfun = statfun, nPerms = nPerms)
	} else {
		mySetScore = NA
	}
	aResult = list(setStat = mySetStat, setScore = mySetScore)
	return(aResult)
}

specScore = function(dataMatrix, classMatrix, grp, setStat, statfun = snr, nPerms = DFT.PERMS){
	pSetStat = spotPerms(dataMatrix, classMatrix, grp, statfun, nPerms)
	nPerms = dim(pSetStat)[2]
	myUniscore = vector(length = length(setStat))
	for (i in 1:length(setStat)){
		myUniscore[i] = max(sum(abs(pSetStat[i,])>=abs(setStat[i]))/nPerms, 1/nPerms) 
	}
	return(myUniscore)
}

phenoPerms = function(dataMatrix, classMatrix, grp, setStat, statfun, nPerms = DFT.PERMS){
	grpPerms = mkPerms(grp, nPerms)
	nPerms = dim(grpPerms)[2]
	pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
	for (i in 1:nPerms){
		pGrp = grp[grpPerms[,i]]
		pStat[,i] = setStats(dataMatrix, classMatrix, pGrp, statfun)
	}	
	setScore = vector(length = length(setStat))
	for (i in 1:length(setScore)){
			setScore[i] = max(sum(abs(pStat[i,])>=abs(setStat[i]))/nPerms, 1/nPerms) 
	}
	return(setScore)
}

spotPerms = function(dataMatrix, classMatrix, grp, statfun, nPerms = DFT.PERMS){
	spotPerms = mkPerms(1:dim(dataMatrix)[2], nPerms)
	nPerms = dim(spotPerms)[2]
	pStat = matrix(nrow = dim(classMatrix)[2], ncol = nPerms)
	for (i in 1:nPerms){
		pMatrix  = dataMatrix[, spotPerms[,i]]
		pStat[,i] = setStats(pMatrix, classMatrix, grp, statfun)
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

showResults <- function(properties=properties, folder=folder){
	aRunFile = file.path(folder, "runData.RData")
	load(aRunFile)
	statOption 	= kStatType
	scoreOption = kScoreMethod
	idxTest = grep(statOption, c('SNR', 'Difference'))
	if (length(idxTest) == 0){
		stop("Invalid value for parameter Kinase statistic")
	}
	if (idxTest == 1){
		statfun = snr
	}
	if (idxTest == 2){
		statfun = delta
	}
	
	# optim plot
	if (length(aScanResult)> 1){
		anOptim = ldply(aScanResult, .fun = function(aList)return(data.frame(mx.rank = aList$maxRank, m.cor = aList$mr)))
		op = ggplot(anOptim, aes(x = mx.rank, y = m.cor)) + geom_point(colour = "blue") + geom_line()
		op = op + xlab("Maximum substrate rank") + ylab("Correlation with grouping") + scale_x_continuous(breaks=anOptim$mx.rank[1]:anOptim$mx.rank[length(anOptim$mx.rank)])
	}

	
	# score plot
	optIdx = which.max(anOptim$m.cor)
	scoreFrame = aScanResult[[optIdx]]$aResult
	scoreFrame$clrScore = scoreFrame$pSpecificityScore # clr for graph output
	scoreFrame$clrScore[scoreFrame$clrScore > 2] = 2
	scoreFrame$clrScore = 0.5 * scoreFrame$clrScore
	sp = ggplot(scoreFrame) 
	sp = sp + geom_bar(aes(x = reorder(kinaseName, finalScore), y = Normalized.Kinase.Stat, fill = clrScore), stat = "identity")

	sp = sp + scale_fill_gradientn(name = "Specificity Score",space = "rgb", 
                                       colours = c("black", "white", "red"), 
                                       values = c(0,0.65, 1),
                                       breaks = c(0,0.65, 1),
                                       labels = c(0,0.65, 1)*2)
									   
	sp = sp + xlab("Kinase Top List") + ylab("Normalized kinase statistic")
	sp = sp + coord_flip()
	sp = sp + theme(axis.text.y = element_text(size = 10)) 
	
	# kinase volcano
	kv = ggplot(scoreFrame, aes(x = Normalized.Kinase.Stat, y = finalScore, label = kinaseName, size = Peptides.in.Set, colour = clrScore))
	kv = kv + scale_colour_gradientn(name = "Specificity Score",space = "rgb", 
                                       colours = c("black", "white", "red"), 
                                       values = c(0,0.65, 1),
                                       breaks = c(0,0.65, 1),
                                       labels = c(0,0.65, 1)*2)
	kv = kv + theme(panel.background = element_rect(fill = 'lightblue', colour = 'black'))
	kv = kv + geom_text() + scale_size(range = c(3, 10)) + xlab("Normalized kinase statistic") + ylab("Overall score")
	
	# kinase basis plot
	perSpotStats  = statfun(aScanResult[[optIdx]]$dataMatrix, aScanResult[[optIdx]]$grp)
	aStat = matrix(nrow = nrow(aScanResult[[optIdx]]$classMatrix), ncol = ncol(aScanResult[[optIdx]]$classMatrix), data = perSpotStats^2)
	aStat[!(aScanResult[[optIdx]]$classMatrix>0)] = NA
	colnames(aStat) = colnames(aScanResult[[optIdx]]$classMatrix)
	rownames(aStat) = rownames(aScanResult[[optIdx]]$classMatrix)
	nPeptides = dim(aStat)[1]
	nKinases = dim(aStat)[2]
	aStatBasis = melt(aStat)
	colnames(aStatBasis) = c("Peptide", "Kinase.Name","Stat.Squared")
	class.order = t(matrix(nrow = dim(aStat)[2], ncol = dim(aStat)[1], data = scoreFrame$finalScore))
	spot.order = matrix(nrow = dim(aStat)[1], ncol = dim(aStat)[2], data = perSpotStats^2) 
	aStatBasis = data.frame(aStatBasis, class.order = as.vector(class.order), spot.order = as.vector(spot.order))
	sb = ggplot(aStatBasis, aes(x = reorder(Peptide, spot.order), y = reorder(Kinase.Name, class.order) )) + geom_tile(aes(fill = Stat.Squared)) 
	sb = sb + labs(x = "Peptide", y = "Kinase")
	sb = sb + theme(axis.text.y = element_text(size = 6) , axis.text.x = element_text(size = 6, angle = 45, hjust = 1))
	sb = sb + scale_fill_gradient(high = "green")

	# advanced plot
	if(length(aScanResult) > 1){
		scanFrame = ldply(aScanResult, function(aList)return(aList$aResult))
		# format a df for creating the rank history plot in ggplot
		scanFrame$pSpecificityScore = -log10(scanFrame$specificityScore)
		scanFrame$pSensitivityScore = -log10(scanFrame$sensitivityScore)
		scanFrame = ddply(scanFrame, ~kRank, .fun = rankByScore) 
		cntFrame = ddply(scanFrame, ~kinaseName, .fun = rcount, scoreFrame, clip = 20)

		rPlot = ggplot(cntFrame, aes(x = reorder(kinaseName, finalScore), y = rank)) + geom_point(aes(size = count, colour = count)) + coord_flip()
		rPlot = rPlot + scale_y_continuous(limits = c(0.6, 20.4), breaks = 1:20) + ylab("Relative Score") + xlab("Kinase Name") 
		rPlot = rPlot + scale_size(range = c(2, 10), breaks = c(1,2,4,8))
	}
	

	scrPlotHelpTxt = paste("This plot shows kinases ranked according to their score in the upstream kinase analysis, where the difference between the",  
							grpLabel,"=", lvGrp[1],"and",grpLabel,"=",lvGrp[2], "samples is analyzed using sets of peptides that are predicted to be substrates of the respective kinases.", 
							"\nThe length of the bars shows the Kinase statistic, a measure for the change between the groups (",statOption,").",
							"\nA positive value of the kinase statistic indicates that the associated kinase activity was higher in the ", grpLabel, "=", lvGrp[2], " group." ,
							"\nThe color of the bars indicates the specificity of the kinase set",
							"\nThe kinases are ranked based on ", scoreOption, ".")
							
	basisPlotHelpTxt = paste("Peptide to kinase matrix used as the basis for the upstream kinase scoring. The green tiles indicate peptides that are included in the set for the respective kinases",
							 "The brightness of the marker indicates a signal to noise value that is calculated for each peptide individually for the difference between the",
							 grpLabel,"=", lvGrp[1],"and",grpLabel,"=",lvGrp[2], "samples. The matrix contains",
							 nPeptides,"different peptides that are substrates for one or more of",nKinases,"kinases." 
							 )
	
	optimPlotHelpTxt = paste(	"For constructing the peptide to kinase matrix for the upstream kinase scoring it is necessary to decide how many putative kinases to include for each substrate",
								"For this, the data separation between the",
								grpLabel,"=", lvGrp[1],"and",grpLabel,"=",lvGrp[2],
								"was analyzed, after including kinases that are in the top", anOptim$mx.rank[1],
								"up to including the top", anOptim$mx.rank[length(aScanResult)],".",
								"Using the top", anOptim$mx.rank[optIdx], "resulted in optimal data separation, and this setting was used for constructing the Upstream Kinase Bassis and calculating the  Upstream Kinase Score.") 
								
	scanPlotHelpTxt = paste(	"The plot below shows how the relative scoring (rank) of kinases as a function of the number of kinases included for each peptide.", 
								"This can be used to check if the final result depends strongly on the analysis settings.",
								"Each point indicates what rank was attained how often (size) during subsequent analyses. In this graph a rank of twenty indicates twenty or higher.") 
								
	volcanoPlotHelpTxt = paste(	"The Kinase Volcano Plot shows for each kinase on the x-axis the Kinase statistic (", statOption,") , a measure for the difference between ",grpLabel,"=",lvGrp[1], " and ",grpLabel,"=",lvGrp[2] 
								,", and on the y-axis the score used for ranking the kinases: ", scoreOption
								,". \nA positive value of the kinase statistic indicates that the associated kinase activity was higher in the ", grpLabel, "=", lvGrp[2], " group." 
								, sep = "")				
	
	aShinyList = list(ui = 
		pageWithSidebar(
			# Application title
			headerPanel("Upstream Kinase Analysis"),
			# Sidebar with controls 
			sidebarPanel(
				helpText("Putative upstream kinase scoring using kinase-substrate predictions."),
				h5("For more information:", a("help on PamCloud", href="https://pamcloud.pamgene.com/wiki/Wiki.jsp?page=Background%20to%20the%20Upstream%20Kinase%20Analysis%20PamApp") ),
				tags$hr(),
				h5(Color = "DarkRed", "Make sure to press the Done button before returning to Bionavigator!"),
				actionButton("stop", "Done")
			),

  #	 mainpanel
			mainPanel( 
				tabsetPanel(
					tabPanel( "Upstream Kinase Score",helpText(scrPlotHelpTxt), plotOutput("scorePlot", height = nKinases *12)),
					tabPanel( "Kinase Volcano", helpText(volcanoPlotHelpTxt), plotOutput("volcanoPlot", height = 700)),
					tabPanel("Upstream Kinase Basis", helpText(basisPlotHelpTxt), plotOutput("basisPlot", height = max(700,nKinases*9))),
					tabPanel("Advanced", helpText(optimPlotHelpTxt), plotOutput("optimPlot", width = 600), helpText(scanPlotHelpTxt), plotOutput("scanResults", height = nKinases*20))
				),
				textOutput("done")
			)		
		)
	,
	server = function(input, output, session) {
		session$onSessionEnded(function() {
			stopApp();
		})	

		output$scorePlot <- renderPlot({
			print(sp)
		})
		
		output$volcanoPlot <- renderPlot({
			print(kv)
		})
		
		output$basisPlot <- renderPlot({
			print(sb)
		})
		
		output$optimPlot <- renderPlot({
			print(op)
		})
		
		output$scanResults <- renderPlot({
			print(rPlot)
		})
		
		done <- reactive({return(input$stop>0)})
		
		output$done = renderText({
			if(done()){
				stopApp()
				return("Done")
			} else {
				return("")
			}
		})
	})	

	aTry = try(runApp(aShinyList), silent = FALSE)
	if(inherits(aTry, "try-error")){
		#handle the case of an already running shiny app
		stopApp()
		tkmessageBox(message = "Could not start the Shiny service, please retry")
	}

}

rankByScore = function(aFrame){
	# showResults helper function
	aFrame = data.frame(aFrame, setRank = rank( -aFrame$finalScore , ties.method = "max"))
	return(aFrame)
}

rcount = function(aFrame, scores, clip = 20){
  # show results helper function
  aFrame$setRank[aFrame$setRank>clip] = clip
  rankAsFactor = as.factor(aFrame$setRank)
  rankSummary = summary(rankAsFactor)

  kinaseMatch = scores$kinaseName %in% droplevels(aFrame$kinaseName)

  if(any(kinaseMatch)){
	combScore = scores$finalScore[kinaseMatch]
  } else {
	combScore = 0
  }
  aResult = data.frame(rank = as.numeric(names(rankSummary)), count = rankSummary, finalScore = combScore)
  return(aResult)
}