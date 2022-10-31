plotDEXSeq_RDS <- function (object, geneID, file_path, FDR = 0.1, fitExpToVar = "condition", 
    norCounts = FALSE, expression = TRUE, splicing = FALSE, displayTranscripts = FALSE, 
    names = FALSE, legend = FALSE, color = NULL, color.samples = NULL, 
    transcriptDb = NULL, additionalAnnotation = NULL, maxRowsMF = 2400, 
    ...) 
{
    stopifnot(is(object, "DEXSeqResults") | is(object, "DEXSeqDataSet"))
    if (!fitExpToVar %in% colnames(object@modelFrameBM)) {
        stop(sprintf("The value of the parameter fitExpToVar,'%s', is not a column name of the 'colData' DataFrame from the DEXSeqDataSet object.", 
            fitExpToVar))
    }
    op <- sum(c(expression, splicing, norCounts))
    if (op == 0) {
        stop("Please indicate what would you like to plot\n")
    }
    if (!is.null(transcriptDb)) {
        stopifnot(is(transcriptDb, "TxDb"))
    }
    if (!is.null(additionalAnnotation)) {
        stopifnot(is(additionalAnnotation, "GRangesList"))
    }
    if (is(object, "DEXSeqResults")) {
        sampleData <- object@sampleData
        genomicData <- object$genomicData
        rt <- which(object$groupID == geneID)
        count <- t(t(object$countData[rt, ])/sampleData$sizeFactor)
        transcripts <- object$transcripts[rt]
        each <- object$padj[rt]
    }
    else {
        sampleData <- sampleAnnotation(object)
        genomicData <- rowRanges(object)
        rt <- which(mcols(object)$groupID == geneID)
        transcripts <- genomicData$transcripts[rt]
        mcols(genomicData) <- NULL
        count <- featureCounts(object, normalized = TRUE)[rt, 
            ]
        each <- rep(1, length.out = length(rt))
    }
    if (sum(count) == 0) {
        warning("No read counts falling in this gene, there is nothing to plot.")
        return()
    }
    if (FDR > 1 | FDR < 0) {
        stop("FDR has to be a numeric value between 0 - 1")
    }
    rango <- seq(along = rt)
    intervals <- (0:nrow(count))/nrow(count)
    numcond <- length(unique(sampleData[[fitExpToVar]]))
    numexons <- nrow(count)
    exoncol <- ifelse(each <= FDR, "#F219ED", "#CCCCCC")
    exoncol[is.na(exoncol)] <- "white"
    colorlines <- ifelse(each <= FDR, "#F219ED60", "#B3B3B360")
    colorlines[is.na(colorlines)] <- "#B3B3B360"
    colorlinesB <- ifelse(each <= FDR, "#9E109B", "#666666")
    colorlinesB[is.na(colorlinesB)] <- "#666666"
    if (length(unlist(start(genomicData))) > 0) {
        sub <- data.frame(start = start(genomicData[rt]), end = end(genomicData[rt]), 
            chr = as.character(seqnames(genomicData[rt])), strand = as.character(strand(genomicData[rt])))
        rownames(sub) <- rownames(object)[rt]
        if (!is.null(additionalAnnotation)) {
            additionalHits <- findOverlaps(additionalAnnotation, 
                range(genomicData[rt]))
            additionalAnnotation <- additionalAnnotation[queryHits(additionalHits)]
            if (length(additionalAnnotation) == 0) {
                additionalAnnotation <- NULL
            }
        }
        rel <- (data.frame(sub$start, sub$end)) - min(sub$start)
        rel <- rel/max(rel[, 2])
        trans <- unique(unlist(transcripts))
        trans <- trans[!is.na(trans)]
        numberOfTrans <- length(trans) + length(additionalAnnotation)
        if ((displayTranscripts & !is.null(unlist(transcripts))) | 
            !is.null(additionalAnnotation)) {
            if (numberOfTrans > 40) {
                warning("This gene contains more than 40 transcripts annotated, only the first 40 will be plotted\n")
            }
            if (!displayTranscripts) {
                numberOfTrans <- numberOfTrans - length(trans)
            }
            mat <- seq_len(3 + min(numberOfTrans, 40))
            hei <- c(8, 1, 1.5, rep(1.5, min(numberOfTrans, 40)))
        }
        else {
            mat <- 1:3
            hei <- c(5, 1, 1.5)
        }
        if (op > 1) {
            hei <- c(rep(hei[1], op - 1), hei)
            mat <- c(mat, length(mat) + seq(along = op))
        }
        hei <- c(hei, 0.2)
        mat <- c(mat, length(mat) + 1)
        layout(matrix(mat), heights = hei)
        par(mar = c(2, 4, 4, 2))
    }
    else if (op > 1) {
        par(mfrow = c(op, 1))
    }
    if (is.null(color)) {
        if (numcond < 10) {
            color <- suppressWarnings(brewer.pal(numcond, "Set1")[seq_len(numcond)])
        }
        else {
            color <- rgb(colorRamp(brewer.pal(5, "Set1"))(seq(0, 
                1, length.out = numcond)), maxColorValue = 255, 
                alpha = 175)
        }
    }
    names(color) <- sort(unique(as.character((sampleData[[fitExpToVar]]))))
    if (expression | splicing) {
        if (!is(object, "DEXSeqResults")) {
            stop("To visualize beta estimates, please provide a DEXSeqResults object")
        }
        effects <- getEffectsForGene(geneID, object, maxRowsMF, 
            fitExpToVar)
        if (is.null(effects[["splicing"]])) {
            return()
        }
        if (!all(rownames(sub) %in% rownames(effects[["splicing"]]))) {
            return()
        }
    }
    if (expression) {
        coeff <- effects[["expression"]]
        coeff <- coeff[rownames(sub), ]
        coeff <- exp(coeff)
        ylimn <- c(0, max(coeff, na.rm = TRUE))
        coeff <- vst(coeff, object)
        drawPlot(matr = coeff, ylimn, object, intervals, rango, 
            textAxis = "Expression", rt = rt, color = rep(color[colnames(coeff)], 
                each = numexons), colorlines = colorlines, ...)
    }
    if (splicing) {
        coeff <- effects[["splicing"]]
        coeff <- coeff[rownames(sub), ]
        coeff <- exp(coeff)
        ylimn <- c(0, max(coeff, na.rm = TRUE))
        coeff <- vst(coeff, object)
        plot_object = list(coeff, ylimn, object, intervals, rango, rt, numexons, colorlines)
        names(plot_object) = c("coeff", "ylimn", "object", "intervals", "rango", "rt", "numexons", "colorlines")
        saveRDS(plot_object, file_path)
        drawPlot(matr = coeff, ylimn, object, intervals, rango, 
            textAxis = "Exon usage", rt = rt, color = rep(color[colnames(coeff)], 
                each = numexons), colorlines = colorlines, ...)
    }
    if (norCounts) {
        ylimn <- c(0, max(count, na.rm = TRUE))
        count <- vst(count, object)
        if (is.null(color.samples)) {
            colorcounts <- rep(color[as.character(sampleData[[fitExpToVar]])], 
                each = numexons)
        }
        else {
            colorcounts <- rep(color.samples, each = numexons)
        }
        drawPlot(matr = count, ylimn, object, intervals, rango, 
            textAxis = "Normalized counts", rt = rt, color = colorcounts, 
            colorlines = colorlines, ...)
    }
    if (length(unlist(start(genomicData))) > 0) {
        par(mar = c(0, 4, 0, 2))
        plot.new()
        segments(apply((rbind(rel[rango, 2], rel[rango, 1])), 
            2, median), 0, apply(rbind(intervals[rango], intervals[rango + 
            1] - ((intervals[rango + 1] - intervals[rango]) * 
            0.2)), 2, median), 1, col = colorlinesB)
        par(mar = c(1.5, 4, 0, 2))
        drawGene(min(sub$start), max(sub$end), tr = sub, exoncol = exoncol, 
            names, trName = "Gene model", cex = 0.8)
        if (length(unlist(transcripts)) > 0) {
            i <- 1
            if (displayTranscripts) {
                for (i in seq_len(min(length(trans), 40))) {
                  logicexons <- sapply(transcripts, function(x) {
                    length(which(x == trans[i]))
                  })
                  tr <- reduce(IRanges(sub$start[logicexons == 
                    1], sub$end[logicexons == 1]))
                  if (is.null(transcriptDb)) {
                    tr <- as.data.frame(tr)[, c("start", "end")]
                    drawGene(min(sub$start), max(sub$end), tr = tr, 
                      exoncol = "black", names, trName = trans[i], 
                      cex = 0.8)
                  }
                  else {
                    codingRanges <- select(transcriptDb, keys = trans[i], 
                      columns = c("CDSSTART", "CDSEND"), keytype = "TXNAME")
                    if (is.na(any(codingRanges$CDSSTART))) {
                      tr <- as.data.frame(tr)[, c("start", "end")]
                      drawGene(min(sub$start), max(sub$end), 
                        tr = tr, exoncol = NULL, names, trName = trans[i], 
                        cex = 0.8, miny = 0.25, maxy = 0.75)
                    }
                    else {
                      codingRanges <- IRanges(codingRanges$CDSSTART, 
                        codingRanges$CDSEND)
                      utrRanges <- setdiff(tr, codingRanges)
                      drawGene(min(sub$start), max(sub$end), 
                        tr = as.data.frame(codingRanges)[, c("start", 
                          "end")], exoncol = "black", names, 
                        trName = trans[i], cex = 0.8, drawNames = FALSE, 
                        drawIntronLines = FALSE)
                      if (length(utrRanges) > 0) {
                        drawGene(min(sub$start), max(sub$end), 
                          tr = as.data.frame(utrRanges)[, c("start", 
                            "end")], exoncol = NULL, names, trName = trans[i], 
                          cex = 0.8, drawNames = FALSE, drawIntronLines = FALSE, 
                          newPanel = FALSE, miny = 0.25, maxy = 0.75)
                      }
                      drawGene(min(sub$start), max(sub$end), 
                        tr = as.data.frame(tr)[, c("start", "end")], 
                        exoncol = "black", names, trName = trans[i], 
                        cex = 0.8, newPanel = FALSE, drawExons = FALSE)
                    }
                  }
                }
            }
            if (!is.null(additionalAnnotation)) {
                for (j in seq_along(additionalAnnotation)) {
                  tr <- as.data.frame(additionalAnnotation[[j]])[, 
                    c("start", "end")]
                  drawGene(min(sub$start), max(sub$end), tr = tr, 
                    exoncol = "darkred", names, trName = names(additionalAnnotation)[j], 
                    cex = 0.8, introncol = "darkred")
                  i <- i + 1
                  if (i > 40) 
                    break
                }
            }
        }
        axis(1, at = round(seq(min(sub$start), max(sub$end), 
            length.out = 10)), labels = round(seq(min(sub$start), 
            max(sub$end), length.out = 10)), pos = 0, lwd.ticks = 0.2, 
            padj = -0.7, ...)
    }
    if (legend) {
        mtext(paste(geneID, unique(sub$strand)), side = 3, adj = 0.25, 
            padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
        posforlegend <- seq(0.7, 0.9, length.out = numcond)
        for (i in seq(along = color)) {
            mtext(names(color[i]), side = 3, adj = posforlegend[i], 
                padj = 1.5, line = 0, outer = TRUE, col = color[i], 
                ...)
        }
    }
    else {
        mtext(paste(geneID, unique(sub$strand)), side = 3, adj = 0.5, 
            padj = 1.5, line = 0, outer = TRUE, cex = 1.5)
    }
}

arrangeCoefs <- function( frm, mf, mm = model.matrix( frm, mf ), fit = NULL, insertValues = TRUE ) {

   if( any( attr( mm, "contrasts" ) != "contr.treatment" ) )
      stop( "Can only deal with standard 'treatment' contrasts." )   # Do I need this?
   if( is.null(fit) & insertValues )
      stop( "If fit==NULL, returnCoefValues must be FALSE" )
   if( !is.null(fit) )
      stopifnot( all( colnames(mm) == names(coefficients(fit)) ) )

   fctTbl <- attr( terms(frm), "factors" )

   coefIndicesList <- 
   lapply( seq_len(ncol(fctTbl)), function( fctTblCol ) {
      termName <- colnames(fctTbl)[ fctTblCol ]
      varsInTerm <- stringr::str_split( termName, stringr::fixed(":") )[[1]] 
      stopifnot( all( fctTbl[ varsInTerm, fctTblCol ] == 1 ) )
      stopifnot( sum( fctTbl[ , fctTblCol ] ) == length( varsInTerm ) )
      coefNames <- colnames(mm)[ attr( mm, "assign" ) == fctTblCol ]
      lvlTbl <- stringr::str_match( coefNames, 
         stringr::str_c( "^", stringr::str_c( varsInTerm, "([^:]*)", collapse=":" ), "$" ) )[ , -1, drop=FALSE ]
      stopifnot( ncol(lvlTbl) == length( varsInTerm ) )
      stopifnot( nrow(lvlTbl) == length( coefNames ) )
      if( !all( sapply( varsInTerm, function(v) is.factor(mf[[v]]) | is.character(mf[[v]]) ) ) )
         stop( "Non-factor in model frame" )

      varLevels <- lapply( varsInTerm, function(v) levels( factor( mf[[v]] ) ) ) 
      coefIndices <- array( NA_character_, dim = sapply( varLevels, length ), dimnames = varLevels )
      names( dimnames( coefIndices ) ) <- varsInTerm

      for( i in seq_len( nrow(lvlTbl) ) )
         coefIndices <- do.call( `[[<-`, c( quote(coefIndices), as.list( lvlTbl[ i, ] ), coefNames[i] ) )

      coefIndices
   } )
   names( coefIndicesList ) <- colnames( fctTbl )

   if( attr( terms(frm), "intercept" ) ) {
      a <- array( c( `(Intercept)` = "(Intercept)" ) )
      dimnames(a) <- list( `(Intercept)` = c( "(Intercept)" ) )
      coefIndicesList <- c( list( `(Intercept)` = a ), coefIndicesList )
   }

   if( !insertValues )
      ans <- coefIndicesList
   else
      ans <- lapply( coefIndicesList, function(coefIndices) {
         a <- ifelse( is.na(coefIndices), 0, coefficients(fit)[ coefIndices ] )
         attr( a, "variables" ) <- attr( coefIndices, "variables" )
         a } )
      
   lapply( ans, function(x) 
      if( is.array(x) ) 
         x 
      else { 
         y <- array( x, dim=length(x) )
         attr( y, "variables" ) <- attr( x, "variables" )
         dimnames(y) <- list( names(x) )
         y } )
}

apply2 <- function( X, MARGIN, FUN, ... ) {
   if( length(MARGIN) > 0 ) 
      apply( X, MARGIN, FUN, ... ) 
   else 
      FUN( X, ... ) }

balanceExons <- function( coefs, dispersions ) {
   stopifnot( any( sapply( coefs, function(x) 
      identical( names(dimnames(x)), "(Intercept)" ) ) ) )
   termsWithExon <- sapply( coefs, function(x) "exon" %in% names(dimnames(x)) )
   meanMainEffect <- sum( sapply( coefs[!termsWithExon], mean, na.rm=TRUE ) )
   meanExonEffects <- rowSums( sapply( coefs[termsWithExon], function(x) 
      apply2( x, "exon", mean, na.rm=TRUE ) ) )

   meanExonFittedValue <- exp( meanMainEffect + meanExonEffects )

   exonWeights <-  1 / ( dispersions + 1 / meanExonFittedValue )

   shifts <- lapply( coefs[termsWithExon], function(x) { 
      nonExonDims <- which(  names(dimnames(x)) != "exon" )
      list(
         vars = names(dimnames(x))[ nonExonDims ],
         wmeans = apply2( x, nonExonDims, weighted.mean, exonWeights) ) } )

   lapply( coefs, function(x) {
      nonExonVars <- names(dimnames(x))[ names(dimnames(x)) != "exon" ]
      if( identical( nonExonVars, "(Intercept)" ) )
         whichShift <- which( sapply( shifts, function(xx) length( xx$vars ) == 0 ) )
      else
         whichShift <- which( sapply( shifts, function(xx) identical( xx$vars, nonExonVars ) ) )
      if( length( whichShift ) == 0 )
         return( x )
      if( length( whichShift ) > 1 )
         stop( "Confused about selecting shift." )
      if( "exon" %in% names(dimnames(x)) )
         x - shifts[[ whichShift ]]$wmeans
      else
         x + shifts[[ whichShift ]]$wmeans
    } )
}         


fitAndArrangeCoefs <- function( frm = count ~ condition * exon, balanceExons = TRUE, mf, fitExpToVar, geneID)
{
   if( length(levels(mf$exon)) <= 1 )
      return( NULL )
   mm <- model.matrix( frm, mf )
   fit <- tryCatch(
       {
           glmnb.fit(mm, mf$count, dispersion=mf$dispersion,
                     offset=log(mf$sizeFactor), tol=0.01)
       },
       error=function(cond){
           message( sprintf("Fit for gene/exon %s failed, coefficients for this gene won't show up.", geneID ) )
           return(NULL)
       },
       warning=function(cond){
           message( sprintf("Fit for gene/exon %s threw the next warning(s): %s", geneID, unique(cond$message) ) )
       } )
   
   if( is.null( fit ) ){
       return( NULL )
   }

   if( is( mf[[fitExpToVar]], "numeric" ) ){
       coefs <- fit$coefficients
       attributes(coefs)$fitType <- "numeric"
   }else{
       coefs <- arrangeCoefs( frm, mf, mm, fit )
       if( balanceExons ) {
           coefs <- balanceExons( coefs, tapply( mf$dispersion, mf$exon, `[`, 1 ) )
       }
       attributes(coefs)$fitType <- "factor"
   }
   coefs
}


getEffectsForPlotting <- function( coefs, groupingVar = "condition", averageOutExpression=FALSE, frm, mf )
{
    if( attributes(coefs)$fitType == "factor" ){
        groupingExonInteraction <- which( sapply( coefs, function(x) 
            all( c( groupingVar, "exon") %in% names(dimnames(x)) ) & length(dim(x)) == 2 ) ) 
        fittedValues <- coefs[[ groupingExonInteraction ]]
        if( names(dimnames(fittedValues))[1] == "exon" )
            fittedValues <- t( fittedValues )
        stopifnot( identical( names(dimnames(fittedValues)), c( groupingVar, "exon" ) ) )
        for( x in coefs[ -groupingExonInteraction ] ) {
            if( all( c( groupingVar, "exon") %in% names(dimnames(x)) ) )
                stop( "Cannot yet deal with third-order terms." )
            if( !any( c( groupingVar, "exon") %in% names(dimnames(x)) ) ) {
                fittedValues <- fittedValues + mean( x )
            } else if( averageOutExpression & identical( names(dimnames(x)), groupingVar ) ) {
                fittedValues <- fittedValues + mean( x )
            } else if( groupingVar %in% names(dimnames(x)) ) {
                groupMeans <- apply2( x, groupingVar, mean )
                stopifnot( identical( names(groupMeans), dimnames(fittedValues)[[1]] ) )
                fittedValues <- fittedValues + groupMeans
            } else if( "exon" %in% names(dimnames(x)) ) {
                exonMeans <- apply2( x, "exon", mean )
                fittedValues <- t( t(fittedValues) + exonMeans )
            } else {
                print( x )
                stop( "Unexpected term encountered." )
            }
        }
        return( fittedValues )
    }else{
        stopifnot( "(Intercept)" %in% names(coefs) )
        stopifnot( "exonthis" %in% names(coefs) )
        allVars <- all.vars(frm)
        continuousVar <- allVars[!allVars %in% c("count", "exon")]
        interactionCoefName <- paste0( continuousVar, ":exonthis" )
        stopifnot( interactionCoefName %in% names(coefs) )
        mf[[continuousVar]]
        predictors <- unique( mf[[continuousVar]] )
        if( averageOutExpression ){
            fittedValues <- coefs["exonthis"] + coefs[interactionCoefName]*predictors
        }else{
            fittedValues <- coefs["(Intercept)"] + coefs["exonthis"] +
                (coefs[continuousVar] + coefs[interactionCoefName])*predictors
        }
        fittedValues <- matrix(fittedValues, ncol=1)
        colnames(fittedValues) <- "this"
#        rownames(fittedValues) <- sprintf("%s=%s", continuousVar, predictors)
        rownames(fittedValues) <- as.character(predictors)
        return( fittedValues )
    }
}


modelFrameSM <- function(object)
{
    mfSmall <- as.data.frame( colData(object) )
    mfSmall$exon <- factor( mfSmall$exon, levels=c("others", "this") )
    ##mfSmall$exon <- relevel( mfSmall$exon, "others" )
    mfSmall$dispersion <- NA
    mfSmall$count <- NA
    mfSmall
}

getEffectsForGeneBM <- function(geneID, groups, notNAs, countsAll,
                                disps, features, mf, frm, numsamples,
                                fitExpToVar, averageOutExpression=TRUE)
{
    rt <- groups %in% geneID & notNAs
    if( sum(rt) < 2 ){ return(NULL) }
    countsThis <- countsAll[rt,]
    rownames(countsThis) <- gsub("\\S+:", "", rownames(countsThis))
    dispsThis <- disps[rt]
    names(dispsThis) <- features[rt]
    numexons <- sum(rt)
    newMf <- mf[as.vector( sapply( split( seq_len(nrow(mf)), mf$sample ), "[", seq_len( numexons ) ) ),]
    newMf$exon <- factor( rep( features[rt], numsamples ) )
    for (i in seq_len(nrow(newMf))) {
       newMf[i, "dispersion"] <- dispsThis[as.character(newMf[i, "exon"])]
       newMf[i, "count"] <- countsThis[as.character(newMf[i, "exon"]), as.character(newMf[i, "sample"])]
    }
    newMf <- droplevels(newMf)
    coefficients <- fitAndArrangeCoefs( frm, balanceExons = TRUE, mf=newMf, fitExpToVar=fitExpToVar, geneID=geneID)
    if (is.null(coefficients)){
       return(coefficients)
    }
    ret <- t( getEffectsForPlotting(coefficients, averageOutExpression = averageOutExpression, 
        groupingVar = fitExpToVar, frm=frm, mf=newMf))
    rownames(ret) <- paste(geneID, rownames(ret), sep = ":")
    return(ret)
}

getEffectsForExonsSM <- function(index, frm, countsAll, disps,
                                 mfSmall, averageOutExpression=TRUE,
                                 fitExpToVar)
{
    mfSmall$count <- countsAll[index,]
    mfSmall$dispersion <- disps[index]
    coefs <- fitAndArrangeCoefs( frm, mf=mfSmall, balanceExons=FALSE, fitExpToVar=fitExpToVar, rownames(countsAll)[index])
    if( is.null(coefs) ){
        return(NULL)
    }
    getEffectsForPlotting(coefs,
            averageOutExpression=averageOutExpression,
            groupingVar=fitExpToVar, frm, mfSmall)[,"this"]
}

getEffectsForGene <- function( geneID, object, maxRowsMF, fitExpToVar)
{
    rt <- object$groupID %in% geneID
    sampleData <- object@sampleData
    if( is(object@sampleData[[fitExpToVar]], "numeric") ){
        maxRowsMF <- 0
    }
    numsamples <- nrow(object@sampleData)
    numexons <- sum(rt)
    featuresInGene <- object$featureID[rt]
    dispersions <- object$dispersion[rt]
    dispersions[is.na(dispersions)] <- 1e-08
    frm <- as.formula(paste("count ~", fitExpToVar, "* exon"))
    bigFlag <- numsamples*numexons < maxRowsMF
    if( bigFlag ){
        mf <- object@modelFrameBM
        mf <- mf[as.vector(sapply(split(seq_len(nrow(mf)), mf$sample), 
            "[", seq_len(numexons))), ]
        mf$exon <- factor(rep(featuresInGene, nrow(sampleData)))
        counts <- object$countData[rt,]
        rownames(counts) <- gsub("\\S+:", "", rownames(counts))
        names(dispersions) <- object$featureID[rt]
        for (i in seq_len(nrow(mf))) {
            mf[i, "dispersion"] <-
                dispersions[as.character(mf[i, "exon"])]
            mf[i, "count"] <-
                counts[as.character(mf[i, "exon"]), as.character(mf[i, "sample"])]
        }
        mf <- droplevels(mf)
        coefs <- fitAndArrangeCoefs(frm, balanceExons=TRUE, mf=mf, fitExpToVar=fitExpToVar, geneID)
        if( is.null(coefs ) ){
            return(NULL)
        }
        splicing <- t(getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=TRUE, frm=frm, mf=mf))
        expression <- t(getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=FALSE, frm=frm, mf=mf))
        rownames(splicing) <- sprintf("%s:%s", geneID, rownames(splicing))
        rownames(expression) <- rownames(splicing)
        list( expression=expression, splicing=splicing )
    }else{
        mf <- object@sampleData
        mf <- rbind( data.frame(mf, exon="this"), data.frame(mf, exon="others"))
        mf$exon <- factor( mf$exon, levels=c("others", "this") )
        ##mf$exon <- relevel( mf$exon, "others" )
        countsThis <- object$countData[rt,]
        countsOthers <- sapply( rownames( countsThis ),
                               function(x){
                                   colSums(countsThis[!rownames(countsThis) %in% x,,drop=FALSE])
                               })
        countsOthers <- t(countsOthers)
        stopifnot(all(rownames(countsThis) ==  rownames(countsOthers)))
        effects <- lapply( seq_len(numexons), function(x){
                   mf$count <- c( countsThis[x,], countsOthers[x,])
                   mf$dispersion <- dispersions[x]
                   coefs <- fitAndArrangeCoefs(frm, balanceExons=FALSE, mf=mf, fitExpToVar=fitExpToVar, geneID)
                   if( is.null(coefs) ){
                       return(NULL)
                   }
                   splicing <- getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=TRUE, frm=frm, mf=mf)[,"this"]
                   expression <- getEffectsForPlotting( coefs, groupingVar=fitExpToVar, averageOutExpression=FALSE, frm=frm, mf=mf)[,"this"]
                   list(splicing=splicing, expression=expression)
               })
        names(effects) <- rownames(object)[rt]
        splicing <- t(sapply(effects, "[[", "splicing"))
        expression <- t( sapply(effects, "[[", "expression" ))
        list( expression=expression, splicing=splicing )
    }
}