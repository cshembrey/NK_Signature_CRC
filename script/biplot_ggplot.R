
##############################################################################
## -- Author: Momeneh (Sepideh) Foroutan -- 
## -- First generated on: 24th of Sep 2019
## -- Last updated on: 24th of Sep 2019
##############################################################################

# biplot(): This function takes two expression data (e.g. TCGA and blood) 
#           and generates biplot comparing expression of genes in the two data sets

# this script also includes helper functions, mostly adopted from helper functions in singscore package.


##------------------ INPUTS:
## -- data.expr: a log scaled expression data matrix (e.g. TCGA CRC): This is a data set in
#     which we would like to have expression of candidate genes
#     Genes as row.names and samples in columns
#     The data should not have extra column of information, only gene expression

## -- data.noexpr: a log scaled expression data matrix (e.g. blood): This is a data set in
#     which we would NOT like to have expression of candidate genes
#     Genes as row.names and samples in columns
#     The data should not have extra column of information, only gene expression

## -- data.names: character vector of length two, describing the names of the
#     two data sets we are comparing; by default it is c("data.expr", "data.noexpr")
#     Note that the ata expr data name should be named first.

## -- data.expr.th: a value between 0 to 1, which is the percentile of the genes across
#     samples in the data.expr 

## -- data.noexpr.th: a value between 0 to 1, which is the percentile of the genes across
#     samples in the data.noexpr 

## -- annot: a NAMED vector, whose names are the same as row names in the two expression data.
#     This is used to colour data points and currently takes up to 8 groups 

## -- annot.name: a character, for the legend title

## -- point.labels: a NAMED vector, whose names are the same as row names in the two expression data.
#     This is used to label data points

## -- text.size: plot text size. Default is 1.2.


##------------------ OUTPUT:
## -- a biplot comparing expression of genes in the two data sets


##----------------- NOTE:
## The first three functions have been adopted from singscore codes.
## You can add lines using geom_vline(xintercept = 2) or geom_hline(yintercept = 2) later

################################## Helper functions #################################

##-------------- getTheme
getTheme <- function(rl = 1.2) {
  current_theme = ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.border = element_rect(colour = 'black', fill = NA),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = rel(rl) * 1.1),
      axis.text = element_text(size = rel(rl)),
      plot.title = element_text(size = rel(rl * 1.2)),
      strip.background = element_rect(fill = NA, colour = 'black'),
      strip.text = element_text(size = rel(rl)),
      legend.text = element_text(size = rel(rl)),
      legend.title = element_text(size = rel(rl), face = 'italic'),
      legend.position = 'bottom',
      legend.direction = 'horizontal'
    )
  
  return(current_theme)
}

MyCols <- c("#440154FF", "#482677FF", "#33638DFF", "#287D8EFF", "#20A387FF", "#55C667FF",
            "#95D840FF", "#FDE725FF")
##-------------- getColorScale
getColorScale <- function(annot) {
  #specify a discrete scale for categoricals
  if (is.factor(annot)) {
    #specify a discrete scale for categoricals
    if (length(levels(annot)) > 8) {
      warning('Too many levels of the annotation, using default ggplot2 colours')
      return(NULL)
    } else{
      return(scale_color_viridis_d(begin = 0, end = 0.65))
    }
  }
  
  #specify a continous scale for numerics
  if (is.numeric(annot)) {
    return(ggplot2::scale_colour_viridis_c())
  }
  
  return(NULL)
}

##-------------- processAnnotation
processAnnotation <- function(df, annot) {
  #return vector of empty strings if nothing is specified
  if (is.null(annot))
    annot = rep('', nrow(df))
  
  #characters will either be column names or character annotations
  if (is.character(annot) && length(annot) == 1 && annot %in% colnames(df)) {
    #a column name has been specified, extract the annotation
    annot = df[[annot]]
  }
  
  if (is.character(annot)) {
    #convert char annot to factor
    annot = as.factor(annot)
  }
  
  #check length of annotation matches number of observations
  stopifnot(length(annot) == nrow(df))
  
  #do nothing if numeric
  return(annot)
}





################################## biplot function #################################


biplot <- function(
  data.expr = logRPKM[1:200, 1:10],
  data.noexpr = logRPKM[1:200, 11:19],
  data.names = c("data.expr", "data.noexpr"),
  data.expr.th = 0.25,
  data.noexpr.th = 0.75,
  annot = NULL,
  annot.name = "Gene type",
  point.labels = NULL,
  text.size = 1.2) {
  
  
  data.expr   <- data.frame(data.expr)
  data.noexpr <- data.frame(data.noexpr)
  
  ## Remove any genes that have NA values (if any)
  data.expr   <- data.expr[complete.cases(data.expr),]
  data.noexpr <- data.noexpr[complete.cases(data.noexpr),]
  
  ## check for the genes in two data sets
  # if ( ! nrow(data.expr) == nrow(data.noexpr)){
  commonGenes <-
    as.character(intersect(row.names(data.expr), row.names(data.noexpr)))
  data.expr <- data.expr[commonGenes,]
  data.noexpr <- data.noexpr[commonGenes,]
  # }
  
  ## Add thresholds to two data based on the arguments
  data.expr$th <-
    apply(data.expr, 1, quantile, probs = data.expr.th)
  
  data.noexpr$th <-
    apply(data.noexpr, 1, quantile, probs = data.noexpr.th)
  
  if (!is.null(point.labels)) {
    point.labels <- point.labels[commonGenes]
  }
  
  #if no gene labels are provided
  if (is.null(point.labels)) {
    point.labels = ""
  } else{
    if (length(point.labels) != nrow(data.expr))
      stop(
        "point.labels must contain the same number of labels with the number of genes in the data"
      )
  }

  
  if (!is.null(annot)) {
    annot <- annot[commonGenes]
  }
  
  annot <- processAnnotation (data.expr, annot)
  

  ## generate a data with values and annots
  newdata <- cbind(data.expr$th, data.noexpr$th)
  newdata <- data.frame(newdata)
  colnames(newdata) <- c("data.expr", "data.noexpr")
  newdata$Labels <- point.labels
  newdata$Class <- annot
  
  ## set x and y axes limits
  minXY <- min(newdata[, 1:2])
  maxXY <- max(newdata[, 1:2])
  xylim <- c(minXY, maxXY)
  
  p0 <-
    ggplot2::ggplot(newdata, aes(x = data.noexpr, y = data.expr)) +
    geom_point() +
    xlim(xylim) +
    ylim(xylim) +
    xlab(data.names[2]) +
    ylab(data.names[1]) +
    getTheme(text.size) 

 
  p1 <- p0 +
    geom_point(
      aes(text = Labels, color = Class),
      shape = 21,
      fill = "white",
      size = 4,
      stroke = 4,
      data = newdata
    ) + 
    getColorScale(annot) +
    labs(colour = annot.name) +
    getTheme(text.size) 
   
  if (all(annot %in% '')) {
    p1 <-  p1 + guides(colour = FALSE)
  }
  
  p1 <- p1 +
    ggrepel::geom_label_repel(
      data = newdata, 
      mapping = aes(label = Labels, color = Class),
      size = 6, show.legend = FALSE, hjust = -0.3, vjust = 0.3)
  
  
  return(p1)
}



#---------------------------- Example of use:

# point.labels <- rownames(data.expr)
# names(point.labels) <- rownames(data.expr)
# 
# annot <- c(rep("A", nrow(data.expr)/2), rep("B", nrow(data.expr)/2))
# names(annot) <-  rownames(data.expr)

# pp <- biplot(
#   data.expr = logRPKM[1:50, 1:10],
#   data.noexpr = logRPKM[1:50, 11:19],
#   data.names = c("data.expr", "data.noexpr"),
#   data.expr.th = 0.25,
#   data.noexpr.th = 0.75,
#   annot = annot,
#   annot.name = "Gene type",
#   point.labels = point.labels,
#   text.size = 1.2)

