library(tidyverse); library(reshape2)
library(pracma); library(grid); library(pwr)
library(randomForest); library(e1071); library(caret)

base_directory = '/Users/Mo/Desktop/lab/metabolomics/metabolomics_2017_submit' # REPLACE IF DATA SAVED ELSEWHERE
setwd(base_directory)

# make function to get stats on each metabolite
multi.fun = function(x) {
  c(mean = mean(x, na.rm =T), median = median(x, na.rm =T), var = var(x, na.rm =T), sd = sd(x, na.rm =T)) }

# preprocess dataset function to remove undetected mets
remove_NAs = function(raw) { # raw = raw file, unprocessed
  # if entire row is nas, remove 
  raw_naR = as.matrix(raw[rowSums(is.na(raw))!=(dim(raw)[2]), ]) #dim = 375  17
  return(raw_naR) }

# preprocess dataset function to scale to median, and imput missing values
scale_to_median = function(raw_naR) { # raw = raw file or output of remove_NAs
  stats = {}; storage.mode(raw_naR) = "numeric"
  for (i in 1:nrow(raw_naR)) {
    (stats = rbind(stats, multi.fun(as.numeric(raw_naR[i,])))) }
  row.names(stats) = row.names(raw_naR)
  data_for_scale = as.data.frame(merge(raw_naR,stats, by = 'row.names'))
  rownames(data_for_scale) = data_for_scale$Row.names
  data_for_scale = data_for_scale[,2:dim(data_for_scale)[2]]
  
  for (i in 1:nrow(data_for_scale)) { # scale with median met abundance
    for (j in 1:(ncol(data_for_scale)-4)) { # all columns except stats columns
      data_for_scale[i,j] = 
        ( data_for_scale[i,j] / data_for_scale$median[i] ) } }
  data_scaled = data_for_scale[,1:(ncol(data_for_scale)-4)] # remove stats columns
  return(data_scaled) }

# preprocess dataset function to imput missing values
imputing_missing = function(data_scaled, scaled) { # data = raw, raw_naR, or data_scaled file, scaled = 1/0 for if scaled or not
  remove_NAs = data_scaled  #replace NAs with minimum metabolite value
  for (i in 1:nrow(remove_NAs)) {
    min1 = min(remove_NAs[i,],na.rm= T)
    if (scaled == 1) { min2 = min1 } else if (scaled == 0) { min2 = min1/2
    } else {print('scaled must be a binary 1/0 for if data is already scaled')}
    for (j in 1:ncol(remove_NAs)) {
      if (is.na(remove_NAs[i,j])) { remove_NAs[i,j] = min2 }  }  }
  return(remove_NAs) }

# Function for generating multipanel plots
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  # see http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
  plots = c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout = matrix(seq(1, cols * ceiling(numPlots/cols)),ncol = cols, 
                    nrow = ceiling(numPlots/cols)) }
  if (numPlots==1) { print(plots[[1]]) } 
  else { # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx = as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col)) } } }

# Change rows in raw dataset 
prep_raw = function(raw) {
  rownames(raw) = gsub(" ",'',raw[,1]) # add row names (metabolites), removing spaces from strings
  raw = raw[,2:dim(raw)[2]] # remove column that was moved to row names
  return(raw) }

# Filter sample identifiers by what is in dataset
prep_sample_identifiers = function(raw) {
  setwd(base_directory) # open file containing sample measures
  identifiers = read.csv("sample_data_2017.csv", sep = ",", header = F, 
                         na.strings=c('',"NA"), stringsAsFactors = F)
  rownames(identifiers) = identifiers$V1 # add row names
  identifiers$V1 = NULL # remove column that was moved to row names
  identifiers = identifiers[,which(identifiers[3,]%in%colnames(raw))] # keep samples in raw
  return(identifiers) }

# Normalize metabolite abundances by sample measures (dna, protein abundance, or parasite number)
normalize = function(input,identifiers_file,type) {
  # input = input file (matrix) with rows as mets and columns as samples
  # identifiers_file = sample data with rows as sample measurements, columns as samples
  # type = type of normalization (dna, protein, para)
  if (strcmp(type, "dna")) {
    x = 6 
  } else if (strcmp(type, "protein")) {
    x = 5 
  } else if (strcmp(type, "para")) {
    x = 9 
  } else { print("ERROR, WRONG TYPE INPUT")}
  pre_norm1 = input
  for (i in 1:dim(pre_norm1)[1]) { # divide each metabolite abundance by it's sample value (dna, protein, or parasite #)
    for (j in 1:dim(pre_norm1)[2]) {
      (pre_norm1[i,j] = (as.numeric(pre_norm1[i,j])/as.numeric(identifiers_file[x,j]))) }}
  norm_file = pre_norm1
  return(norm_file) }


## remove characters from colnames if they will throw errors in classifiers
# use always for consistency
remove_characters_from_colnames = function(dataframe) {
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[.]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[(]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[)]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[']", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[)]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[:]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[;]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[[]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[]]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[*]", "_predicted", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[/]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("-", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("+", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[+]", "", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("[_][_]", "_", x) })
  colnames(dataframe) = lapply(colnames(dataframe), function(x) { gsub("__", "_", x) })
  
  for (i in 1:ncol(dataframe)) { # MAKE SURE STRINGS DONT START WITH NUMBERS 
    st = substr(colnames(dataframe)[i],0,1) # first character in string
    n_st = as.numeric(st)
    if (!is.na(n_st)) { colnames(dataframe)[i] = paste('N',colnames(dataframe)[i], sep = '')} } 
  
  return(dataframe)
}

custom_theme_pca = theme(
  plot.title = element_text(colour="black",size=9),
  axis.text.x = element_text(colour="black",size=6,angle=0,hjust=.5,vjust=.5,face="plain"), 
  axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0,face="plain"), 
  axis.title.x = element_text(colour="black",size=9,angle=0,hjust=.5,vjust=0,face="plain"), 
  axis.title.y = element_text(colour="black",size=9,angle=90,hjust=.5,vjust=.5,face="plain"),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = "NA", colour = "black"),
  axis.line = element_blank(),
  panel.spacing.x = unit(2, "lines"),
  legend.title=element_blank())

custom_theme_dist = theme(
  axis.text.y = element_text(colour="black",size=6), 
  axis.title.y = element_text(colour="black",size=9,angle=90,hjust=.5,vjust=0,face="plain"),
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = "NA", colour = "NA"),
  axis.line = element_blank(),
  panel.spacing.x = unit(0.05, "lines"),
  legend.title=element_blank())

custom_theme2 = theme(
  axis.text.x = element_text(colour="black",size=6), 
  axis.text.y = element_blank(), axis.title.y = element_blank(),
  axis.ticks.y = element_blank(), axis.ticks.x = element_blank(),
  axis.title.x = element_text(colour="black",size=9),
  panel.background = element_blank(),
  panel.grid.major = element_line(colour = "grey"),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = "NA", colour = "black"),
  axis.line = element_blank(),
  panel.spacing.x = unit(2, "lines"),
  plot.margin = unit(c(5,5,5,5),"mm"))

custom_theme3 = theme(
  plot.title = element_text(size = 9),
  axis.text.x = element_text(colour="black",size=6,angle=0,hjust=.5,vjust=.5,face="plain"), 
  axis.text.y = element_text(colour="black",size=6,angle=0,hjust=1,vjust=0,face="plain"), 
  axis.title.x = element_text(colour="black",size=9,angle=0,hjust=.5,vjust=0,face="plain"), 
  axis.title.y = element_text(colour="black",size=9,angle=90,hjust=.5,vjust=.5,face="plain"),
  panel.background = element_blank(),
  panel.grid.major = element_line(colour = "grey"),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(fill = "NA", colour = "black"),
  axis.line = element_blank(),
  panel.spacing.x = unit(2, "lines"))

'%notin%' <- function(x,y)!('%in%'(x,y))

median_center_scale_data <- function(x) {
  # x is matrix or dataframe with rows as samples and columns as mets
  sd_x = apply(x,2,sd,na.rm=TRUE)
  x_centered = apply(x, 2, function(y) (y - median(y,na.rm=TRUE)))
  x_centered_scaled = sweep(x_centered, 2, sd_x, "/")
  indx = setdiff((which(is.na(x_centered_scaled))),which(is.na(x_centered)))
  if (any(is.na(x_centered_scaled))) { 
    x_centered_scaled[indx] = 0}
  if (any(is.na(x_centered_scaled))) { warning('sd == 0')}
  return(x_centered_scaled)}

mean_center_scale_data <- function(x) {
  # x is matrix or dataframe with rows as samples and columns as mets
  sd_x = apply(x,2,sd,na.rm=TRUE)
  x_centered = apply(x, 2, function(y) (y - mean(y,na.rm=TRUE)))
  x_centered_scaled = sweep(x_centered, 2, sd_x, "/")
  indx = setdiff((which(is.na(x_centered_scaled))),which(is.na(x_centered)))
  if (any(is.na(x_centered_scaled))) { 
    x_centered_scaled[indx] = 0}
  if (any(is.na(x_centered_scaled))) { warning('sd == 0')}
  return(x_centered_scaled)}

# Get upper triangle of the correlation matrix
# modified from http://www.sthda.com/english/wiki/
# ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat) }

# also from www.sthda.com
reorder_cormat <- function(cormat){ # Use correlation between variables as distance
  dd <- as.dist(cormat); hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order] }


# make df values numeric instead of character or factor
make_numeric = function(df) { # doesn't really matter what columns or rows are
  indx = sapply(df, is.factor)
  df[indx] = lapply(df[indx], function(x) as.character(x)) 
  indx = sapply(df, is.character)
  df[indx] = lapply(df[indx], function(x) as.numeric(x)) 
  return(df) }

rM = base::rowMeans
