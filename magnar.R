#' Set Working Directory
#' 
#' @param dir Working directory where io files are stored
#' @export
mgwas_wd <- function(dir) {
  setwd(dir)
}

#' Read Genetic Data
#' 
#' @param file File path containing the genetic data
#' @return data Vector containing the genotype mapped to specific species that have that genotype
#' @export
parse_haplo <-function(file) {
  
  table <- as.matrix(read.table(file, header=F, sep="|")) #now it will separate into two listsls
  
  haplo <- table[,1]
  haplo_list <- strsplit(haplo, split = " ")
  haplo_matrix <- matrix(unlist(haplo_list), ncol=length(haplo_list[[1]]), byrow=T)
  haplo_matrix <- t(haplo_matrix)
  
  return(haplo_matrix)
}

parse_names <-function(file) {
  
  table <- as.matrix(read.table(file, header=T, sep="|"))
  
  names <- table[,2]
  names_list <- strsplit(names, split = ",")
  
  return(names_list)
}

merge_tx <-function(file, haplo_matrix) {
  
  tx <- read.table(file, header=T)
  
  haplo_tx <- merge(tx, haplo_matrix, by.y="V1", by.x="Treatment", all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx$tgl_fly>0))   
  
  return(haplo_tx)
}

create_haplo_matrix <-function(haplo_file, var_file) {
  
  haplo_matrix <- parse_haplo(haplo_file)
  haplo_tx <- merge_tx(var_file, haplo_matrix)
  
  return(haplo_tx)
}

#'Model keys:
#'lme:f1+f2|r1+r2 (DONE) 
#'lme:f1|r1/r2    (DONE)
#'lme:f1|r1+r2    (DONE)
#'lme:f1|r1       (DONE)
#'lm:f1           (DONE)
#'wx:f1           (DONE)
#'
#'@author Hayden Smith
#'@description Something about the function...
#'@usage Something about the usage
#'@param Something about the param
#'@examples Some sore of example
#'
analyze_OrthoMCL <- function(mcl_file, var_file, merge_var, species_name, model, resp=NULL, fix2=NULL, rndm1=NULL, rndm2=NULL) {
  
  cat(" _____ _____ _____ _____ _____ _____ _ _ _ _____ _____ \n")
  cat("|     |  _  |   __|   | |  _  |     | | | |  _  | __  |\n")
  cat("| | | |     |  |  | | | |     | | | | | | |     |    -|\n")
  cat("|_|_|_|__|__|_____|_|___|__|__|_|_|_|_____|__|__|__|__|\n")
  
  haplo_file <- "haplo_data.txt"
  
  py_submit <- paste(c("python convert_OrthoMCL.py", mcl_file, haplo_file), collapse=" ")
  system(py_submit)
  
  if(model == "lme:f1+f2|r1+r2")   #DONE
    mtrx <- analyze.ffrr(haplo_file, var_file, merge_var, species_name, resp, fix2, rndm1, rndm2)
  else if(model == "lme:f1|r1/r2") #DONE
    mtrx <- analyze.frr.div(haplo_file, var_file, merge_var, species_name, resp, rndm1, rndm2)
  else if(model == "lme:f1|r1+r2") #DONE
    mtrx <- analyze.frr.plus(haplo_file, var_file, merge_var, species_name, resp, rndm1, rndm2)
  else if(model == "lme:f1|r1")    #DONE
    mtrx <- analyze.fr(haplo_file, var_file, merge_var, species_name, resp, rndm1)
  else if(model == "lm:f1")        #DONE
    mtrx <- analyze.f(haplo_file, var_file, merge_var, species_name, resp)
  else if(model == "wx:f1")        #DONE
    mtrx <- analyze.wilcox(haplo_file, var_file, merge_var, species_name, resp)
  else 
    print("Could not find a correct match for your model declaration")
  
  rm_submit <- paste("rm", haplo_file, sep=" ")
  system(rm_submit)
  return(mtrx)
}

analyze.ffrr <- function(haplo_file, var_file, merge_var, species_name, resp_var, fix_var2, rndm1_eff, rndm2_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=9)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$fix2 <- haplo_tx[,which(colnames(haplo_tx)==fix_var2)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + fix2 + (1|rndm1) + (1|rndm2), data=sub),T)
    l1.1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l1.2 <- try(summary(l1.1),T)
    p_val1 <- try(l1.2$test$pvalues[1],T)
    
    l2.1 <- try(glht(lmm, mcp(fix2="Tukey")),T)
    l2.2 <- try(summary(l2.1),T)
    p_val2 <- "variable not a factor"
    try(p_val2 <- l2.2$test$pvalues[1],T)
    
    p_val2_corrected <- "ibid"
    try(p_val2_corrected <- p_val2*num_pdg,T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val1, (p_val1*num_pdg), p_val2, p_val2_corrected, mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=9, ] 
  colnames(output_clean) <- c("pdg", "p-val1", "corrected_p-val1", "p-val2", "corrected_p-val2", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.frr.plus <- function(haplo_file, var_file, merge_var, species_name, resp_var, rndm1_eff, rndm2_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm1) + (1|rndm2), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.frr.div <- function(haplo_file, var_file, merge_var, species_name, resp_var, rndm1_eff, rndm2_eff) {
 
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm1/rndm2), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}
  
analyze.fr <- function(haplo_file, var_file, merge_var, species_name, resp_var, rndm_eff) {
  
  library(lme4)
  library(multcomp)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm <- haplo_tx[,which(colnames(haplo_tx)==rndm_eff)]
  
  sub$rndm <- factor(sub$rndm)
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear mixed model
    lmm <- try(lmer(resp_var~fix + (1|rndm), data=sub),T)
    l1 <- try(glht(lmm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.f <- function(haplo_file, var_file, merge_var, species_name, resp_var) {
  
  library(multcomp)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear model
    lm <- try(lm(resp_var~fix, data=sub),T)
    l1 <- try(glht(lm, mcp(fix="Tukey")),T)
    l2 <- try(summary(l1),T)
    p_val <- try(l2$test$pvalues[1],T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

analyze.wilcox <- function(haplo_file, var_file, merge_var, species_name, resp_var) {
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #number of phylogenetic distribution groups
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=7)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  
  cat("\tPercent Complete:\n\t")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    #getting column names for expiramental fixed effect
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[, which(colnames(haplo_tx)==name)]
    sub$species <- haplo_tx[, which(colnames(haplo_tx)==species_name)]
    #View(sub)
    
    #fitting the linear model
    wix <- try(wilcox.test(resp_var~fix, data=sub),T)
    p_val <- try(wix$p.value,T)
    
    #calculating meta-data
    mean_calc <- aggregate(resp_var~fix, sub, mean) #matrix with mean_contain and mean_missing
    mean_calc$fix <- as.numeric(as.character(mean_calc$fix))
    mean_calc2 <- mean_calc[order(mean_calc$fix, decreasing=F),]
    mean_contain <- mean_calc2[2,2]
    mean_missing <- mean_calc2[1,2]
    
    taxa <- aggregate(as.numeric(as.character(fix))~species, sub, mean)
    colnames(taxa) <- c("taxa_name", "presence") 
    
    taxa_mtrx1 <- droplevels(subset(taxa, taxa$presence==1))
    taxa_contain <- paste(taxa_mtrx1$taxa_name, collapse="|")
    
    taxa_mtrx0 <- droplevels(subset(taxa, taxa$presence==0))
    taxa_missing <- paste(taxa_mtrx0$taxa_name, collapse="|")
    
    #adding data to output matrix
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg), mean_contain, mean_missing, taxa_contain, taxa_missing),T)
      out_cnt <- out_cnt + 1
    }
    
    #printing progress
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  cat("Finished!\n")
  
  cat("Cleaning Final Output")
  output_clean <- output[rowSums(is.na(output))!=7, ] 
  colnames(output_clean) <- c("pdg", "p-val", "corrected_p-val", "mean_pdgContain", "mean_pdgMiss", "taxa_contain", "taxa_miss")
  
  return(output_clean)
}

print_mcl <- function(mtrx, filename) {
  write.tsv(mtrx, filename)
}
  
join_repset <- function(reps_file, mcl_mtrx) {
  
  library("seqinr")
  
  reps_fa <- read.fasta(file=reps_file, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
  fa_mtrx <- matrix(, nrow=length(reps_fa), ncol=4)
  colnames(fa_mtrx) <- c("pdg", "id", "annot", "seq")
  
  for(i in 1:length(reps_fa)) {
    info <- strsplit(getAnnot(reps_fa[[i]]), split = "\t")
    fa_mtrx[i,] <- c(info[[1]][1], info[[1]][2], info[[1]][3], reps_fa[[i]][1])
  }
  
  mcl_reps <- merge(mcl_mtrx, fa_mtrx, by="pdg", all=F)
  
  return(mcl_reps)
}

download_packages <- function() {
  
  try(install.packages("lme4", dependencies=T), F)
  
  try(install.packages("multcomp", dependencies=T), F)
  
  try(install.packages("seqinr", dependencies=T), F)
  
}

format_MCLfastas <- function(fa_dir, filename, genbnk_id) {
  
  library("seqinr")
  
  outfile <- paste(c(fa_dir, filename), collapse="")
  
  files <- dir(fa_dir, pattern = ".fasta")
  if(filename %in% files) {
    file.remove(outfile)
  }
  
  for(i in 1:length(files)) {
    
    if(files[i] != filename) {
      
      abs_path <- paste(c(fa_dir, files[i]), collapse="")
      cat(abs_path)
      cat("\n")
      
      id <- strsplit(files[i], split="\\.")
      
      reps_fa <- read.fasta(file=abs_path, as.string=T, forceDNAtolower=F, seqonly=F, strip.desc=T)
      
      for(j in 1:length(reps_fa)) {
        
        info <- strsplit(getAnnot(reps_fa[[j]]), split="\\|")
        
        mcl_info <- paste(c(id[[1]][1], info[[1]][genbnk_id]), collapse="|")
        
        write.fasta(reps_fa[[j]][1], mcl_info, outfile, open="a")
      }
    }
  }
  
  return(outfile)
}

pick_repseq <- function(mcl_file, all_fa) {
  
  data <- as.matrix(read.delim(mcl_file, header=F, sep=" ", row.names=1))
  data <- t(data)
  colnames(data) <- sub(":", "", colnames(data))
  
  View(data)
  
  #randomly picks rep sequences from mcl_file and formatted MCLfastas from formated_MCLfastas (that I uploaded to MCL) OR de novo MCL 
}


  
  
  
  
  
  




