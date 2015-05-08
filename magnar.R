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

analyze_OrthoMCL <- function(mcl_file, var_file, resp_var, rndm1_eff, rndm2_eff, merge_var) {
  
  haplo_file <- "haplo_data.txt"
  
  py_submit <- paste(c("python convert_OrthoMCL.py", mcl_file, haplo_file), collapse=" ")
  system(py_submit)
  
  return(analyze_mtrx(haplo_file, var_file, resp_var, rndm1_eff, rndm2_eff, merge_var))
}

#'rndm1_eff is handled different from rndm2_eff
#'
#'
analyze_mtrx <- function(haplo_file, var_file, resp_var, rndm1_eff, rndm2_eff, merge_var) {
  
  cat(" _____ _____ _____ _____ _____ _____ \n")
  cat("|     |  _  |   __|   | |  _  | __  |\n")
  cat("| | | |     |  |  | | | |     |    -|\n")
  cat("|_|_|_|__|__|_____|_|___|__|__|__|__|\n\n")
 
  library(nlme)
  
  cat("Importing Data\n")
  haplo_mtrx <- parse_haplo(haplo_file)
  
  haplo_names <- parse_names(haplo_file)
  
  tx <- read.table(var_file, header=T)
  
  cat("Merging Files\n")
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #phylogenetic distribution group (move this?)
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  
  cat("Running Analysis\n")
  output <- matrix(, nrow=count, ncol=3)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  cat("Percent Complete:\n")
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[,which(colnames(haplo_tx)==name)]
    
    #fitting the linear mixed model
    lmm <- try(lme(resp_var~fix, random=~1|rndm1/rndm2, data=sub),T)
    p_val <- try(summary(lmm)$tTable[2,5], T)
    
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg)),T)
      out_cnt <- out_cnt + 1
    }
    
    if(((i-1) %% 100) == 0) 
      cat(paste(round(((i-1)/num_pdg*100), digits=2), "%__", sep=""))
  }
  
  cat("Finished!\n")
  return(output)
}
  
  
  
  
  
  
  
  
  
  




