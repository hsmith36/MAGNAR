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
  
  # if TAG < 0, throw out the value and should check for negative values
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx$tgl_fly>0))   
  
  return(haplo_tx)
}

create_haplo_matrix <-function(haplo_file, var_file) {
  
  haplo_matrix <- parse_haplo(haplo_file)
  haplo_tx <- merge_tx(var_file, haplo_matrix)
  
  return(haplo_tx)
}

#'rndm1_eff is handled different from rndm2_eff
#'
#'
#haplo_file <- haplo_50
#var_file <- var_file
#resp_var <- "tgl_fly"
#rndm1_eff <- "Treatment"
#rndm2_eff <- "Experiment"
#merve_var <- "Merger"
analyze_mtrx <- function(haplo_file, var_file, resp_var, rndm1_eff, rndm2_eff, merge_var) {
  
  print("Start")
  library(nlme)
  
  haplo_mtrx <- parse_haplo(haplo_file)
  #View(haplo_mtrx)
  
  haplo_names <- parse_names(haplo_file)
  #print(haplo_names)
  
  tx <- read.table(var_file, header=T)
  
  #View(tx)
  #print(colnames(tx))
  
  haplo_tx <- merge(tx, haplo_mtrx, by.y="V1", by.x=merge_var, all=F)
  haplo_tx <- droplevels(subset(haplo_tx, haplo_tx[resp_var] > 0))
  
  #View(haplo_tx)
  
  num_pdg <- dim(haplo_mtrx)[2] - 1 #phylogenetic distribution group (move this?)
  #print(num_pdg)
  
  count <- 0
  for(i in haplo_names) {
    for(j in haplo_names[i])
      count <- count + 1
  }
  #print(count)
  output <- matrix(, nrow=count, ncol=3)
  
  sub <- matrix(, nrow=0, ncol=0)
  sub$resp_var <- haplo_tx[,which(colnames(haplo_tx)==resp_var)]
  sub$rndm1 <- haplo_tx[,which(colnames(haplo_tx)==rndm1_eff)]
  sub$rndm2 <- haplo_tx[,which(colnames(haplo_tx)==rndm2_eff)]
  
  sub$rndm1 <- factor(sub$rndm1)
  sub$rndm2 <- factor(sub$rndm2)
  
  #Print percentage complete
  out_cnt <- 1
  for(i in 2:(num_pdg + 1)) { 
    
    name <- paste("V", i, sep="")
    sub$fix <- haplo_tx[,which(colnames(haplo_tx)==name)]
    print(name)
    
    #fitting the linear mixed model
    lmm <- try(lme(resp_var~fix, random=~1|rndm1/rndm2, data=sub),T)
    p_val <- try(summary(lmm)$tTable[2,5], T)
    
    for( j in haplo_names[[i-1]]) {
      
      try(output[out_cnt,] <- c(j, p_val, (p_val*num_pdg)),T)
      out_cnt <- out_cnt + 1
    }
  }
  
  return(output)
}
  
  
  
  
  
  
  
  
  
  




