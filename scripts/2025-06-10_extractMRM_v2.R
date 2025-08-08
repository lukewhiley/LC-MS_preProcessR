library(tidyverse)

#read in delim data
msData <- read.table(
  file = file.choose(),
  header = FALSE, 
  sep = " ", 
  col.names = paste0("V",seq_len(20)), 
  fill = TRUE) %>% 
  as_tibble() %>%
  add_column(
    idx = c(1:nrow(.)),
    .before = 1
    )


#extract function list
functionList <- msData %>%
  add_column(Function = paste0(.$V1, " ", .$V2), .before = "V1") %>%
  filter(grepl("Function", Function)) %>%
  filter(!grepl(":", Function)) 

#extract required data

FunctionLoopOut <- NULL
for(idxFunction in functionList$Function){
  
  if(which(functionList$Function == idxFunction) < nrow(functionList)){
    loopIndex <- c(
      (functionList$idx[which(functionList$Function == idxFunction)]+1):
        (functionList$idx[which(functionList$Function == idxFunction)+1]-1)
    )
  }
    
  
  if(which(functionList$Function == idxFunction) == nrow(functionList)){
    loopIndex <- c(
      (functionList$idx[which(functionList$Function == idxFunction)]+1):
        nrow(msData)
      )
  }

  
  #create tibble
  #if(which(functionList$Function == idxFunction)==1){
 loopTibble <- msData[loopIndex,] 
 
 #getMRM
 mrmRowIndex <- which(rowSums(loopTibble == ">")>0)
 mrmColIndex <- c(which(colSums(loopTibble == ">")>0)-1,which(colSums(loopTibble == ">")>0)+1)
 mrmTibble <- loopTibble[mrmRowIndex,mrmColIndex]
 
 #getCompoundName
 compoundName <- loopTibble[nrow(loopTibble),] %>%
   select(-idx)
 compoundName <- compoundName[, which(grepl(unique(mrmTibble[,1]) %>% as.character() %>% gsub("\\..*", "", .), compoundName[1,]))] %>%
   gsub(",.*", "", .) %>%
   as.character() 
 
 if(length(compoundName) == 0){
   blankSpace <- which(loopTibble[nrow(loopTibble),] == "")
   compoundName <- paste0(loopTibble[nrow(loopTibble),-blankSpace] %>% select(-idx), collapse = ".")
   #compoundName <- compoundName %>% sub(" .*", "", .)
 }

 
 #molName
 molName <- loopTibble[nrow(loopTibble),] %>%
   select(-idx)
 molName <- molName[, (which(grepl(unique(mrmTibble[,1]) %>% as.character() %>% gsub("\\..*", "", .), molName[1,]))+1)] %>%
   as.character()
 molName <- molName %>% gsub("_.*", "", .) %>% gsub("-.*", "", .)
 
 if(grepl("(.?)SIL",compoundName)){molName <- "SIL"}
 if(grepl("(.?)QC",compoundName)){molName <- "QC"}
 
 if(length(molName)==0){
   blankSpace <- which(loopTibble[nrow(loopTibble),] == "")
   molName <- paste0(loopTibble[nrow(loopTibble),-blankSpace] %>% select(-idx), collapse = "-") 
     
   }
 
 
 #molName = "mcx-n"
 
 print(paste(idxFunction,molName, compoundName))
 
 #getCharge
 charge <- loopTibble[which(grepl("Ionization", loopTibble$V1)),]
 charge <- charge[which(grepl("ES", charge))] %>% as.character()
 if(charge == "ES+"){charge <- 1} else if(charge == "ES-"){charge <- -1} 
 
 #getTime
 rt <- loopTibble[which(grepl("Start", loopTibble$V1)),]
 startRT <- rt[which(grepl("to", rt))-1] %>% as.numeric()
 endRT <- rt[which(grepl("to", rt))+1] %>% as.numeric()
 explicitRT <- median(c(startRT, endRT))
 

 mrmTibbleOut <- NULL
 for(idxMRM in 1:nrow(mrmTibble)){
   mrmTibbleOut <- bind_rows(
     mrmTibbleOut,
     tibble(
       `Molecule List Name` = molName,
       `Molecule Name` = compoundName,
       `Precursor Mz` = mrmTibble[idxMRM,1] %>% as.numeric(),
       `Product Mz` = mrmTibble[idxMRM,2] %>% as.numeric(),
       `Explicit Retention Time` = explicitRT %>% round(2),
       `Precursor Charge` = charge,
       `Product Charge` = charge
     )
   )
}
 
 FunctionLoopOut <- bind_rows(
   FunctionLoopOut,
   mrmTibbleOut
 )
}

write_csv(FunctionLoopOut,
          file.choose(new=TRUE))



    
  