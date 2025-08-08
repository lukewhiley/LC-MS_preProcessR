library(tidyverse)

#read in delim ms method ata from a waters .exp file
msDataMethod <- read.table(
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
functionList <- bind_cols(
  msDataMethod %>% filter(grepl("FUNCTION", V1)) %>% unite(mrmFunction, V1, V2) %>% select(idx, mrmFunction),
  tibble(
    metabolite = msDataMethod %>% filter(grepl("(?i)CompoundName_1", V1)) %>%  
      mutate(V2 = gsub(",", "-", V2)) %>%
      mutate(V1 = gsub(".*,", "_ _,", V1)) %>%
      paste0("_1_", .$V1, "_2_", .$V2, "_3_", .$V3, "_4_") %>%
      sub(".*,", "", .) %>%
      sub("_2__3_", "", .),
    rtStart = msDataMethod  %>% filter(grepl("(?i)FunctionStartTime", V1))  %>% .$V1 %>% gsub(".*,", "", .),
    rtEnd = msDataMethod %>% filter(grepl("(?i)FunctionEndTime", V1)) %>% .$V1 %>% gsub(".*,", "", .),
    polarity = msDataMethod %>% filter(grepl("(?i)FunctionPolarity", V1)) %>%.$V1 %>% gsub(".*,", "", .),
    mrm1_q1 = msDataMethod %>% filter(grepl("(?i)SIRMass1", V1)) %>% .$V1  %>% gsub(".*,", "", .) %>% as.numeric(),
    mrm1_q3 = msDataMethod %>% filter(grepl("(?i)SIRMass_2_1", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric(),
    mrm1_ce = msDataMethod %>% filter(grepl("(?i)collisionenergy\\(V\\)_1", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric(),
    mrm2_q1 = NA,
    mrm2_q3 = NA,
    mrm2_ce = NA,
    mrm3_q1 = NA,
    mrm3_q3 = NA,
    mrm3_ce = NA
  )
  ) %>% 
  # create experiment channel
  add_column(
    .before = "metabolite",
    method = sub("_2_.*", "", .$metabolite),
    experiment = sub(".*_2_", "", .$metabolite) %>% sub("_3_.*", "", .)
  ) %>%
  mutate(metabolite = sub(".*_3_", "", metabolite)) %>%
  mutate(metabolite = sub("_4_.*", "", metabolite)) %>%
  mutate(metabolite = sub("\\.\\[.*", "", metabolite))
  #mutate(experiment = sub("_3_.*", "", experiment)) %>%
  #mutate(experiment = sub(".*-QC", "QC", experiment)) %>%
  #mutate(metabolite = sub("-QC", "", metabolite)) %>%
  #mutate(experiment = sub("_EXP", "", experiment)) %>%
  #mutate(experiment = sub("\\(.*", "", experiment))
  
#Combine experiment + metabolite into compound name for skyline
#Create RT window for skyline
functionList <- functionList %>% 
  unite("Precursor Name", experiment, metabolite,
        sep = " ", remove = FALSE) %>% 
  mutate(
    `Explicit Retention Time`     = (as.numeric(rtStart) + as.numeric(rtEnd)) / 2,
    `Explicit Retention Time Window` = as.numeric(rtEnd) - as.numeric(rtStart)
  ) %>%
  relocate(`Precursor Name`, .after = metabolite)

#tidy up method file by adding blank space as experiment where extraction of info has an error

for(idx_metabolite in functionList$metabolite){
  functionList[["experiment"]][which(functionList[["experiment"]]==idx_metabolite)] <- ""
}

#add experiment type for each SIL
for(idx_sil in functionList$metabolite[grepl("-SIL", functionList$metabolite)]){
  experiment <- functionList[["experiment"]][which(functionList[["metabolite"]]==gsub("-SIL.*","",idx_sil))]
  if(length(experiment) ==1){
    functionList[["experiment"]][which(functionList[["metabolite"]]==idx_sil)] <- experiment
  }
}


#add Q2 and Q3 (missmatch in length so can't do simple string extract)
idxStart = functionList$idx[1]
for(idxStart in functionList$idx){
  
  #find end index of function data
  idxEnd <- functionList$idx[which(functionList$idx == idxStart)+1]
  
  #extract data specific for function
  msDataLoop <- msDataMethod %>%
    filter(idx >= idxStart & idx <= idxEnd)
  
  #extract mrm2/mrm3 data
  mrm2_q1 <- msDataLoop %>% filter(grepl("(?i)SIRMass2", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  mrm2_q3 <- msDataLoop %>% filter(grepl("(?i)SIRMass_2_2", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  mrm2_ce <- msDataLoop %>% filter(grepl("(?i)collisionenergy\\(V\\)_2", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  mrm3_q1 <- msDataLoop %>% filter(grepl("(?i)SIRMass3", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  mrm3_q3 <- msDataLoop %>% filter(grepl("(?i)SIRMass_2_3", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  mrm3_ce <- msDataLoop %>% filter(grepl("(?i)collisionenergy\\(V\\)_3", V1)) %>% .$V1 %>% gsub(".*,", "", .) %>% as.numeric()
  
  #fill columns
  if(length(mrm2_q1)==1){functionList$mrm2_q1[which(functionList$idx == idxStart)] <- mrm2_q1}
  if(length(mrm2_q3)==1){functionList$mrm2_q3[which(functionList$idx == idxStart)] <- mrm2_q3}
  if(length(mrm2_ce)==1){functionList$mrm2_ce[which(functionList$idx == idxStart)] <- mrm2_ce}
  if(length(mrm3_q1)==1){functionList$mrm3_q1[which(functionList$idx == idxStart)] <- mrm3_q1}
  if(length(mrm3_q3)==1){functionList$mrm3_q3[which(functionList$idx == idxStart)] <- mrm3_q3}
  if(length(mrm3_ce)==1){functionList$mrm3_ce[which(functionList$idx == idxStart)] <- mrm3_ce}
}

functionList <- functionList %>%
  arrange(experiment, metabolite) %>%
  add_column(
    .before = 1,
    metaboliteIdx = 1:nrow(.)
  )

#turn long for rob
functionListLong <- bind_cols(
  #get q1
functionList %>%
  select(!contains("mrm"), contains("_q1")) %>%
  pivot_longer(cols = contains("_q1"), names_to = "mrm", values_to = "q1") %>%
  mutate(mrm = sub("_.*", "", mrm)) %>%
  mutate(mrm = sub("mrm", "", mrm)) %>%
  mutate(mrm = as.numeric(mrm)),
  #get q3
functionList %>%
  select(contains("_q3")) %>%
  pivot_longer(cols = contains("_q3"), names_to = "mrm", values_to = "q3") %>%
  select(q3),
#get ce
functionList %>%
  select(contains("_ce")) %>%
  pivot_longer(cols = contains("_ce"), names_to = "mrm", values_to = "ce") %>%
  select(ce)
) %>%
  filter(!is.na(q1))


#export csv
write_csv(functionListLong, file.choose(new = TRUE))
