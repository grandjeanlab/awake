
library(tidyverse)

meta_data<-'/home/traaffneu/joagra/code/awake/assets/tables/rat_metadata.tsv'
folder_data<-'/project/4180000.36/awake/bids/rat'


df<-read_tsv(meta_data) %>% select(rodent.ds, rodent.sub, rodent.session, rodent.run, MRI.TR)

df$scan.exist<-NA 
#list all .nii.gz in  folder_data
files<-list.files(folder_data, pattern = "bold.nii.gz", recursive = TRUE, full.names = TRUE)

no_match<-c()

#for every file in files, get the name of the file and split it by "_"
#assign the values to variables 'sub', 'ses', 'run'

for(i in files){
 #get the basename of i
  base<-basename(i)
  #split the basename by "_"
  base<-strsplit(base, "_")[[1]]
  #find which element of base is the subject, session and run
  sub<-base[grep("sub", base)] %>% 
    gsub("sub-0", "", .) %>% 
    as.numeric()
  ses<-base[grep("ses", base)] %>% 
    gsub("ses-", "", .) %>%
    as.numeric()
  run<-base[grep("run", base)] %>% 
    gsub("run-", "", .) %>%
    as.numeric()
  
  #find which row in df matches the sub, ses, and run. 
  #if there is a match, assign 1 to scan.exist

  df_id<-which(df$rodent.sub==sub & df$rodent.session==ses & df$rodent.run==run)
  if(length(df_id)==0){
    print(paste0("No match for ", i))
    no_match<-c(no_match, i)
    next
  }
  df$scan.exist[df_id]<-1
  tr<-df$MRI.TR[df_id]
  #replace the nii.gz with json
  json<-gsub(".nii.gz", ".json", i)
  json_content<-paste0('{"RepetitionTime":', tr, '}')
  write(json_content, json)
}

#list all df rows that have scan.exist==NA
missing<-df %>% filter(is.na(scan.exist))
print(missing)

#save missing to ../tables/mouse_missing.tsv
write_tsv(missing, 'assets/tables/mouse_missing.tsv')
#save no_match to ../tables/mouse_no_match.txt
write(no_match, 'assets/tables/mouse_no_match.txt')

