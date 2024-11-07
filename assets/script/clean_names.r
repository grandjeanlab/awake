data_folder<-'/project/4180000.36/awake/bids/rat/'

#list all the files ending with /anat/*.nii.gz in the data_folder

anat_files<-list.files(data_folder, recursive = TRUE, pattern = '*T2w.nii.gz', full.names = TRUE)

#loop over the files and rename them
for(i in anat_files){
  #if the file name contains the string 'run' remove the string 'run' the following characters until the next '_'
  new_name<-gsub('run[-0-9]+_','',i)
  #move i to new_name
  file.rename(i,new_name)
}


#list all the files ending with *bold.nii.gz in the data_folder
bold_files<-list.files(data_folder, recursive = TRUE, pattern = '*bold.nii.gz', full.names = TRUE)

#loop over the files and rename them
for(i in bold_files){
  #if the file name contains the string 'run' remove the string 'run' the following characters until the next '_'
  new_name<-gsub('task-[A-Za-z]+_','',i)
  #move i to new_name
  file.rename(i,new_name)
}

bold_files<-list.files(data_folder, recursive = TRUE, pattern = '*bold.nii.gz', full.names = TRUE)

#loop over the files and rename them
for(i in bold_files){
  #if the file name contains the string 'run' remove the string 'run' the following characters until the next '_'
  new_name<-gsub('run-0','run-',i)
  #move i to new_name
  file.rename(i,new_name)
}


bold_files<-list.files(data_folder, recursive = TRUE, pattern = '*bold.nii.gz', full.names = TRUE)

#loop over the files and rename them
for(i in bold_files){
  #if the file name contains the string 'run' remove the string 'run' the following characters until the next '_'
  new_name<-gsub('_bold','_task-rest_bold',i)
  #move i to new_name
  file.rename(i,new_name)
}

