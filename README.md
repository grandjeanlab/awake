# AwakeRodent, the multi-center, multi-species, multi-modality study!

## 1.1 Executive summary
Awake whole-brain functional imaging in rodents is the new frontier for translational and biomedical research. In previous task-free functional MRI dataset comparisons in mice [1] and rats [2], awake datasets showed fewer patterns of plausible functional connectivity compared to anesthetized counterparts. Similar to our previous comparisons, we posit that there are laboratory and protocol differences yielding varying outcomes in rodent awake functional imaging. We want to exploit this feature to identify practices that enhance the detection of plausible functional networks in awake rodents and establish best practices to potentiate our community. 

[1] https://pubmed.ncbi.nlm.nih.gov/31614221/   
[2] https://www.biorxiv.org/content/10.1101/2022.04.27.489658v1   

## 1.2 Goal
To aggregate awake rodent functional datasets in one dataset collection
To describe reference population parameters distributions (e.g. motion, tSNR) 
To provide evidence-based recommendations for awake functional imaging 
To compare functional Magnetic Resonance Imaging and functional Ultrasound
To foster collaborations and discussion within the community

## 1.3 Collaborative model
Every laboratory is invited to participate with one or multiple datasets (see 1.4 Data specification). Invitation to collaborate is made through email, announcements at conferences (ESMI, ISMRM, ESMRMB), and on Twitter (@grandjeanlab). See the list of currently listed laboratories (section 2). 

The collaborative model is as follows: 
- Every contributing laboratory can nominate any number of collaborators who contributed to the dataset (e.g., acquisition, funding, management, etc..), split between junior and senior collaborators. 
- Joanes Grandjean (JG, joanes.grandjean@radboudumc.nl) will put the dataset collection together and perform the primary analysis. 
- Every collaborator is invited to further contribute to the preregistration proposal, the analysis, and manuscript preparation.
- Additional collaborators are invited by JG to cover additional expertise (open science, analysis, editorial). 

In the ensuing journal publication, the author list will be as such: 
junior with an extra contribution, [junior collaborator in alphabetical order], [senior collaborator in alphabetical order], senior with an extra contribution. JG reserves the right to nominate the first and last authors. 

## 1.4 Data specification
Inclusion criteria: Individual datasets each consist of n >= 4 task-free functional MRI or functional ultrasound scan obtained in either mice or rats, without genetic or experimental manipulation, any strain, any gender, any age, any weight. The datasets have not undergone preprocessing, beyond data conversion. Datasets are provided with meta-data (see section 3.3.2 Measured variables). MRI data are additionally provided with anatomical scans.   
   
Exclusion criteria: The dataset must not require additional processing beyond that performed within RABIES. Scans must be acquired with a sampling rate >= 0.25 Hz (TR = 4s) and at least cover the following field-of-view: 
FOV mice:  AP 1 to -2 mm, LR 3.5 to -3.5 mm 
FOV rats: AP 1.5 to -3 mm, LR 5 to -5 mm 

The AP boundaries correspond to the most anterior part of the corpus callosum to the beginning of the hippocampus and are justified as this is where the regions of interest will be placed. 

The material transfer will be performed under a tacit e-mail agreement, or if the data provider host institution requires it, a Material Transfer Agreement will be signed between the provider host institution and the Radboud University Medical Center, Nijmegen, The Netherlands (JG’s host institution). 

## 1.5 Deliverables
OSF.io preregistered study
Publicly available collective dataset on openneuro.org consisting of individual lab dataset (if not already available elsewhere) 
Publicly available code to replicate the study on github.com
Biorxiv preprint
Journal publication

## 1.6 Timeline
November 2022: Information sessions
December 2022: Pre-registration. Study onset. 
December 2023: Manuscript preparation. 

## 1.7 Communication
Important milestones: Group email (bcc). 
Preregistration preparation: This google doc template. Every collaborator is invited to comment/edit.
Analysis update: Public Github repository github.com/grandjealab/awake. Every collaborator is invited to comment/raise an issue. 
Manuscript preparation: Google doc. 

## 1.8 Data storage
Short-term data storage on google drive/dropbox/surfdrive to transfer data. 
Mid-term data storage at the Donders Institute high-performance computer in BIDS format
Long-term data storage at the openneuro.org repository in BIDS format, publicly available when the first preprint is submitted. 

## 1.9 Data processing
Data will be preprocessed using RABIES (https://github.com/CoBrALab/RABIES), a BIDS-based software based on the fMRIprep pipeline. Data will be co-registered into either the DSURQE space (mice) or SIGMA space (rats). Ultrasound data will be registered to a vascular template in the corresponding spaces [1,2]. The first-level analysis will be performed within RABIES. Analysis, plotting, and statistics will be done with Python. 
For a detailed procedure, see section 4 Detailed analysis.

[1] https://www.nature.com/articles/s41596-021-00548-8   
[2] https://github.com/nerf-common/whole-brain-fUS   

## 1.10 Retraction
At any time and without justification, collaborators can decide to be removed from the study. Their dataset will be deleted and their results not used in the final results. This needs to be communicated to JG. 
