# global analysis
Joanes Grandjean
Invalid Date

``` python
%autoindent 
```

``` python
import polars as pl
#import plotting modules
import matplotlib.pyplot as plt
import seaborn as sns


pl.Config(
    tbl_formatting="MARKDOWN",
    tbl_hide_column_data_types=True,
    tbl_hide_dataframe_shape=True,
    tbl_rows=100,

)

analysis_list = [ "gsr1", "gsr2", "gsr3","wmcsf1", "wmcsf2", "wmcsf3", "aCompCor1", "aCompCor2", "aCompCor3" ]  

cat_index = ["Specific", "Non-specific", "Spurious", "No"]

seed_ref = "s1_r"
seed_specific = "s1_l"
seed_unspecific = "aca_r"
seed_thalamus = "vpm_r"

roi_list=["s1","thal"]
```

## Mouse analysis.

``` python
rodent_list = ["mouse", "rat"]
for rodent in rodent_list:
#rodent='rat'
  print("#### NOW DOING " + rodent + " ####")
  df = pl.read_csv("../assets/tables/"+rodent+"_metadata_process.tsv", separator="\t")
  df_summary=pl.read_csv("../assets/tables/"+rodent+"_summary_processed.tsv", separator="\t")

  print("summary of the data that we collected")
  print("we processed " + str(df_summary["rodent.ds"].count()) + " datasets") 
  print("totalling "+ str(df_summary["total_run"].sum()) +" runs")
  print("from "+str(df_summary["total_animal"].sum()) +" mice")
  print("the smallest dataset had "+ str(df_summary["total_run"].min())+" runs") 
  print("the largest dataset had "+str(df_summary["total_run"].max())+" runs")
  print("we could processed "+str(df_summary["total_included"].sum()) + "/" + str(df_summary["total_run"].sum())+ " runs.")

  print("below is a summary of the data included per dataset")
#to add the summary of the data included
  print(df_summary.select("rodent.ds", "total_run", "total_animal", "total_included").sort(by="rodent.ds"))

  print("information about sex ratio")
  print("the datasets contained "+str(df_summary["male"].sum())+ " males and " +str(df_summary["female"].sum())+ " females")
  print("that corresponds to "+str(round(100*df_summary["female"].sum()/(df_summary["female"].sum()+df_summary["male"].sum()),2))+"% females ")
  print("information about animal handling")
  print(str((df_summary["headplate"] == 'y').sum()) + " datasets used headplates")
  print(str((df_summary["restrained"] == 'y').sum()) + " datasets used body restraining")
  print(str((df_summary["anesthesia"].is_in(['Isoflurane', 'Sevoflurane'])).sum()) + " datasets used anesthesia before acquisition")
  print(str((df_summary["exp.gender"] == 'm').sum()) + " datasets were collected by males, " + str((df_summary["exp.gender"] == 'f').sum())+ " by women")

  print(df_summary.select("rodent.ds", "headplate", "restrained", "anesthesia","exp.gender").sort(by="rodent.ds"))

  print("information about the scanner and sequence")
  print("lowest field strength was " + str(df_summary["field_strength"].min()) + "T")
  print("highest field strength was " + str(df_summary["field_strength"].max()) + "T")
  print(df_summary.select("rodent.ds", "field_strength", "sequence", "TE").sort(by="rodent.ds"))

#first, let's extract some infomation about motion and summarize it per dataset
  print("#### MOTION ANALYISIS ####")
  print("mean fd across all "+rodent+" datasets")
  print(df.select("fd.mean").mean())
  print("mean fd per dataset")
  print(df_summary.select("rodent.ds", "fd.mean").sort(by="rodent.ds"))
#now we run the analysis per denoising style, we extract the number of dropped frames, the s1-s1, s1-aca, and s1-thal correlations. finally we estimate connectivity specificity
  print("Number of dropped frames for each dataset and denoising method")
  print("This corresponds to the following rabies flags for mice")
  print("#gsr1: --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true")
  print("#gsr2: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true")
  print("#gsr3: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false")
  print("#### DENOISE ANALYSIS ####")
  print("dropped frames per dataset")
  print(df_summary.select("rodent.ds", "total.frames", "dropped.frames.gsr1", "dropped.frames.gsr2", "dropped.frames.gsr3").sort(by="rodent.ds"))
  print("#### FC specifiticy analysis ####")
  for analysis in analysis_list:
    print(" ")
    print("overall FC specificity for "+ analysis)
    print(df["s1.cat."+analysis].value_counts().with_columns(pl.col("count")/pl.sum("count")).sort(by='count', descending=True))

  print("for gsr processing, heavy scrubbing lead to  "+str((df["s1.cat.gsr1"]=="Specific").sum()) + "/" +str((df["s1.cat.gsr1"]=="Specific").count())  + " Specific scans")
  print("for gsr processing, light scrubbing lead to  "+str((df["s1.cat.gsr3"]=="Specific").sum()) + "/" +str((df["s1.cat.gsr3"]=="Specific").count())  + " Specific scans")

  print("for aCompCor processing, heavy scrubbing lead to  "+str((df["s1.cat.aCompCor1"]=="Specific").sum()) + "/" +str((df["s1.cat.aCompCor1"]=="Specific").count())  + " Specific scans")
  print("for aCompCor processing, light scrubbing lead to  "+str((df["s1.cat.aCompCor3"]=="Specific").sum()) + "/" +str((df["s1.cat.aCompCor3"]=="Specific").count())  + " Specific scans")

  print("s1-s1 specificity analysis per dataset (only showing specific values)")
  print(df_summary[["rodent.ds"]+["s1."+analysis+".Specific" for analysis in analysis_list]].sort(by="rodent.ds"))

  print("it seems that less scrubbing lead to better outcomes")
  print("for gsr processing, "+str((df_summary["s1.gsr1.Specific"] < df_summary["s1.gsr3.Specific"]).sum())+"/"+str(df_summary["s1.gsr1.Specific"].count())+" dataset performed better with less scrubbing")
  print("for aCompCor processing, "+str((df_summary["s1.aCompCor1.Specific"] < df_summary["s1.aCompCor3.Specific"]).sum())+"/"+str(df_summary["s1.aCompCor1.Specific"].count())+" dataset performed better with less scrubbing")
  print("for comparing gsr to aCompCor processing, "+str((df_summary["s1.gsr3.Specific"] < df_summary["s1.aCompCor3.Specific"]).sum())+"/"+str(df_summary["s1.aCompCor3.Specific"].count())+" dataset performed better with aCompCor compared to gsr")
  print("for comparing gsr to aCompCor processing, "+str((df_summary["s1.gsr3.Specific"] > df_summary["s1.aCompCor3.Specific"]).sum())+"/"+str(df_summary["s1.aCompCor3.Specific"].count())+" dataset performed worst with aCompCor compared to gsr")
```

    #### NOW DOING mouse ####
    summary of the data that we collected
    we processed 13 datasets
    totalling 887 runs
    from 116 mice
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 744/887 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included |
    |-----------|-----------|--------------|----------------|
    | 1001      | 4         | 4            | 4              |
    | 1002      | 8         | 8            | 8              |
    | 1003      | 12        | 12           | 12             |
    | 1004      | 21        | 5            | 21             |
    | 1005a     | 13        | 4            | 13             |
    | 1006      | 20        | 10           | 20             |
    | 1007      | 14        | 7            | 14             |
    | 1008      | 36        | 9            | 34             |
    | 1009      | 107       | 5            | 34             |
    | 1011      | 51        | 17           | 48             |
    | 3001      | 112       | 6            | 112            |
    | 3002      | 10        | 10           | 10             |
    | 3003      | 479       | 19           | 414            |
    information about sex ratio
    the datasets contained 691 males and 84 females
    that corresponds to 10.84% females 
    information about animal handling
    11 datasets used headplates
    8 datasets used body restraining
    4 datasets used anesthesia before acquisition
    8 datasets were collected by males, 2 by women
    | rodent.ds | headplate | restrained | anesthesia | exp.gender |
    |-----------|-----------|------------|------------|------------|
    | 1001      | y         | n          | n          | m          |
    | 1002      | y         | n          | n          | m          |
    | 1003      | y         | y          | n          | m          |
    | 1004      | y         | y          | n          | f          |
    | 1005a     | y         | y          | Isoflurane | f          |
    | 1006      | y         | n          | Isoflurane | m          |
    | 1007      | y         | y          | n          | m          |
    | 1008      | n         | y          | n          | m          |
    | 1009      | y         | y          | Isoflurane | m          |
    | 1011      | y         | y          | Isoflurane | m          |
    | 3001      | y         | null       | null       | null       |
    | 3002      | y         | y          | n          | null       |
    | 3003      | null      | null       | null       | null       |
    information about the scanner and sequence
    lowest field strength was 7.0T
    highest field strength was 15.2T
    | rodent.ds | field_strength | sequence | TE     |
    |-----------|----------------|----------|--------|
    | 1001      | 9.4            | GE-EPI   | 15.0   |
    | 1002      | 7.0            | GE-EPI   | 16.385 |
    | 1003      | 7.0            | GE-EPI   | 15.0   |
    | 1004      | 11.7           | GE-EPI   | 11.0   |
    | 1005a     | 9.4            | SE-EPI   | 12.5   |
    | 1006      | 11.7           | GE-EPI   | 0.014  |
    | 1007      | 9.4            | GE-EPI   | 12.0   |
    | 1008      | 9.4            | GE-EPI   | 12.0   |
    | 1009      | 9.4            | GE-EPI   | 15.0   |
    | 1011      | 15.2           | GE-EPI   | 11.7   |
    | 3001      | 9.4            | SE-EPI   | 18.398 |
    | 3002      | 7.0            | GE-EPI   | 15.0   |
    | 3003      | 9.4            | SE-EPI   | 18.398 |
    #### MOTION ANALYISIS ####
    mean fd across all mouse datasets
    | fd.mean  |
    |----------|
    | 0.036873 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 1001      | 0.007517 |
    | 1002      | 0.032146 |
    | 1003      | 0.011058 |
    | 1004      | 0.065948 |
    | 1005a     | 0.010155 |
    | 1006      | 0.06656  |
    | 1007      | 0.016836 |
    | 1008      | 0.033217 |
    | 1009      | 0.03189  |
    | 1011      | 0.011393 |
    | 3001      | 0.035982 |
    | 3002      | 0.022321 |
    | 3003      | 0.040797 |
    Number of dropped frames for each dataset and denoising method
    This corresponds to the following rabies flags for mice
    #gsr1: --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true
    #gsr2: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true
    #gsr3: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false
    #### DENOISE ANALYSIS ####
    dropped frames per dataset
    | rodent.ds | total.frames | dropped.frames.gsr1 | dropped.frames.gsr2 | dropped.frames.gsr3 |
    |-----------|--------------|---------------------|---------------------|---------------------|
    | 1001      | 600.0        | 42.25               | 42.25               | null                |
    | 1002      | 900.0        | 112.5               | 75.125              | null                |
    | 1003      | 1800.0       | 157.0               | 157.0               | null                |
    | 1004      | 334.0        | 129.238095          | 59.095238           | 5.9                 |
    | 1005a     | 600.0        | 65.384615           | 62.461538           | null                |
    | 1006      | 618.0        | 216.55              | 77.5                | 18.866667           |
    | 1007      | 600.0        | 39.714286           | 39.714286           | null                |
    | 1008      | 400.0        | 83.823529           | 58.764706           | null                |
    | 1009      | 397.058824   | 69.323529           | 52.617647           | 5.111111            |
    | 1011      | 360.0        | 40.208333           | 40.208333           | null                |
    | 3001      | 200.0        | 38.142857           | 28.821429           | 4.404255            |
    | 3002      | 1920.0       | 400.2               | 388.0               | 3.0                 |
    | 3003      | 200.0        | 55.074879           | 41.342995           | 4.953608            |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Spurious     | 0.473118 |
    | No           | 0.233871 |
    | Specific     | 0.19086  |
    | Non-specific | 0.102151 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Spurious     | 0.461022 |
    | No           | 0.251344 |
    | Specific     | 0.200269 |
    | Non-specific | 0.087366 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Spurious     | 0.443548 |
    | No           | 0.259409 |
    | Specific     | 0.211022 |
    | Non-specific | 0.086022 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.466398 |
    | No            | 0.217742 |
    | Specific      | 0.186828 |
    | Non-specific  | 0.129032 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.451613 |
    | No            | 0.245968 |
    | Specific      | 0.194892 |
    | Non-specific  | 0.107527 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Spurious      | 0.431452 |
    | No            | 0.25     |
    | Specific      | 0.209677 |
    | Non-specific  | 0.108871 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.479839 |
    | Specific         | 0.225806 |
    | No               | 0.186828 |
    | Non-specific     | 0.107527 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.461022 |
    | Specific         | 0.228495 |
    | No               | 0.193548 |
    | Non-specific     | 0.116935 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Spurious         | 0.407258 |
    | No               | 0.251344 |
    | Specific         | 0.229839 |
    | Non-specific     | 0.111559 |
    for gsr processing, heavy scrubbing lead to  142/744 Specific scans
    for gsr processing, light scrubbing lead to  157/744 Specific scans
    for aCompCor processing, heavy scrubbing lead to  168/744 Specific scans
    for aCompCor processing, light scrubbing lead to  171/744 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 1001      | 0.5       | 0.5       | 0.75      | … | 0.5       | 0.5       | 0.5       | 0.5      |
    | 1002      | 0.5       | 0.375     | 0.625     | … | 0.75      | 0.75      | 0.625     | 0.875    |
    | 1003      | 0.0       | 0.0       | 0.083333  | … | 0.083333  | 0.083333  | 0.083333  | 0.083333 |
    | 1004      | 0.047619  | 0.142857  | 0.238095  | … | 0.142857  | 0.047619  | 0.095238  | 0.095238 |
    | 1005a     | 0.307692  | 0.307692  | 0.076923  | … | 0.153846  | 0.461538  | 0.461538  | 0.307692 |
    | 1006      | 0.2       | 0.3       | 0.15      | … | 0.4       | 0.3       | 0.4       | 0.45     |
    | 1007      | 0.142857  | 0.142857  | 0.0       | … | 0.0       | 0.142857  | 0.142857  | 0.142857 |
    | 1008      | 0.235294  | 0.235294  | 0.382353  | … | 0.323529  | 0.205882  | 0.235294  | 0.117647 |
    | 1009      | 0.235294  | 0.264706  | 0.235294  | … | 0.235294  | 0.176471  | 0.235294  | 0.352941 |
    | 1011      | 0.166667  | 0.166667  | 0.291667  | … | 0.208333  | 0.208333  | 0.270833  | 0.291667 |
    | 3001      | 0.303571  | 0.3125    | 0.267857  | … | 0.214286  | 0.303571  | 0.294643  | 0.223214 |
    | 3002      | 0.5       | 0.5       | 0.7       | … | 0.7       | 0.1       | 0.0       | 0.1      |
    | 3003      | 0.149758  | 0.154589  | 0.161836  | … | 0.178744  | 0.207729  | 0.198068  | 0.21256  |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 8/13 dataset performed better with less scrubbing
    for aCompCor processing, 6/13 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 6/13 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 5/13 dataset performed worst with aCompCor compared to gsr
    #### NOW DOING rat ####
    summary of the data that we collected
    we processed 5 datasets
    totalling 421 runs
    from 130 mice
    the smallest dataset had 7 runs
    the largest dataset had 291 runs
    we could processed 190/421 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included |
    |-----------|-----------|--------------|----------------|
    | 2001      | 16        | 16           | 16             |
    | 2002a     | 7         | 7            | 7              |
    | 2002b     | 10        | 10           | 10             |
    | 2003      | 97        | 8            | 97             |
    | 4001      | 291       | 89           | 60             |
    information about sex ratio
    the datasets contained 413 males and 8 females
    that corresponds to 1.9% females 
    information about animal handling
    1 datasets used headplates
    5 datasets used body restraining
    5 datasets used anesthesia before acquisition
    1 datasets were collected by males, 3 by women
    | rodent.ds | headplate | restrained | anesthesia  | exp.gender |
    |-----------|-----------|------------|-------------|------------|
    | 2001      | n         | y          | Isoflurane  | f          |
    | 2002a     | n         | y          | Sevoflurane | f          |
    | 2002b     | n         | y          | Isoflurane  | m          |
    | 2003      | y         | y          | Isoflurane  | f          |
    | 4001      | n         | y          | Isoflurane  | null       |
    information about the scanner and sequence
    lowest field strength was 7T
    highest field strength was 7T
    | rodent.ds | field_strength | sequence | TE |
    |-----------|----------------|----------|----|
    | 2001      | 7              | GE-EPI   | 17 |
    | 2002a     | 7              | GE-EPI   | 18 |
    | 2002b     | 7              | SE-EPI   | 45 |
    | 2003      | 7              | GE-EPI   | 15 |
    | 4001      | 7              | GE-EPI   | 15 |
    #### MOTION ANALYISIS ####
    mean fd across all rat datasets
    | fd.mean  |
    |----------|
    | 0.044321 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 2001      | 0.03643  |
    | 2002a     | 0.052678 |
    | 2002b     | 0.044433 |
    | 2003      | 0.04796  |
    | 4001      | 0.039548 |
    Number of dropped frames for each dataset and denoising method
    This corresponds to the following rabies flags for mice
    #gsr1: --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true
    #gsr2: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true
    #gsr3: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false
    #### DENOISE ANALYSIS ####
    dropped frames per dataset
    | rodent.ds | total.frames | dropped.frames.gsr1 | dropped.frames.gsr2 | dropped.frames.gsr3 |
    |-----------|--------------|---------------------|---------------------|---------------------|
    | 2001      | 600.0        | 52.8125             | 20.9375             | 7.0                 |
    | 2002a     | 1500.0       | 426.714286          | 267.142857          | 5.666667            |
    | 2002b     | 750.0        | 143.1               | 97.0                | 6.0                 |
    | 2003      | 90.0         | 21.092784           | 9.185567            | 10.5                |
    | 4001      | 600.0        | 103.85              | 79.716667           | 13.151515           |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Specific     | 0.410526 |
    | Spurious     | 0.226316 |
    | Non-specific | 0.226316 |
    | No           | 0.136842 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Specific     | 0.410526 |
    | Spurious     | 0.231579 |
    | Non-specific | 0.221053 |
    | No           | 0.136842 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Specific     | 0.436842 |
    | Spurious     | 0.257895 |
    | Non-specific | 0.210526 |
    | No           | 0.094737 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Non-specific  | 0.326316 |
    | Specific      | 0.310526 |
    | Spurious      | 0.247368 |
    | No            | 0.115789 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Specific      | 0.336842 |
    | Non-specific  | 0.278947 |
    | Spurious      | 0.257895 |
    | No            | 0.126316 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Specific      | 0.368421 |
    | Non-specific  | 0.305263 |
    | Spurious      | 0.231579 |
    | No            | 0.094737 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Specific         | 0.363158 |
    | Non-specific     | 0.310526 |
    | Spurious         | 0.205263 |
    | No               | 0.121053 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Specific         | 0.357895 |
    | Non-specific     | 0.331579 |
    | Spurious         | 0.194737 |
    | No               | 0.115789 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Specific         | 0.384211 |
    | Non-specific     | 0.326316 |
    | Spurious         | 0.178947 |
    | No               | 0.110526 |
    for gsr processing, heavy scrubbing lead to  78/190 Specific scans
    for gsr processing, light scrubbing lead to  83/190 Specific scans
    for aCompCor processing, heavy scrubbing lead to  69/190 Specific scans
    for aCompCor processing, light scrubbing lead to  73/190 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 2001      | 0.625     | 0.625     | 0.625     | … | 0.4375    | 0.625     | 0.5       | 0.5625   |
    | 2002a     | 0.714286  | 0.571429  | 0.285714  | … | 0.142857  | 0.285714  | 0.285714  | 0.285714 |
    | 2002b     | 0.3       | 0.4       | 0.3       | … | 0.3       | 0.3       | 0.5       | 0.4      |
    | 2003      | 0.402062  | 0.391753  | 0.463918  | … | 0.391753  | 0.329897  | 0.298969  | 0.350515 |
    | 4001      | 0.35      | 0.366667  | 0.383333  | … | 0.35      | 0.366667  | 0.4       | 0.4      |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 2/5 dataset performed better with less scrubbing
    for aCompCor processing, 3/5 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 2/5 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 2/5 dataset performed worst with aCompCor compared to gsr
