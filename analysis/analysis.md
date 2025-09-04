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

from scipy.stats import chi2_contingency

def chi2_continuous(df, cat1, cat2, rodent, split=4, labels = ["lowest","low","high","highest"]):
  print(f'looking at {cat2} effect in {cat1} in {rodent}')
  chi2_table = df.select(cat1, cat2).with_columns(pl.col(cat2).qcut(split, labels=labels).alias('quartiles')).drop_nulls()
  total_count = chi2_table.height/100
  chi2_table = chi2_table.select(cat1,'quartiles').group_by([cat1,'quartiles']).agg(pl.len()).sort(by='quartiles').pivot(index=cat1, on='quartiles', values='len').sort(by=cat1)
  print(chi2_table.with_columns(pl.all().exclude(chi2_table.columns[0])/total_count))
  chi2_table = chi2_table.select(pl.all().exclude(chi2_table.columns[0]))
  q, p, dof, expect = chi2_contingency(chi2_table)
  print(f'the effect of {cat2} on {cat1} in {rodent} is q =  {round(q,2)} with p-value = {round(p,5)}, dof = {dof}')

def chi2_categorical(df, cat1, cat2, rodent):
  print(f'looking at {cat2} effect in {cat1} in {rodent}')
  chi2_table = df.select(cat1, cat2).drop_nulls()
  total_count = chi2_table.height/100
  chi2_table = chi2_table.select(cat1,cat2).group_by([cat1,cat2]).agg(pl.len()).sort(by=cat2).pivot(index=cat1, on=cat2, values='len').sort(by=cat1)
  print(chi2_table.with_columns(pl.all().exclude(chi2_table.columns[0])/total_count))
  chi2_table = chi2_table.select(pl.all().exclude(chi2_table.columns[0]))
  q, p, dof, expect = chi2_contingency(chi2_table)
  print(f'the effect of {cat2} on {cat1} in {rodent} is q =  {round(q,2)} with p-value = {round(p,5)}, dof = {dof}')

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
#rodent='mouse'
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
  print(str((df_summary["exp.gender"] == 'm').sum()) + " datasets were collected by men, " + str((df_summary["exp.gender"] == 'f').sum())+ " by women")

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

  print('doing statistical analysis')
  print('first with gsr processed data')
  chi2_continuous(df, 's1.cat.gsr3', 'fd.mean', rodent)
  chi2_continuous(df, 's1.cat.gsr3', 's1.gsrcov.l.gsr3', rodent)
  chi2_continuous(df, 's1.cat.gsr3', 's1.tsnr.l', rodent)

  chi2_categorical(df, 's1.cat.gsr3', 'main.experimenter.gender', rodent)
  chi2_categorical(df, 's1.cat.gsr3', 'rodent.sex', rodent)
  chi2_categorical(df, 's1.cat.gsr3', 'head-plate', rodent)
  chi2_categorical(df, 's1.cat.gsr3', 'body.restrained', rodent)

  print('second with aCompCor processed data') 
  chi2_continuous(df, 's1.cat.aCompCor3', 'fd.mean', rodent)
  chi2_continuous(df, 's1.cat.aCompCor3', 's1.gsrcov.l.gsr3', rodent)
  chi2_continuous(df, 's1.cat.aCompCor3', 's1.tsnr.l', rodent)

  chi2_categorical(df, 's1.cat.aCompCor3', 'main.experimenter.gender', rodent)
  chi2_categorical(df, 's1.cat.aCompCor3', 'rodent.sex', rodent)
  chi2_categorical(df, 's1.cat.aCompCor3', 'head-plate', rodent)
  chi2_categorical(df, 's1.cat.aCompCor3', 'body.restrained', rodent)
```

    #### NOW DOING mouse ####
    summary of the data that we collected
    we processed 18 datasets
    totalling 1338 runs
    from 171 mice
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 1151/1338 runs.
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
    | 1012      | 26        | 26           | 21             |
    | 1013      | 12        | 3            | 12             |
    | 1014      | 4         | 4            | 4              |
    | 3001      | 112       | 6            | 112            |
    | 3002      | 10        | 10           | 10             |
    | 3003      | 479       | 19           | 414            |
    | 3004      | 54        | 9            | 54             |
    | 3005      | 355       | 13           | 316            |
    information about sex ratio
    the datasets contained 1076 males and 150 females
    that corresponds to 12.23% females 
    information about animal handling
    14 datasets used headplates
    10 datasets used body restraining
    6 datasets used anesthesia before acquisition
    10 datasets were collected by men, 8 by women
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
    | 1012      | n         | null       | Isoflurane | f          |
    | 1013      | y         | n          | no         | f          |
    | 1014      | y         | y          | null       | f          |
    | 3001      | y         | null       | null       | m          |
    | 3002      | y         | y          | n          | f          |
    | 3003      | null      | null       | null       | m          |
    | 3004      | n         | y          | Isoflurane | f          |
    | 3005      | y         | n          | null       | f          |
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
    | 1012      | 9.4            | GE-EPI   | 6.0    |
    | 1013      | 9.4            | GE-EPI   | 14.0   |
    | 1014      | 7.0            | GE-EPI   | 17.0   |
    | 3001      | 9.4            | SE-EPI   | 18.398 |
    | 3002      | 7.0            | GE-EPI   | 15.0   |
    | 3003      | 9.4            | SE-EPI   | 18.398 |
    | 3004      | 11.7           | GE-EPI   | 15.0   |
    | 3005      | 9.4            | SE-EPI   | null   |
    #### MOTION ANALYISIS ####
    mean fd across all mouse datasets
    | fd.mean  |
    |----------|
    | 0.037095 |
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
    | 1012      | 0.017715 |
    | 1013      | 0.028279 |
    | 1014      | 0.005115 |
    | 3001      | 0.035982 |
    | 3002      | 0.022321 |
    | 3003      | 0.040797 |
    | 3004      | 0.03515  |
    | 3005      | 0.046262 |
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
    | 1012      | 300.0        | 26.714286           | 25.428571           | 4.333333            |
    | 1013      | 480.0        | 61.333333           | 56.333333           | null                |
    | 1014      | 150.0        | 13.0                | 12.75               | null                |
    | 3001      | 200.0        | 38.142857           | 28.821429           | 4.404255            |
    | 3002      | 1920.0       | 400.2               | 388.0               | 3.0                 |
    | 3003      | 200.0        | 55.074879           | 41.342995           | 4.953608            |
    | 3004      | null         | 51.296296           | 39.150943           | 7.75                |
    | 3005      | 200.0        | 52.234177           | 38.547468           | 4.965986            |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Spurious     | 0.454764 |
    | No           | 0.239392 |
    | Specific     | 0.21297  |
    | Non-specific | 0.092874 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Spurious     | 0.445156 |
    | No           | 0.252202 |
    | Specific     | 0.218575 |
    | Non-specific | 0.084067 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Spurious     | 0.423539 |
    | No           | 0.273018 |
    | Specific     | 0.224179 |
    | Non-specific | 0.079263 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.443555 |
    | No            | 0.238591 |
    | Specific      | 0.204163 |
    | Non-specific  | 0.113691 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.436349 |
    | No            | 0.254604 |
    | Specific      | 0.210568 |
    | Non-specific  | 0.098479 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Spurious      | 0.397118 |
    | No            | 0.265813 |
    | Specific      | 0.225781 |
    | Non-specific  | 0.111289 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.417134 |
    | No               | 0.232986 |
    | Specific         | 0.232186 |
    | Non-specific     | 0.117694 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.413931 |
    | No               | 0.236189 |
    | Specific         | 0.231385 |
    | Non-specific     | 0.118495 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Spurious         | 0.359488 |
    | No               | 0.292234 |
    | Specific         | 0.233787 |
    | Non-specific     | 0.114492 |
    for gsr processing, heavy scrubbing lead to  266/1249 Specific scans
    for gsr processing, light scrubbing lead to  280/1249 Specific scans
    for aCompCor processing, heavy scrubbing lead to  290/1249 Specific scans
    for aCompCor processing, light scrubbing lead to  292/1249 Specific scans
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
    | 1012      | 0.333333  | 0.285714  | 0.333333  | … | 0.380952  | 0.428571  | 0.428571  | 0.380952 |
    | 1013      | 0.166667  | 0.166667  | 0.25      | … | 0.166667  | 0.083333  | 0.166667  | 0.083333 |
    | 1014      | 0.25      | 0.25      | 0.0       | … | 0.0       | 0.75      | 0.5       | 0.0      |
    | 3001      | 0.303571  | 0.3125    | 0.267857  | … | 0.214286  | 0.303571  | 0.294643  | 0.223214 |
    | 3002      | 0.5       | 0.5       | 0.7       | … | 0.7       | 0.1       | 0.0       | 0.1      |
    | 3003      | 0.149758  | 0.154589  | 0.161836  | … | 0.178744  | 0.207729  | 0.198068  | 0.21256  |
    | 3004      | 0.592593  | 0.62963   | 0.611111  | … | 0.555556  | 0.407407  | 0.388889  | 0.37037  |
    | 3005      | 0.21519   | 0.212025  | 0.253165  | … | 0.272152  | 0.231013  | 0.224684  | 0.246835 |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 11/18 dataset performed better with less scrubbing
    for aCompCor processing, 7/18 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 7/18 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 8/18 dataset performed worst with aCompCor compared to gsr
    doing statistical analysis
    first with gsr processed data
    looking at fd.mean effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | lowest    | low      | high      | highest  |
    |--------------|-----------|----------|-----------|----------|
    | No           | 6.810897  | 6.971154 | 7.211538  | 6.330128 |
    | Non-specific | 1.923077  | 1.842949 | 1.442308  | 2.724359 |
    | Specific     | 4.166667  | 6.25     | 5.128205  | 6.891026 |
    | Spurious     | 12.099359 | 9.935897 | 11.217949 | 9.054487 |
    the effect of fd.mean on s1.cat.gsr3 in mouse is q =  22.36 with p-value = 0.00781, dof = 9
    looking at s1.gsrcov.l.gsr3 effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | lowest    | low       | high     | highest   |
    |--------------|-----------|-----------|----------|-----------|
    | No           | 7.234727  | 6.189711  | 6.028939 | 7.877814  |
    | Non-specific | 1.768489  | 1.768489  | 2.491961 | 1.92926   |
    | Specific     | 5.707395  | 4.581994  | 6.993569 | 5.144695  |
    | Spurious     | 10.289389 | 12.459807 | 9.485531 | 10.048232 |
    the effect of s1.gsrcov.l.gsr3 on s1.cat.gsr3 in mouse is q =  19.52 with p-value = 0.02114, dof = 9
    looking at s1.tsnr.l effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | lowest   | low       | high     | highest  |
    |--------------|----------|-----------|----------|----------|
    | No           | 7.371795 | 6.730769  | 8.573718 | 4.647436 |
    | Non-specific | 1.682692 | 1.842949  | 2.163462 | 2.24359  |
    | Specific     | 7.13141  | 4.246795  | 4.567308 | 6.490385 |
    | Spurious     | 8.814103 | 12.179487 | 9.695513 | 11.61859 |
    the effect of s1.tsnr.l on s1.cat.gsr3 in mouse is q =  38.46 with p-value = 1e-05, dof = 9
    looking at main.experimenter.gender effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | f         | m         |
    |--------------|-----------|-----------|
    | No           | 10.08807  | 17.213771 |
    | Non-specific | 3.202562  | 4.723779  |
    | Specific     | 10.888711 | 11.529223 |
    | Spurious     | 11.929544 | 30.424339 |
    the effect of main.experimenter.gender on s1.cat.gsr3 in mouse is q =  34.21 with p-value = 0.0, dof = 3
    looking at rodent.sex effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | f        | m         |
    |--------------|----------|-----------|
    | No           | 2.814424 | 27.00088  |
    | Non-specific | 1.671064 | 5.452946  |
    | Specific     | 4.573439 | 17.414248 |
    | Spurious     | 5.101143 | 35.971856 |
    the effect of rodent.sex on s1.cat.gsr3 in mouse is q =  22.21 with p-value = 6e-05, dof = 3
    looking at head-plate effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | n        | y         |
    |--------------|----------|-----------|
    | No           | 2.155689 | 21.676647 |
    | Non-specific | 1.796407 | 7.54491   |
    | Specific     | 6.347305 | 19.161677 |
    | Spurious     | 2.754491 | 38.562874 |
    the effect of head-plate on s1.cat.gsr3 in mouse is q =  44.1 with p-value = 0.0, dof = 3
    looking at body.restrained effect in s1.cat.gsr3 in mouse
    | s1.cat.gsr3  | n         | y         |
    |--------------|-----------|-----------|
    | No           | 15.954416 | 11.111111 |
    | Non-specific | 4.558405  | 3.846154  |
    | Specific     | 13.390313 | 11.680912 |
    | Spurious     | 17.378917 | 22.079772 |
    the effect of body.restrained on s1.cat.gsr3 in mouse is q =  10.8 with p-value = 0.01284, dof = 3
    second with aCompCor processed data
    looking at fd.mean effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | lowest    | low      | high      | highest  |
    |------------------|-----------|----------|-----------|----------|
    | No               | 12.259615 | 5.288462 | 5.849359  | 5.849359 |
    | Non-specific     | 2.564103  | 3.365385 | 2.483974  | 3.044872 |
    | Specific         | 4.887821  | 6.730769 | 5.208333  | 6.570513 |
    | Spurious         | 5.288462  | 9.615385 | 11.458333 | 9.535256 |
    the effect of fd.mean on s1.cat.aCompCor3 in mouse is q =  92.43 with p-value = 0.0, dof = 9
    looking at s1.gsrcov.l.gsr3 effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | lowest   | low       | high     | highest   |
    |------------------|----------|-----------|----------|-----------|
    | No               | 7.073955 | 10.048232 | 5.22508  | 6.993569  |
    | Non-specific     | 2.652733 | 3.054662  | 3.456592 | 2.33119   |
    | Specific         | 5.546624 | 5.546624  | 6.913183 | 5.305466  |
    | Spurious         | 9.726688 | 6.350482  | 9.405145 | 10.369775 |
    the effect of s1.gsrcov.l.gsr3 on s1.cat.aCompCor3 in mouse is q =  40.18 with p-value = 1e-05, dof = 9
    looking at s1.tsnr.l effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | lowest   | low       | high      | highest   |
    |------------------|----------|-----------|-----------|-----------|
    | No               | 5.528846 | 5.048077  | 8.012821  | 10.657051 |
    | Non-specific     | 3.205128 | 2.644231  | 2.083333  | 3.525641  |
    | Specific         | 6.410256 | 5.448718  | 4.887821  | 6.650641  |
    | Spurious         | 9.855769 | 11.858974 | 10.016026 | 4.166667  |
    the effect of s1.tsnr.l on s1.cat.aCompCor3 in mouse is q =  90.05 with p-value = 0.0, dof = 9
    looking at main.experimenter.gender effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | f         | m         |
    |------------------|-----------|-----------|
    | No               | 8.566853  | 20.656525 |
    | Non-specific     | 5.044035  | 6.405124  |
    | Specific         | 9.127302  | 14.251401 |
    | Spurious         | 13.370697 | 22.578062 |
    the effect of main.experimenter.gender on s1.cat.aCompCor3 in mouse is q =  12.53 with p-value = 0.00576, dof = 3
    looking at rodent.sex effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | f        | m         |
    |------------------|----------|-----------|
    | No               | 5.189094 | 26.649077 |
    | Non-specific     | 2.902375 | 6.772208  |
    | Specific         | 4.485488 | 18.997361 |
    | Spurious         | 1.583113 | 33.421284 |
    the effect of rodent.sex on s1.cat.aCompCor3 in mouse is q =  59.84 with p-value = 0.0, dof = 3
    looking at head-plate effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | n        | y         |
    |------------------|----------|-----------|
    | No               | 3.353293 | 27.305389 |
    | Non-specific     | 3.832335 | 9.341317  |
    | Specific         | 3.832335 | 20.598802 |
    | Spurious         | 2.035928 | 29.700599 |
    the effect of head-plate on s1.cat.aCompCor3 in mouse is q =  37.47 with p-value = 0.0, dof = 3
    looking at body.restrained effect in s1.cat.aCompCor3 in mouse
    | s1.cat.aCompCor3 | n         | y         |
    |------------------|-----------|-----------|
    | No               | 12.678063 | 22.364672 |
    | Non-specific     | 4.558405  | 6.125356  |
    | Specific         | 13.817664 | 10.541311 |
    | Spurious         | 20.22792  | 9.68661   |
    the effect of body.restrained on s1.cat.aCompCor3 in mouse is q =  49.15 with p-value = 0.0, dof = 3
    #### NOW DOING rat ####
    summary of the data that we collected
    we processed 7 datasets
    totalling 445 runs
    from 138 mice
    the smallest dataset had 5 runs
    the largest dataset had 291 runs
    we could processed 212/445 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included |
    |-----------|-----------|--------------|----------------|
    | 2001      | 16        | 16           | 16             |
    | 2002a     | 7         | 7            | 7              |
    | 2002b     | 10        | 10           | 10             |
    | 2003      | 97        | 8            | 97             |
    | 2004      | 19        | 5            | 19             |
    | 2005      | 5         | 3            | 3              |
    | 4001      | 291       | 89           | 60             |
    information about sex ratio
    the datasets contained 423 males and 22 females
    that corresponds to 4.94% females 
    information about animal handling
    2 datasets used headplates
    7 datasets used body restraining
    6 datasets used anesthesia before acquisition
    3 datasets were collected by men, 3 by women
    | rodent.ds | headplate | restrained | anesthesia  | exp.gender |
    |-----------|-----------|------------|-------------|------------|
    | 2001      | n         | y          | Isoflurane  | f          |
    | 2002a     | n         | y          | Sevoflurane | f          |
    | 2002b     | n         | y          | Isoflurane  | m          |
    | 2003      | y         | y          | Isoflurane  | f          |
    | 2004      | y         | y          | n           | m          |
    | 2005      | n         | y          | Isoflurane  | m          |
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
    | 2004      | null           | GE-EPI   | 15 |
    | 2005      | 7              | GE-EPI   | 17 |
    | 4001      | 7              | GE-EPI   | 15 |
    #### MOTION ANALYISIS ####
    mean fd across all rat datasets
    | fd.mean  |
    |----------|
    | 0.041551 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 2001      | 0.03643  |
    | 2002a     | 0.052678 |
    | 2002b     | 0.044433 |
    | 2003      | 0.04796  |
    | 2004      | 0.016877 |
    | 2005      | 0.022399 |
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
    | 2004      | 400.0        | 80.105263           | 79.157895           | null                |
    | 2005      | 2400.0       | 285.666667          | 179.666667          | 3.0                 |
    | 4001      | 600.0        | 103.85              | 79.716667           | 13.151515           |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Specific     | 0.391509 |
    | Spurious     | 0.235849 |
    | Non-specific | 0.207547 |
    | No           | 0.160377 |
    | null         | 0.004717 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Specific     | 0.391509 |
    | Spurious     | 0.240566 |
    | Non-specific | 0.20283  |
    | No           | 0.160377 |
    | null         | 0.004717 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Specific     | 0.410377 |
    | Spurious     | 0.259434 |
    | Non-specific | 0.198113 |
    | No           | 0.127358 |
    | null         | 0.004717 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Non-specific  | 0.316038 |
    | Specific      | 0.292453 |
    | Spurious      | 0.259434 |
    | No            | 0.132075 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Specific      | 0.316038 |
    | Non-specific  | 0.268868 |
    | Spurious      | 0.268868 |
    | No            | 0.141509 |
    | null          | 0.004717 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Specific      | 0.353774 |
    | Non-specific  | 0.287736 |
    | Spurious      | 0.240566 |
    | No            | 0.113208 |
    | null          | 0.004717 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Specific         | 0.334906 |
    | Non-specific     | 0.287736 |
    | Spurious         | 0.235849 |
    | No               | 0.136792 |
    | null             | 0.004717 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Specific         | 0.334906 |
    | Non-specific     | 0.311321 |
    | Spurious         | 0.221698 |
    | No               | 0.127358 |
    | null             | 0.004717 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Specific         | 0.363208 |
    | Non-specific     | 0.29717  |
    | Spurious         | 0.193396 |
    | No               | 0.141509 |
    | null             | 0.004717 |
    for gsr processing, heavy scrubbing lead to  83/211 Specific scans
    for gsr processing, light scrubbing lead to  87/211 Specific scans
    for aCompCor processing, heavy scrubbing lead to  71/211 Specific scans
    for aCompCor processing, light scrubbing lead to  77/211 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 2001      | 0.625     | 0.625     | 0.625     | … | 0.4375    | 0.625     | 0.5       | 0.5625   |
    | 2002a     | 0.714286  | 0.571429  | 0.285714  | … | 0.142857  | 0.285714  | 0.285714  | 0.285714 |
    | 2002b     | 0.3       | 0.4       | 0.3       | … | 0.3       | 0.3       | 0.5       | 0.4      |
    | 2003      | 0.402062  | 0.391753  | 0.463918  | … | 0.391753  | 0.329897  | 0.298969  | 0.350515 |
    | 2004      | 0.210526  | 0.210526  | 0.210526  | … | 0.263158  | 0.105263  | 0.157895  | 0.210526 |
    | 2005      | 0.333333  | 0.333333  | 0.0       | … | 0.0       | 0.0       | 0.0       | 0.0      |
    | 4001      | 0.35      | 0.366667  | 0.383333  | … | 0.35      | 0.366667  | 0.4       | 0.4      |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 2/7 dataset performed better with less scrubbing
    for aCompCor processing, 4/7 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 2/7 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 2/7 dataset performed worst with aCompCor compared to gsr
    doing statistical analysis
    first with gsr processed data
    looking at fd.mean effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | lowest   | low       | high     | highest  |
    |--------------|----------|-----------|----------|----------|
    | No           | 6.161137 | 2.369668  | 2.843602 | 1.421801 |
    | Non-specific | 2.369668 | 5.687204  | 6.161137 | 5.687204 |
    | Specific     | 8.056872 | 12.796209 | 10.42654 | 9.952607 |
    | Spurious     | 8.530806 | 3.791469  | 5.687204 | 8.056872 |
    the effect of fd.mean on s1.cat.gsr3 in rat is q =  19.33 with p-value = 0.02251, dof = 9
    looking at s1.gsrcov.l.gsr3 effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | lowest   | low      | high      | highest   |
    |--------------|----------|----------|-----------|-----------|
    | No           | 1.421801 | 6.635071 | 3.791469  | 0.947867  |
    | Non-specific | 7.582938 | 2.843602 | 4.739336  | 4.739336  |
    | Specific     | 10.42654 | 7.109005 | 10.900474 | 12.796209 |
    | Spurious     | 5.21327  | 8.530806 | 5.687204  | 6.635071  |
    the effect of s1.gsrcov.l.gsr3 on s1.cat.gsr3 in rat is q =  23.81 with p-value = 0.00462, dof = 9
    looking at s1.tsnr.l effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | lowest   | low      | high      | highest   |
    |--------------|----------|----------|-----------|-----------|
    | No           | 4.739336 | 2.843602 | 1.895735  | 3.317536  |
    | Non-specific | 3.791469 | 6.635071 | 6.161137  | 3.317536  |
    | Specific     | 7.582938 | 10.42654 | 12.322275 | 10.900474 |
    | Spurious     | 8.530806 | 5.21327  | 4.739336  | 7.582938  |
    the effect of s1.tsnr.l on s1.cat.gsr3 in rat is q =  12.01 with p-value = 0.21296, dof = 9
    looking at main.experimenter.gender effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | f         | m        |
    |--------------|-----------|----------|
    | No           | 3.97351   | 6.622517 |
    | Non-specific | 19.86755  | 3.97351  |
    | Specific     | 37.748344 | 4.635762 |
    | Spurious     | 17.880795 | 5.298013 |
    the effect of main.experimenter.gender on s1.cat.gsr3 in rat is q =  21.33 with p-value = 9e-05, dof = 3
    looking at rodent.sex effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | f        | m         |
    |--------------|----------|-----------|
    | No           | 2.369668 | 10.42654  |
    | Non-specific | 1.421801 | 18.483412 |
    | Specific     | 3.791469 | 37.440758 |
    | Spurious     | 1.421801 | 24.64455  |
    the effect of rodent.sex on s1.cat.gsr3 in rat is q =  4.01 with p-value = 0.26038, dof = 3
    looking at head-plate effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | n         | y         |
    |--------------|-----------|-----------|
    | No           | 6.635071  | 6.161137  |
    | Non-specific | 8.530806  | 11.374408 |
    | Specific     | 18.009479 | 23.222749 |
    | Spurious     | 11.848341 | 14.218009 |
    the effect of head-plate on s1.cat.gsr3 in rat is q =  0.66 with p-value = 0.88351, dof = 3
    looking at body.restrained effect in s1.cat.gsr3 in rat
    | s1.cat.gsr3  | y         |
    |--------------|-----------|
    | No           | 12.796209 |
    | Non-specific | 19.905213 |
    | Specific     | 41.232227 |
    | Spurious     | 26.066351 |
    the effect of body.restrained on s1.cat.gsr3 in rat is q =  0.0 with p-value = 1.0, dof = 0
    second with aCompCor processed data
    looking at fd.mean effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | lowest   | low       | high      | highest  |
    |------------------|----------|-----------|-----------|----------|
    | No               | 6.635071 | 4.265403  | 1.421801  | 1.895735 |
    | Non-specific     | 2.369668 | 7.109005  | 10.900474 | 9.478673 |
    | Specific         | 7.582938 | 11.374408 | 8.530806  | 9.004739 |
    | Spurious         | 8.530806 | 1.895735  | 4.265403  | 4.739336 |
    the effect of fd.mean on s1.cat.aCompCor3 in rat is q =  33.69 with p-value = 0.0001, dof = 9
    looking at s1.gsrcov.l.gsr3 effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | lowest   | low      | high     | highest   |
    |------------------|----------|----------|----------|-----------|
    | No               | 1.421801 | 6.161137 | 5.21327  | 1.421801  |
    | Non-specific     | 10.42654 | 3.317536 | 7.582938 | 8.530806  |
    | Specific         | 9.004739 | 8.056872 | 8.530806 | 10.900474 |
    | Spurious         | 3.791469 | 7.582938 | 3.791469 | 4.265403  |
    the effect of s1.gsrcov.l.gsr3 on s1.cat.aCompCor3 in rat is q =  24.16 with p-value = 0.00406, dof = 9
    looking at s1.tsnr.l effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | lowest   | low       | high      | highest  |
    |------------------|----------|-----------|-----------|----------|
    | No               | 6.161137 | 1.895735  | 2.843602  | 3.317536 |
    | Non-specific     | 5.21327  | 8.056872  | 11.848341 | 4.739336 |
    | Specific         | 7.109005 | 10.900474 | 8.530806  | 9.952607 |
    | Spurious         | 6.161137 | 4.265403  | 1.895735  | 7.109005 |
    the effect of s1.tsnr.l on s1.cat.aCompCor3 in rat is q =  23.88 with p-value = 0.00449, dof = 9
    looking at main.experimenter.gender effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | f         | m        |
    |------------------|-----------|----------|
    | No               | 3.97351   | 6.622517 |
    | Non-specific     | 31.788079 | 3.311258 |
    | Specific         | 29.801325 | 5.298013 |
    | Spurious         | 13.907285 | 5.298013 |
    the effect of main.experimenter.gender on s1.cat.aCompCor3 in rat is q =  23.12 with p-value = 4e-05, dof = 3
    looking at rodent.sex effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | f        | m         |
    |------------------|----------|-----------|
    | No               | 1.895735 | 12.322275 |
    | Non-specific     | 0.947867 | 28.909953 |
    | Specific         | 4.265403 | 32.227488 |
    | Spurious         | 1.895735 | 17.535545 |
    the effect of rodent.sex on s1.cat.aCompCor3 in rat is q =  4.0 with p-value = 0.26099, dof = 3
    looking at head-plate effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | n         | y         |
    |------------------|-----------|-----------|
    | No               | 9.004739  | 5.21327   |
    | Non-specific     | 10.900474 | 18.957346 |
    | Specific         | 18.483412 | 18.009479 |
    | Spurious         | 6.635071  | 12.796209 |
    the effect of head-plate on s1.cat.aCompCor3 in rat is q =  8.85 with p-value = 0.03131, dof = 3
    looking at body.restrained effect in s1.cat.aCompCor3 in rat
    | s1.cat.aCompCor3 | y         |
    |------------------|-----------|
    | No               | 14.218009 |
    | Non-specific     | 29.85782  |
    | Specific         | 36.492891 |
    | Spurious         | 19.43128  |
    the effect of body.restrained on s1.cat.aCompCor3 in rat is q =  0.0 with p-value = 1.0, dof = 0
