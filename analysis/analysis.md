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
#import stats modules
from scipy.stats import chi2
from pingouin import ttest
import statsmodels.formula.api as smf

import matplotlib
#set matplotlib rendering to TkAgg
matplotlib.use('TkAgg')


def makeswarmplot(title, x, y, hue=None):
    ax = plt.figure(figsize=(10,5))
    ax = sns.stripplot(x=x, y=y, hue=hue)
    plt.title(title)
    return ax

def makeviolinplot(df_to_plot):
  plt.figure(figsize=(10,5))
  ax = sns.violinplot(x=df_to_plot['value'], y=df_to_plot['cont_variable'], hue=df_to_plot['value'],  inner=None)
  ax = sns.stripplot(x=df_to_plot['value'], y=df_to_plot['cont_variable'],  dodge=True, color='black', alpha=0.3)
  return ax

pl.Config(
    tbl_formatting="MARKDOWN",
    tbl_hide_column_data_types=True,
    tbl_hide_dataframe_shape=True,
    tbl_rows=100,
)

analysis_list = [ "gsr1", "gsr2", "gsr3","wmcsf1", "wmcsf2", "wmcsf3", "aCompCor1", "aCompCor2", "aCompCor3" ]  

roi_list=["s1","thal"]
```

## Whole analysis.

``` python
rodent_list = ["mouse"]
for rodent in rodent_list:
#rodent='mouse'
  print("#### NOW DOING " + rodent + " ####")
  df = pl.read_csv("../assets/tables/"+rodent+"_metadata_processed.tsv", separator="\t").sort(by='rodent.ds')
  df_summary=pl.read_csv("../assets/tables/"+rodent+"_summary_processed.tsv", separator="\t").sort(by='rodent.ds')

  print("summary of the data that we collected")
  print("we processed " + str(df_summary["rodent.ds"].count()) + " datasets") 
  print("totalling "+ str(df_summary["total_run"].sum()) +" runs")
  print("from "+str(df_summary["total_animal"].sum()) +" animals")
  print("the smallest dataset had "+ str(df_summary["total_run"].min())+" runs") 
  print("the largest dataset had "+str(df_summary["total_run"].max())+" runs")
  print("we could processed "+str(df_summary["total_included"].sum()) + "/" + str(df_summary["total_run"].sum())+ " runs.")

  print("below is a summary of the data included per dataset")
#to add the summary of the data included
  print(df_summary.select("rodent.ds", "total_run", "total_animal", "total_included", "strain"))

  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_summary["rodent.ds"], y=df_summary["total_included"], hue=df_summary["strain"], dodge=False)
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  ax.figure.savefig("../assets/plot/"+rodent+"_run_included_per_dataset.svg")

  print("information about sex ratio")
  print("the datasets contained "+str(df_summary["male"].sum())+ " male runs and " +str(df_summary["female"].sum())+ " female runs")
  print("that corresponds to "+str(round(100*df_summary["female"].sum()/(df_summary["female"].sum()+df_summary["male"].sum()),2))+"% females ")
  print("information about animal handling")
  print(str((df_summary["headplate"] == 'y').sum()) + "/" +str(len(df_summary)) +" datasets used headplates")
  print(str((df_summary["restrained"] == 'y').sum()) + "/" +str(len(df_summary)) +" datasets used body restraining")
  print(str((df_summary["anesthesia"] == 'y').sum()) + "/" +str(len(df_summary)) +" datasets used anesthesia before acquisition")
  print(str((df_summary["exp.gender"] == 'm').sum()) + "/" +str(len(df_summary)) +" datasets were collected by men, " + str((df_summary["exp.gender"] == 'f').sum())+ " by women")
  

  print(df_summary.select("rodent.ds", "headplate", "restrained", "anesthesia","exp.gender", "habituation.days","habituation.min"))

  fig, ax1 = plt.subplots()
  ax2 = ax1.twinx()
  ax = sns.barplot(x=df_summary["rodent.ds"], y=df_summary["habituation.min"], hue=df_summary["headplate"], dodge=False, ax=ax1)
  ax = sns.swarmplot(x=df_summary["rodent.ds"], y=df_summary["habituation.days"], color='black', marker='o', label='habituation days', legend=False, ax=ax2)
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  ax.figure.savefig("../assets/plot/"+rodent+"_habituation_per_dataset.svg")

  print("information about the scanner and sequence")
  print("lowest field strength was " + str(df_summary["field_strength"].min()) + " T")
  print("highest field strength was " + str(df_summary["field_strength"].max()) + " T")
  print(df_summary.select("rodent.ds", "field_strength", "sequence", "TE").sort(by="rodent.ds"))


#first, let's extract some infomation about motion and summarize it per dataset
  print("#### MOTION ANALYISIS ####")
  print("mean fd across all "+rodent+" datasets")
  print(df.select("fd.mean").mean())
  print("mean fd per dataset")
  print(df_summary.select("rodent.ds", "fd.mean").sort(by="rodent.ds"))

  plt.figure(figsize=(10,5))
  ax = makeswarmplot('framewise displacement per dataset', df["rodent.ds"], df["fd.mean"], hue=df['head-plate'])
  ax.figure.savefig("../assets/plot/"+rodent+"_fd_per_dataset.svg")

  df_to_plot = df.select('head-plate','fd.mean').rename({'head-plate':'value','fd.mean':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fd_headplate_violin.svg")

  #count how many rows of df have fd.max > 0.45 (corresponding roughly to 1.5 inplane voxel size).
  print("number of scans with max FD > 0.45 mm: " + str((df['fd.max']>0.45).sum()) + "/" + str(df.height))
  
  t = ttest(df.filter(pl.col('head-plate')=='y')['fd.mean'], df.filter(pl.col('head-plate')=='n')['fd.mean'])
  print('t-test for head-plate mean.fd > no head-plate mean.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')

  plt.figure(figsize=(10,5))
  ax = makeswarmplot('max framewise displacement per dataset', df["rodent.ds"], df["fd.max"], hue=df['head-plate'])
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_per_dataset.svg")

  df_to_plot = df.select('head-plate','fd.max').rename({'head-plate':'value','fd.max':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_headplate_violin.svg")
  
  t = ttest(df.filter(pl.col('head-plate')=='y')['fd.max'], df.filter(pl.col('head-plate')=='n')['fd.max'])
  print('t-test for head-plate max.fd > no head-plate max.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')


  print("sanity test for manuscript, habituation on motion")
  t = ttest(df.filter(pl.col('short.habituation')=='long')['fd.max'], df.filter(pl.col('short.habituation')=='short')['fd.max'])
  print('t-test for long habituation max.fd > short habituation max.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  
  df_to_plot = df.select('short.habituation','fd.max').rename({'short.habituation':'value','fd.max':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_habituation_violin.svg")
  

  print("sanity test for manuscript, sex on motion")
  t = ttest(df.filter(pl.col('rodent.sex')=='m')['fd.max'], df.filter(pl.col('rodent.sex')=='f')['fd.max'])
  print('t-test for male rodent max.fd > female rodent max.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  
  print("sanity test for manuscript, exp gender on motion")
  t = ttest(df.filter(pl.col('main.experimenter.gender')=='m')['fd.max'], df.filter(pl.col('main.experimenter.gender')=='f')['fd.max'])
  print('t-test for man max.fd > woman max.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')

  df_to_plot = df.select('main.experimenter.gender','fd.max').rename({'main.experimenter.gender':'value','fd.max':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_gender_violin.svg")

  print("sanity test for manuscript, exp gender on motion")
  t = ttest(df.filter(pl.col('main.experimenter.gender')=='m')['fd.mean'], df.filter(pl.col('main.experimenter.gender')=='f')['fd.mean'])
  print('t-test for man mean.fd > woman mean.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
 

#let's extract some infomation about tsnr and summarize it per dataset
  print("#### tSNR ANALYISIS ####")
  print("tSNR across all "+rodent+" datasets")
  print(df.select("s1.tsnr.l").mean())

  plt.figure(figsize=(10,5))
  ax = makeswarmplot('S1 tSNR per dataset', df["rodent.ds"], df["s1.tsnr.l"], hue=df['MRI.field.strength'].cast(pl.String))
  ax.figure.savefig("../assets/plot/"+rodent+"_tsnr_per_dataset.svg")

  df_to_plot = df.select('MRI.field.strength','s1.tsnr.l').rename({'MRI.field.strength':'value','s1.tsnr.l':'cont_variable'}).sort(by='value').with_columns(pl.col('value').cast(pl.String))
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_tsnr_field_violin.svg")
 

#now we run the analysis per denoising style, we extract the number of dropped frames, the s1-s1, s1-aca, and s1-thal correlations. finally we estimate connectivity specificity
  print("Number of dropped frames for each dataset and denoising method")
  print("This corresponds to the following rabies flags for mice")
  print("#gsr1: --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true")
  print("#gsr2: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true")
  print("#gsr3: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false")
  print("#### DENOISE ANALYSIS ####")
  print("dropped frames per dataset")
  print(df_summary.select("rodent.ds", "total.frames", "dropped.frames.gsr1", "dropped.frames.gsr2", "dropped.frames.gsr3").sort(by="rodent.ds"))

  #divide drpopped.frames.gsr1 by total.frames in df to get percentage and mean across all rows
  print("mean percentage of dropped frames across all datasets")
  print(df.select(
      (100*pl.col("dropped.frames.gsr1")/pl.col("total.frames")).alias("dropped.frames.gsr1.perc"),
      (100*pl.col("dropped.frames.gsr2")/pl.col("total.frames")).alias("dropped.frames.gsr2.perc"),
      (100*pl.col("dropped.frames.gsr3")/pl.col("total.frames")).alias("dropped.frames.gsr3.perc"),
  ).mean())

  print("mean percentage of dropped frames per dataset")
  print(df_summary.select("rodent.ds",
      (100*pl.col("dropped.frames.gsr1")/pl.col("total.frames")).alias("dropped.frames.gsr1.perc"),
      (100*pl.col("dropped.frames.gsr2")/pl.col("total.frames")).alias("dropped.frames.gsr2.perc"),
      (100*pl.col("dropped.frames.gsr3")/pl.col("total.frames")).alias("dropped.frames.gsr3.perc"),
  ).sort(by="rodent.ds"))


### adding the dice analysis ####
  if rodent=="mouse":
    analysis = "gsr3"
    dr_s1=2
    dr_aca=9
  else:
    analysis = "gsr3"
    dr_s1=5
    dr_aca=7

  #first test if processing change the dice outcome

  print("do dice statistical analysis for s1")
  cols = [c for c in df.columns if c.startswith(f"dr{dr_s1}")]
  df_dice = df.select(["scan", *cols]).unpivot(index=["scan"], on=cols).drop_nulls().to_pandas()
  
  # Fit linear mixed model and reduced model
  full_model= smf.mixedlm("value ~ variable", df_dice, groups="scan").fit()
  print(full_model.summary())

  reduced_model = smf.mixedlm("value ~ 1", data=df_dice, groups="scan").fit()
 
  # carry out likelyhood test ratio
  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_modelwc - reduced_model.df_modelwc

  p_value = chi2.sf(lr_stat, df_diff)

  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")
  
  plt.figure(figsize=(10,10))
  ax = sns.ecdfplot(data=df_dice, x="value", hue="variable")
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_per_processing_ecdf.svg")

  plt.figure(figsize=(10,10))
  ax = sns.violinplot(x="variable", y="value", data=df_dice, inner="box", hue="variable")
  ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_per_processing.svg")

  print("do dice stat analysis with aca compoment")
  cols = [c for c in df.columns if c.startswith(f"dr{dr_aca}")]
  df_dice = df.select(["scan", *cols]).unpivot(index=["scan"], on=cols).drop_nulls().to_pandas()
  
  # Fit linear mixed model and reduced model
  full_model= smf.mixedlm("value ~ variable", df_dice, groups="scan").fit()
  print(full_model.summary())

  reduced_model = smf.mixedlm("value ~ 1", data=df_dice, groups="scan").fit()
 
  # carry out likelyhood test ratio
  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_modelwc - reduced_model.df_modelwc

  p_value = chi2.sf(lr_stat, df_diff)

  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")

  plt.figure(figsize=(10,10))
  ax = sns.ecdfplot(data=df_dice, x="value", hue="variable")
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_per_processing_ecdf.svg")


  plt.figure(figsize=(10,10))
  ax = sns.violinplot(x="variable", y="value", data=df_dice, inner="box", hue="variable")
  ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_per_processing.svg")

 
  #create two new columns in df, one for the dice score for dr_s1 from the column f"dr{dr_s1}.{analysis} and one for the dice score for dr_aca
  df = df.with_columns(pl.col(f"dr{dr_s1}.{analysis}").alias('dice_s1'))
  df = df.with_columns(pl.col(f"dr{dr_aca}.{analysis}").alias('dice_aca'))

  #plot the dice for s1 and aca per dataset
  plt.figure(figsize=(10,10))
  ax = makeswarmplot('Dice ICA S1 per dataset', df["rodent.ds"], df["dice_s1"], hue=df['head-plate'])
  ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_per_dataset.svg")

  plt.figure(figsize=(10,10))
  ax = makeswarmplot('Dice ICA ACA per dataset', df["rodent.ds"], df["dice_aca"], hue=df['head-plate'])
  ax.set_xticklabels(ax.get_xticklabels(), rotation=45)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_per_dataset.svg")

  #run a couple ttest to recapitulate observations about specificity.
  t = ttest(df.filter(pl.col('head-plate')=='y')['dice_s1'], df.filter(pl.col('head-plate')=='n')['dice_s1'])
  print('t-test for head-plate dice_s1 > no head-plate dice_s1') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('head-plate','dice_s1').rename({'head-plate':'value','dice_s1':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_headplate_violin.svg")

  t = ttest(df.filter(pl.col('head-plate')=='y')['dice_aca'], df.filter(pl.col('head-plate')=='n')['dice_aca'])
  print('t-test for head-plate dice_aca > no head-plate dice_aca') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('head-plate','dice_aca').rename({'head-plate':'value','dice_aca':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_headplate_violin.svg")

  t = ttest(df.filter(pl.col('short.habituation')=='long')['dice_s1'], df.filter(pl.col('short.habituation')=='short')['dice_s1'])
  print('t-test for long habituation dice_s1 > short habituation dice_s1') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('short.habituation','dice_s1').rename({'short.habituation':'value','dice_s1':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_habituation_violin.svg")


  t = ttest(df.filter(pl.col('short.habituation')=='long')['dice_aca'], df.filter(pl.col('short.habituation')=='short')['dice_aca'])
  print('t-test for long habituation dice_aca > short habituation dice_aca') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('short.habituation','dice_aca').rename({'short.habituation':'value','dice_aca':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_habituation_violin.svg")


  t = ttest(df.filter(pl.col('fMRI.sequence')=='GE-EPI')['dice_s1'], df.filter(pl.col('fMRI.sequence')=='SE-EPI')['dice_s1'])
  print('t-test for GE-EPI dice_s1 > SE-EPI dice_s1') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('fMRI.sequence','dice_s1').rename({'fMRI.sequence':'value','dice_s1':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_sequence_violin.svg")


  t = ttest(df.filter(pl.col('fMRI.sequence')=='GE-EPI')['dice_aca'], df.filter(pl.col('fMRI.sequence')=='SE-EPI')['dice_aca'])
  print('t-test for GE-EPI dice_aca > SE-EPI dice_aca') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('fMRI.sequence','dice_aca').rename({'fMRI.sequence':'value','dice_aca':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_sequence_violin.svg")


  t = ttest(df.filter(pl.col('main.experimenter.gender')=='m')['dice_s1'], df.filter(pl.col('main.experimenter.gender')=='f')['dice_s1'])
  print('t-test for man dice_s1 > woman dice_s1') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('main.experimenter.gender','dice_s1').rename({'main.experimenter.gender':'value','dice_s1':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_gender_violin.svg")

  t = ttest(df.filter(pl.col('anesthesia.before.acquisition')=='y')['dice_s1'], df.filter(pl.col('anesthesia.before.acquisition')=='n')['dice_s1'])
  print('t-test for anesthesia dice_s1 > no anesthesia dice_s1') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p_val"].item(),5)}, dof={round(t["dof"].item(),2)}')
  df_to_plot = df.select('anesthesia.before.acquisition','dice_s1').rename({'anesthesia.before.acquisition':'value','dice_s1':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_anesthesia_violin.svg")

 
  print("analysis of dice in s1 as a function of mean framewise displacement")
  df = df.with_columns(pl.col("fd.mean").alias("fd_mean"), pl.col("fd.max").alias("fd_max"))

  full_model = smf.ols("dice_s1 ~ fd_mean", data=df).fit()
  reduced_model = smf.ols("dice_s1 ~ 1", data=df).fit()

  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_model - reduced_model.df_model
  p_value = chi2.sf(lr_stat, df_diff)
  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")

  plt.figure(figsize=(10,10))
  ax = sns.scatterplot(x=df["fd_mean"], y=df["dice_s1"])
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_vs_fd_mean.svg")

  full_model = smf.ols("dice_s1 ~ fd_max", data=df).fit()
  reduced_model = smf.ols("dice_s1 ~ 1", data=df).fit()

  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_model - reduced_model.df_model
  p_value = chi2.sf(lr_stat, df_diff)
  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")

  plt.figure(figsize=(10,10))
  ax = sns.scatterplot(x=df["fd_max"], y=df["dice_s1"])
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_s1_vs_fd_max.svg")


  print("analysis of dice in aca as a function of mean framewise displacement")
  full_model = smf.ols("dice_aca ~ fd_mean", data=df).fit()
  reduced_model = smf.ols("dice_aca ~ 1", data=df).fit()

  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_model - reduced_model.df_model
  p_value = chi2.sf(lr_stat, df_diff)
  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")

  plt.figure(figsize=(10,10))
  ax = sns.scatterplot(x=df["fd_mean"], y=df["dice_aca"])
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_vs_fd_mean.svg")

  full_model = smf.ols("dice_aca ~ fd_max", data=df).fit()
  reduced_model = smf.ols("dice_aca ~ 1", data=df).fit()

  lr_stat = 2 * (full_model.llf - reduced_model.llf)
  df_diff = full_model.df_model - reduced_model.df_model
  p_value = chi2.sf(lr_stat, df_diff)
  print(f"LR statistic: {lr_stat}, df: {df_diff}, p-value: {p_value}")

  plt.figure(figsize=(10,10))
  ax = sns.scatterplot(x=df["fd_max"], y=df["dice_aca"])
  ax.figure.savefig("../assets/plot/"+rodent+"_dice_aca_vs_fd_max.svg")

  #finally, plot the dictionary components
  from nilearn.image import iter_img
  from nilearn.plotting import plot_stat_map, show
  from nibabel import load as load_nii

  bg_img = load_nii("../assets/template/"+rodent+"/template.nii.gz")
  dict_img = load_nii("../assets/template/"+rodent+"/dict_"+analysis+".nii.gz")
  for i, cur_img in enumerate(iter_img(dict_img)):
    plot_stat_map(
        cur_img,
        bg_img=bg_img,
        title=f"IC {int(i)}",
        threshold=0.03,
        vmax=0.1,
        vmin=-0.1,
        colorbar=False,
        output_file=f"../assets/plot/dict/{rodent}_dict_{analysis}_IC{int(i)}.svg"
    )
```

    #### NOW DOING mouse ####
    summary of the data that we collected
    we processed 20 datasets
    totalling 1566 runs
    from 198 animals
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 1484/1566 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain      |
    |-----------|-----------|--------------|----------------|-------------|
    | 1001      | 4         | 4            | 4              | C57BL/6     |
    | 1002      | 8         | 8            | 8              | C57BL/6     |
    | 1003      | 12        | 12           | 12             | C57BL/6     |
    | 1004      | 21        | 5            | 21             | C57BL/6     |
    | 1005      | 207       | 13           | 190            | C57BL/6     |
    | 1006      | 20        | 10           | 19             | C57BL/6     |
    | 1007      | 7         | 7            | 7              | C57BL/6     |
    | 1008      | 36        | 9            | 36             | C57BL/6     |
    | 1009      | 107       | 5            | 49             | 129S2/SvPas |
    | 1011      | 51        | 17           | 51             | C57BL/6     |
    | 1012      | 26        | 26           | 24             | C57BL/6     |
    | 1013      | 32        | 8            | 32             | C57BL/6     |
    | 1014      | 4         | 4            | 4              | ICR         |
    | 1015      | 6         | 6            | 3              | C57BL/6     |
    | 1016      | 28        | 7            | 28             | C57BL/6     |
    | 3001      | 112       | 6            | 112            | C57BL/6     |
    | 3002      | 10        | 10           | 10             | C57BL/6     |
    | 3003      | 479       | 19           | 478            | F1 C6/129P  |
    | 3004      | 54        | 9            | 54             | C57BL/6     |
    | 3005      | 342       | 13           | 342            | C57BL/6     |

    information about sex ratio
    the datasets contained 1303 male runs and 151 female runs
    that corresponds to 10.39% females 
    information about animal handling
    17/20 datasets used headplates
    13/20 datasets used body restraining
    11/20 datasets used anesthesia before acquisition
    12/20 datasets were collected by men, 8 by women
    | rodent.ds | headplate | restrained | anesthesia | exp.gender | habituation.day | habituation.min |
    |           |           |            |            |            | s               |                 |
    |-----------|-----------|------------|------------|------------|-----------------|-----------------|
    | 1001      | y         | n          | n          | m          | 7               | 315             |
    | 1002      | y         | n          | n          | m          | 9               | 500             |
    | 1003      | y         | y          | n          | m          | 19              | 695             |
    | 1004      | y         | y          | n          | f          | 10              | 233             |
    | 1005      | y         | y          | y          | f          | 7               | 1320            |
    | 1006      | y         | n          | y          | m          | 7               | 260             |
    | 1007      | y         | y          | n          | m          | 9               | 515             |
    | 1008      | n         | y          | n          | m          | 8               | 420             |
    | 1009      | y         | y          | y          | m          | 10              | 50              |
    | 1011      | y         | y          | y          | m          | 0               | 0               |
    | 1012      | n         | y          | y          | f          | 5               | 150             |
    | 1013      | y         | n          | n          | f          | 5               | 100             |
    | 1014      | y         | y          | y          | f          | 7               | 180             |
    | 1015      | y         | n          | y          | m          | 4               | 300             |
    | 1016      | y         | n          | n          | m          | 5               | 100             |
    | 3001      | y         | y          | y          | m          | 4               | 42              |
    | 3002      | y         | y          | n          | f          | 16              | 432             |
    | 3003      | y         | y          | y          | m          | 4               | 42              |
    | 3004      | n         | y          | y          | f          | 5               | 165             |
    | 3005      | y         | n          | y          | f          | 4               | 42              |

    information about the scanner and sequence
    lowest field strength was 7.0 T
    highest field strength was 15.2 T
    | rodent.ds | field_strength | sequence | TE     |
    |-----------|----------------|----------|--------|
    | 1001      | 9.4            | GE-EPI   | 15.0   |
    | 1002      | 7.0            | GE-EPI   | 16.385 |
    | 1003      | 7.0            | GE-EPI   | 15.0   |
    | 1004      | 11.7           | GE-EPI   | 11.0   |
    | 1005      | 9.4            | GE-EPI   | 12.5   |
    | 1006      | 11.7           | GE-EPI   | 0.014  |
    | 1007      | 9.4            | GE-EPI   | 12.0   |
    | 1008      | 9.4            | GE-EPI   | 12.0   |
    | 1009      | 9.4            | GE-EPI   | 15.0   |
    | 1011      | 15.2           | GE-EPI   | 11.7   |
    | 1012      | 9.4            | GE-EPI   | 6.0    |
    | 1013      | 9.4            | SE-EPI   | 14.0   |
    | 1014      | 7.0            | GE-EPI   | 17.0   |
    | 1015      | 11.7           | GE-EPI   | 14.0   |
    | 1016      | 9.4            | SE-EPI   | 13.0   |
    | 3001      | 9.4            | SE-EPI   | 18.398 |
    | 3002      | 7.0            | GE-EPI   | 15.0   |
    | 3003      | 9.4            | SE-EPI   | 18.398 |
    | 3004      | 11.7           | GE-EPI   | 15.0   |
    | 3005      | 9.4            | SE-EPI   | null   |
    #### MOTION ANALYISIS ####
    mean fd across all mouse datasets
    | fd.mean  |
    |----------|
    | 0.037097 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 1001      | 0.007488 |
    | 1002      | 0.026981 |
    | 1003      | 0.009839 |
    | 1004      | 0.062667 |
    | 1005      | 0.010512 |
    | 1006      | 0.064551 |
    | 1007      | 0.011723 |
    | 1008      | 0.028795 |
    | 1009      | 0.027444 |
    | 1011      | 0.009327 |
    | 1012      | 0.021375 |
    | 1013      | 0.022732 |
    | 1014      | 0.008774 |
    | 1015      | 0.037024 |
    | 1016      | 0.016518 |
    | 3001      | 0.041209 |
    | 3002      | 0.017535 |
    | 3003      | 0.045527 |
    | 3004      | 0.03758  |
    | 3005      | 0.049052 |
    number of scans with max FD > 0.45 mm: 561/1484
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=3.29, p=0.00127, dof=133.66
    t-test for head-plate max.fd > no head-plate max.fd
    t=3.58, p=0.00047, dof=137.55
    sanity test for manuscript, habituation on motion
    t-test for long habituation max.fd > short habituation max.fd
    t=-18.62, p=0.0, dof=915.06
    sanity test for manuscript, sex on motion
    t-test for male rodent max.fd > female rodent max.fd
    t=0.98, p=0.32981, dof=150.97
    sanity test for manuscript, exp gender on motion
    t-test for man max.fd > woman max.fd
    t=4.34, p=2e-05, dof=1341.64
    sanity test for manuscript, exp gender on motion
    t-test for man mean.fd > woman mean.fd
    t=4.01, p=6e-05, dof=1237.51
    #### tSNR ANALYISIS ####
    tSNR across all mouse datasets
    | s1.tsnr.l |
    |-----------|
    | 15.729699 |
    Number of dropped frames for each dataset and denoising method
    This corresponds to the following rabies flags for mice
    #gsr1: --frame_censoring FD_censoring=true,FD_threshold=0.1,DVARS_censoring=true
    #gsr2: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=true
    #gsr3: --frame_censoring FD_censoring=true,FD_threshold=0.5,DVARS_censoring=false
    #### DENOISE ANALYSIS ####
    dropped frames per dataset
    | rodent.ds | total.frames | dropped.frames.gsr1 | dropped.frames.gsr2 | dropped.frames.gsr3 |
    |-----------|--------------|---------------------|---------------------|---------------------|
    | 1001      | 600.0        | 27.75               | 27.75               | 0.0                 |
    | 1002      | 900.0        | 119.25              | 90.375              | 0.5                 |
    | 1003      | 1800.0       | 155.363636          | 155.363636          | 0.0                 |
    | 1004      | 334.0        | 134.52381           | 48.238095           | 3.809524            |
    | 1005      | 989.473684   | 103.547368          | 99.894737           | 0.0                 |
    | 1006      | 612.631579   | 228.842105          | 65.473684           | 16.526316           |
    | 1007      | 600.0        | 39.428571           | 39.428571           | 0.0                 |
    | 1008      | 400.0        | 81.777778           | 58.027778           | 0.0                 |
    | 1009      | 453.061224   | 76.734694           | 63.44898            | 1.244898            |
    | 1011      | 360.0        | 45.568627           | 45.568627           | 0.0                 |
    | 1012      | 300.0        | 28.958333           | 26.666667           | 0.583333            |
    | 1013      | 480.0        | 51.1875             | 48.8125             | 0.0                 |
    | 1014      | 150.0        | 14.0                | 14.0                | 0.0                 |
    | 1015      | 400.0        | 119.333333          | 60.333333           | 0.0                 |
    | 1016      | 400.0        | 80.464286           | 69.357143           | 0.0                 |
    | 3001      | 200.0        | 42.0                | 27.892857           | 2.25                |
    | 3002      | 1920.0       | 398.9               | 388.8               | 0.4                 |
    | 3003      | 200.0        | 59.301053           | 41.242105           | 2.774737            |
    | 3004      | 180.0        | 48.185185           | 32.722222           | 3.240741            |
    | 3005      | 200.0        | 56.511696           | 38.406433           | 2.71345             |
    mean percentage of dropped frames across all datasets
    | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |--------------------------|--------------------------|--------------------------|
    | 23.706435                | 16.713517                | 0.969816                 |
    mean percentage of dropped frames per dataset
    | rodent.ds | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |-----------|--------------------------|--------------------------|--------------------------|
    | 1001      | 4.625                    | 4.625                    | 0.0                      |
    | 1002      | 13.25                    | 10.041667                | 0.055556                 |
    | 1003      | 8.631313                 | 8.631313                 | 0.0                      |
    | 1004      | 40.27659                 | 14.442543                | 1.140576                 |
    | 1005      | 10.464894                | 10.095745                | 0.0                      |
    | 1006      | 37.353952                | 10.687285                | 2.697595                 |
    | 1007      | 6.571429                 | 6.571429                 | 0.0                      |
    | 1008      | 20.444444                | 14.506944                | 0.0                      |
    | 1009      | 16.936937                | 14.004505                | 0.274775                 |
    | 1011      | 12.657952                | 12.657952                | 0.0                      |
    | 1012      | 9.652778                 | 8.888889                 | 0.194444                 |
    | 1013      | 10.664062                | 10.169271                | 0.0                      |
    | 1014      | 9.333333                 | 9.333333                 | 0.0                      |
    | 1015      | 29.833333                | 15.083333                | 0.0                      |
    | 1016      | 20.116071                | 17.339286                | 0.0                      |
    | 3001      | 21.0                     | 13.946429                | 1.125                    |
    | 3002      | 20.776042                | 20.25                    | 0.020833                 |
    | 3003      | 29.650526                | 20.621053                | 1.387368                 |
    | 3004      | 26.769547                | 18.179012                | 1.800412                 |
    | 3005      | 28.255848                | 19.203216                | 1.356725                 |
    do dice statistical analysis for s1

                   Mixed Linear Model Regression Results
    ====================================================================
    Model:                MixedLM     Dependent Variable:     value     
    No. Observations:     13318       Method:                 REML      
    No. Groups:           1480        Scale:                  0.0002    
    Min. group size:      8           Log-Likelihood:         35464.3850
    Max. group size:      9           Converged:              Yes       
    Mean group size:      9.0                                           
    --------------------------------------------------------------------
                              Coef. Std.Err.    z    P>|z| [0.025 0.975]
    --------------------------------------------------------------------
    Intercept                 0.215    0.002 121.297 0.000  0.211  0.218
    variable[T.dr2.aCompCor2] 0.005    0.000  10.085 0.000  0.004  0.005
    variable[T.dr2.aCompCor3] 0.015    0.000  32.481 0.000  0.014  0.016
    variable[T.dr2.gsr1]      0.016    0.000  35.648 0.000  0.015  0.017
    variable[T.dr2.gsr2]      0.021    0.000  47.340 0.000  0.021  0.022
    variable[T.dr2.gsr3]      0.027    0.000  59.622 0.000  0.026  0.028
    variable[T.dr2.wmcsf1]    0.013    0.000  28.038 0.000  0.012  0.014
    variable[T.dr2.wmcsf2]    0.018    0.000  39.083 0.000  0.017  0.019
    variable[T.dr2.wmcsf3]    0.023    0.000  51.236 0.000  0.022  0.024
    scan Var                  0.004    0.014                            
    ====================================================================

    LR statistic: 4645.298795893643, df: 8, p-value: 0.0

    do dice stat analysis with aca compoment

                   Mixed Linear Model Regression Results
    ====================================================================
    Model:                MixedLM     Dependent Variable:     value     
    No. Observations:     13318       Method:                 REML      
    No. Groups:           1480        Scale:                  0.0001    
    Min. group size:      8           Log-Likelihood:         37518.6128
    Max. group size:      9           Converged:              Yes       
    Mean group size:      9.0                                           
    --------------------------------------------------------------------
                              Coef. Std.Err.    z    P>|z| [0.025 0.975]
    --------------------------------------------------------------------
    Intercept                 0.200    0.001 178.730 0.000  0.198  0.202
    variable[T.dr9.aCompCor2] 0.004    0.000  10.968 0.000  0.004  0.005
    variable[T.dr9.aCompCor3] 0.013    0.000  31.339 0.000  0.012  0.013
    variable[T.dr9.gsr1]      0.018    0.000  43.713 0.000  0.017  0.018
    variable[T.dr9.gsr2]      0.022    0.000  55.124 0.000  0.021  0.023
    variable[T.dr9.gsr3]      0.028    0.000  69.621 0.000  0.027  0.029
    variable[T.dr9.wmcsf1]    0.014    0.000  34.923 0.000  0.013  0.015
    variable[T.dr9.wmcsf2]    0.018    0.000  45.648 0.000  0.018  0.019
    variable[T.dr9.wmcsf3]    0.024    0.000  58.949 0.000  0.023  0.025
    scan Var                  0.002    0.006                            
    ====================================================================

    LR statistic: 6054.304441552813, df: 8, p-value: 0.0

    t-test for head-plate dice_s1 > no head-plate dice_s1
    t=-7.74, p=0.0, dof=142.75
    t-test for head-plate dice_aca > no head-plate dice_aca
    t=-2.13, p=0.03499, dof=137.86
    t-test for long habituation dice_s1 > short habituation dice_s1
    t=40.55, p=0.0, dof=633.37
    t-test for long habituation dice_aca > short habituation dice_aca
    t=27.5, p=0.0, dof=719.8
    t-test for GE-EPI dice_s1 > SE-EPI dice_s1
    t=32.42, p=0.0, dof=688.04
    t-test for GE-EPI dice_aca > SE-EPI dice_aca
    t=24.94, p=0.0, dof=709.32
    t-test for man dice_s1 > woman dice_s1
    t=-10.63, p=0.0, dof=1378.73
    t-test for anesthesia dice_s1 > no anesthesia dice_s1
    t=-17.13, p=0.0, dof=194.25
    analysis of dice in s1 as a function of mean framewise displacement
    LR statistic: 585.2830137294677, df: 1.0, p-value: 2.6602754866042515e-129
    LR statistic: 400.49186466612446, df: 1.0, p-value: 4.303898898509776e-89
    analysis of dice in aca as a function of mean framewise displacement
    LR statistic: 368.4429273994301, df: 1.0, p-value: 4.085253211677793e-82
    LR statistic: 259.66042187658786, df: 1.0, p-value: 2.0348680132684025e-58
