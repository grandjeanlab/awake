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
from scipy.stats import chi2_contingency, ttest_ind
from pingouin import ttest

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

def makestackplot(df_melted):
  plt.figure(figsize=(10,5))
  ax = sns.displot(x=df_melted['variable'], hue=df_melted['value'], multiple='stack', hue_order=cat_index, stat='probability')
  ax.ax.invert_yaxis()
  ax.ax.tick_params(axis='x', rotation=90)
  return ax

def makespecificityplot(title, x,y, hue=None):
  ax = sns.jointplot(x=x, y=y, hue=hue)
  ax.fig.suptitle(title)
  ax.fig.subplots_adjust(top=0.9)
  ax.ax_joint.set(xlabel='Specific ROI [r]', ylabel='Unspecific ROI [r]')
  ax.ax_joint.set_xlim(-1, 1)
  ax.ax_joint.set_ylim(-1, 1)
  ax.ax_joint.vlines(0.1,ymin=-1,ymax=1,linestyles='dashed', color='black')
  ax.ax_joint.vlines(-0.1, -0.1,0.1,linestyles='dashed', color='black')
  ax.ax_joint.hlines(-0.1, -0.1,0.1,linestyles='dashed', color='black')
  ax.ax_joint.hlines(0.1, -0.1,xmax=1,linestyles='dashed', color='black')
  ax.ax_marg_x.axvline(x=0.1, color='black')
  ax.ax_marg_y.axhline(y=0.1, color='black')
  return ax

def split_continuous(df, cat1, cat2, split=4, labels = ["lowest","low","high","highest"]):
  split_table = df.select(cat1, cat2).with_columns(pl.col(cat2).qcut(split, labels=labels, allow_duplicates=True).alias('quartiles')).drop_nulls()
  return split_table

def chi2_continuous(df, cat1, cat2, rodent, split=4, labels = ["lowest","low","high","highest"]):
  chi2_table = split_continuous(df, cat1, cat2, split, labels)
  print(f'looking at {cat2} effect in {cat1} in {rodent}')
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

## Whole analysis.

``` python
rodent_list = ["mouse", "rat"]
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
  print(f't={round(t["T"].item(),2)}, p={round(t["p-val"].item(),5)}, dof={round(t["dof"].item(),2)}')

  plt.figure(figsize=(10,5))
  ax = makeswarmplot('max framewise displacement per dataset', df["rodent.ds"], df["fd.max"], hue=df['head-plate'])
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_per_dataset.svg")

  df_to_plot = df.select('head-plate','fd.max').rename({'head-plate':'value','fd.max':'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fdmax_headplate_violin.svg")
  
  t = ttest(df.filter(pl.col('head-plate')=='y')['fd.max'], df.filter(pl.col('head-plate')=='n')['fd.max'])
  print('t-test for head-plate max.fd > no head-plate mean.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p-val"].item(),5)}, dof={round(t["dof"].item(),2)}')



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

  print("#### plot fc categories per denoising method ####")
  df_melted = df.melt(id_vars="rodent.ds", value_vars=['s1.cat.' + x for x in analysis_list])
  plt.figure(figsize=(10,5))
  ax = makestackplot(df_melted)
  ax.figure.savefig("../assets/plot/"+rodent+"_specificity.svg")

  chi2_categorical(df_melted, 'variable','value', 'mouse')

  if rodent=="mouse":
    analysis = "wmcsf3"
  else:
    analysis = "gsr3"

  plt.figure(figsize=(10,5))
  ax = makespecificityplot(rodent+': S1 specificity with ' + analysis, df['s1.specific.'+analysis], df['s1.unspecific.'+analysis])
  ax.figure.savefig("../assets/plot/specific/"+rodent+"_s1_specificity_"+analysis+".svg")

  plt.figure(figsize=(10,5))
  ax = makespecificityplot(rodent+': Thal FC specificity with ' + analysis, df['thal.specific.'+analysis], df['s1.unspecific.'+analysis])
  ax.figure.savefig("../assets/plot/specific/"+rodent+"_thal_specificity_"+analysis+".svg")

  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))

  cat2 = 'fd.mean'
  df_to_plot = split_continuous(df, cat1, cat2)
  df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})

  plt.figure(figsize=(10,5))
  ax = makestackplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_fd.svg")

  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_fd_violin.svg")

  cat2 = 'fd.max'
  df_to_plot = split_continuous(df, cat1, cat2)
  df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})

  plt.figure(figsize=(10,5))
  ax = makestackplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_fdmax.svg")

  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_fdmax_violin.svg")

  cat2='s1.gsrcov.l.'+analysis
  df_to_plot = split_continuous(df, cat1, cat2)
  df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makestackplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_gsrcov.svg")

  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_gsrcov_violin.svg")

  cat2='s1.tsnr.l'
  df_to_plot = split_continuous(df, cat1, cat2)
  df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makestackplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_tsnr.svg")

  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_tsnr_violin.svg")


  cat2='main.experimenter.gender'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_expgender.svg')

  cat2='rodent.sex'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_sex.svg')

  cat2='head-plate'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_headplate.svg')

  cat2='body.restrained'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_restrained.svg')

  cat2='short.habituation'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_shorthabituation.svg')

  cat2='anesthesia.before.acquisition'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_anesthesia.svg')

  cat2='MRI.field.strength'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_MRI.svg')

  cat2='fMRI.sequence'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(10,5))
  ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
  ax.figure.savefig('../assets/plot/'+rodent+'_'+analysis+'_fMRI.svg')

  cat2='habituation.min'
  df_to_plot = split_continuous(df, 'cat', cat2)
  df_to_plot = df_to_plot.rename({'cat':'value', 'quartiles':'variable', cat2:'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_habituation_min.svg")

  cat2='habituation.days'
  df_to_plot = split_continuous(df, 'cat', cat2)
  df_to_plot = df_to_plot.rename({'cat':'value', 'quartiles':'variable', cat2:'cont_variable'})
  plt.figure(figsize=(10,5))
  ax = makeviolinplot(df_to_plot)
  ax.set_xlabel('Connectivity category')
  ax.figure.savefig("../assets/plot/"+rodent+"_"+analysis+"_habituation_days.svg")



  print('doing statistical analysis')
  print('first with '+analysis+' processed data')
  cat2 = 'fd.mean'
  chi2_continuous(df, 'cat', cat2, rodent)
  cat2 = 'fd.max'
  chi2_continuous(df, 'cat', cat2, rodent)
  cat2 = 's1.gsrcov.l.'+analysis
  chi2_continuous(df, 'cat', cat2, rodent)
  cat2 = 's1.tsnr.l'
  chi2_continuous(df, 'cat', cat2, rodent)
  cat2 = 'habituation.min'
  chi2_continuous(df, 'cat', cat2, rodent, 2, ['low','high'])
  cat2 = 'habituation.days'
  chi2_continuous(df, 'cat', cat2, rodent, 2, ['low','high'])

  chi2_categorical(df, 'cat', 'short.habituation', rodent)
  chi2_categorical(df, 'cat', 'main.experimenter.gender', rodent)
  chi2_categorical(df, 'cat', 'rodent.sex', rodent)
  chi2_categorical(df, 'cat', 'head-plate', rodent)
  chi2_categorical(df, 'cat', 'body.restrained', rodent)
  chi2_categorical(df, 'cat', 'anesthesia.before.acquisition', rodent)
  chi2_categorical(df, 'cat', 'MRI.field.strength', rodent)
  chi2_categorical(df, 'cat', 'fMRI.sequence', rodent)
```

    #### NOW DOING mouse ####
    summary of the data that we collected
    we processed 20 datasets
    totalling 1579 runs
    from 198 animals
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 1361/1579 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain      |
    |-----------|-----------|--------------|----------------|-------------|
    | 1001      | 4         | 4            | 4              | C57BL/6     |
    | 1002      | 8         | 8            | 8              | C57BL/6     |
    | 1003      | 12        | 12           | 12             | C57BL/6     |
    | 1004      | 21        | 5            | 21             | C57BL/6     |
    | 1005      | 207       | 13           | 176            | C57BL/6     |
    | 1006      | 20        | 10           | 20             | C57BL/6     |
    | 1007      | 7         | 7            | 7              | C57BL/6     |
    | 1008      | 36        | 9            | 34             | C57BL/6     |
    | 1009      | 107       | 5            | 34             | 129S2/SvPas |
    | 1011      | 51        | 17           | 48             | C57BL/6     |
    | 1012      | 26        | 26           | 21             | C57BL/6     |
    | 1013      | 32        | 8            | 32             | C57BL/6     |
    | 1014      | 4         | 4            | 4              | ICR         |
    | 1015      | 6         | 6            | 6              | C57BL/6     |
    | 1016      | 28        | 7            | 28             | C57BL/6     |
    | 3001      | 112       | 6            | 112            | C57BL/6     |
    | 3002      | 10        | 10           | 10             | C57BL/6     |
    | 3003      | 479       | 19           | 414            | F1 C6/129P  |
    | 3004      | 54        | 9            | 54             | C57BL/6     |
    | 3005      | 355       | 13           | 316            | C57BL/6     |

    information about sex ratio
    the datasets contained 1316 male runs and 151 female runs
    that corresponds to 10.29% females 
    information about animal handling
    17/20 datasets used headplates
    13/20 datasets used body restraining
    11/20 datasets used anesthesia before acquisition
    12/20 datasets were collected by men, 8 by women
    | rodent.ds | headplate | restrained | anesthesia | exp.gender | habituation.day | habituation.min |
    |           |           |            |            |            | s               |                 |
    |-----------|-----------|------------|------------|------------|-----------------|-----------------|
    | 1001      | y         | n          | n          | m          | 7               | 315             |
    | 1002      | y         | n          | n          | m          | 0               | 0               |
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
    | 1005      | 9.4            | GE-EPI   | 11.0   |
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
    | 0.035054 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 1001      | 0.007517 |
    | 1002      | 0.032146 |
    | 1003      | 0.011058 |
    | 1004      | 0.065948 |
    | 1005      | 0.011795 |
    | 1006      | 0.06656  |
    | 1007      | 0.016836 |
    | 1008      | 0.033217 |
    | 1009      | 0.03189  |
    | 1011      | 0.011393 |
    | 1012      | 0.017715 |
    | 1013      | 0.025696 |
    | 1014      | 0.005115 |
    | 1015      | 0.034866 |
    | 1016      | 0.018205 |
    | 3001      | 0.035982 |
    | 3002      | 0.022321 |
    | 3003      | 0.040797 |
    | 3004      | 0.03515  |
    | 3005      | 0.046262 |
    number of scans with max FD > 0.45 mm: 537/1361
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=2.09, p=0.03895, dof=124.45
    t-test for head-plate max.fd > no head-plate mean.fd
    t=2.73, p=0.00729, dof=130.73
    #### tSNR ANALYISIS ####
    tSNR across all mouse datasets
    | s1.tsnr.l |
    |-----------|
    | 12.535792 |
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
    | 1005      | 977.272727   | 109.693182          | 105.227273          | null                |
    | 1006      | 618.0        | 216.55              | 77.5                | 18.866667           |
    | 1007      | 600.0        | 39.714286           | 39.714286           | null                |
    | 1008      | 400.0        | 83.823529           | 58.764706           | null                |
    | 1009      | 397.058824   | 69.323529           | 52.617647           | 5.111111            |
    | 1011      | 360.0        | 40.208333           | 40.208333           | null                |
    | 1012      | 300.0        | 26.714286           | 25.428571           | 4.333333            |
    | 1013      | 480.0        | 57.65625            | 52.96875            | null                |
    | 1014      | 150.0        | 13.0                | 12.75               | null                |
    | 1015      | 400.0        | 100.2               | 64.6                | null                |
    | 1016      | 400.0        | 86.0                | 75.714286           | null                |
    | 3001      | 200.0        | 38.142857           | 28.821429           | 4.404255            |
    | 3002      | 1920.0       | 400.2               | 388.0               | 3.0                 |
    | 3003      | 200.0        | 55.074879           | 41.342995           | 4.953608            |
    | 3004      | 180.0        | 51.296296           | 39.150943           | 7.75                |
    | 3005      | 200.0        | 52.234177           | 38.547468           | 4.965986            |
    mean percentage of dropped frames across all datasets
    | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |--------------------------|--------------------------|--------------------------|
    | 22.421984                | 16.953124                | 2.487458                 |
    mean percentage of dropped frames per dataset
    | rodent.ds | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |-----------|--------------------------|--------------------------|--------------------------|
    | 1001      | 7.041667                 | 7.041667                 | null                     |
    | 1002      | 12.5                     | 8.347222                 | null                     |
    | 1003      | 8.722222                 | 8.722222                 | null                     |
    | 1004      | 38.69404                 | 17.693185                | 1.766467                 |
    | 1005      | 11.224419                | 10.767442                | null                     |
    | 1006      | 35.040453                | 12.540453                | 3.052859                 |
    | 1007      | 6.619048                 | 6.619048                 | null                     |
    | 1008      | 20.955882                | 14.691176                | null                     |
    | 1009      | 17.459259                | 13.251852                | 1.287243                 |
    | 1011      | 11.168981                | 11.168981                | null                     |
    | 1012      | 8.904762                 | 8.47619                  | 1.444444                 |
    | 1013      | 12.011719                | 11.035156                | null                     |
    | 1014      | 8.666667                 | 8.5                      | null                     |
    | 1015      | 25.05                    | 16.15                    | null                     |
    | 1016      | 21.5                     | 18.928571                | null                     |
    | 3001      | 19.071429                | 14.410714                | 2.202128                 |
    | 3002      | 20.84375                 | 20.208333                | 0.15625                  |
    | 3003      | 27.53744                 | 20.671498                | 2.476804                 |
    | 3004      | 28.497942                | 21.750524                | 4.305556                 |
    | 3005      | 26.117089                | 19.273734                | 2.482993                 |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Spurious     | 0.418075 |
    | No           | 0.271124 |
    | Specific     | 0.218957 |
    | Non-specific | 0.091844 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Spurious     | 0.409258 |
    | No           | 0.28288  |
    | Specific     | 0.223365 |
    | Non-specific | 0.084497 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Spurious     | 0.379133 |
    | No           | 0.295371 |
    | Specific     | 0.240265 |
    | Non-specific | 0.085231 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.40191  |
    | No            | 0.270389 |
    | Specific      | 0.214548 |
    | Non-specific  | 0.113152 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.396032 |
    | No            | 0.283615 |
    | Specific      | 0.221161 |
    | Non-specific  | 0.099192 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Spurious      | 0.35856  |
    | No            | 0.287289 |
    | Specific      | 0.250551 |
    | Non-specific  | 0.1036   |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.415136 |
    | Specific         | 0.237325 |
    | No               | 0.22263  |
    | Non-specific     | 0.124908 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.415871 |
    | Specific         | 0.235121 |
    | No               | 0.224835 |
    | Non-specific     | 0.124173 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Spurious         | 0.369581 |
    | No               | 0.271859 |
    | Specific         | 0.23953  |
    | Non-specific     | 0.11903  |
    for gsr processing, heavy scrubbing lead to  298/1361 Specific scans
    for gsr processing, light scrubbing lead to  327/1361 Specific scans
    for aCompCor processing, heavy scrubbing lead to  323/1361 Specific scans
    for aCompCor processing, light scrubbing lead to  326/1361 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 1001      | 0.5       | 0.5       | 0.75      | … | 0.5       | 0.5       | 0.5       | 0.5      |
    | 1002      | 0.5       | 0.375     | 0.625     | … | 0.75      | 0.75      | 0.625     | 0.875    |
    | 1003      | 0.0       | 0.0       | 0.083333  | … | 0.083333  | 0.083333  | 0.083333  | 0.083333 |
    | 1004      | 0.047619  | 0.142857  | 0.238095  | … | 0.142857  | 0.047619  | 0.095238  | 0.095238 |
    | 1005      | 0.215909  | 0.210227  | 0.210227  | … | 0.25      | 0.25      | 0.221591  | 0.227273 |
    | 1006      | 0.2       | 0.3       | 0.15      | … | 0.4       | 0.3       | 0.4       | 0.45     |
    | 1007      | 0.142857  | 0.142857  | 0.0       | … | 0.0       | 0.142857  | 0.142857  | 0.142857 |
    | 1008      | 0.235294  | 0.235294  | 0.382353  | … | 0.323529  | 0.205882  | 0.235294  | 0.117647 |
    | 1009      | 0.235294  | 0.264706  | 0.235294  | … | 0.235294  | 0.176471  | 0.235294  | 0.352941 |
    | 1011      | 0.166667  | 0.166667  | 0.291667  | … | 0.208333  | 0.208333  | 0.270833  | 0.291667 |
    | 1012      | 0.333333  | 0.285714  | 0.333333  | … | 0.380952  | 0.428571  | 0.428571  | 0.380952 |
    | 1013      | 0.125     | 0.15625   | 0.15625   | … | 0.15625   | 0.125     | 0.21875   | 0.125    |
    | 1014      | 0.25      | 0.25      | 0.0       | … | 0.0       | 0.75      | 0.5       | 0.0      |
    | 1015      | 0.333333  | 0.166667  | 0.333333  | … | 0.5       | 0.333333  | 0.166667  | 0.5      |
    | 1016      | 0.321429  | 0.321429  | 0.25      | … | 0.392857  | 0.178571  | 0.25      | 0.25     |
    | 3001      | 0.303571  | 0.3125    | 0.267857  | … | 0.214286  | 0.303571  | 0.294643  | 0.223214 |
    | 3002      | 0.5       | 0.5       | 0.7       | … | 0.7       | 0.1       | 0.0       | 0.1      |
    | 3003      | 0.149758  | 0.154589  | 0.161836  | … | 0.178744  | 0.207729  | 0.198068  | 0.21256  |
    | 3004      | 0.592593  | 0.62963   | 0.611111  | … | 0.555556  | 0.407407  | 0.388889  | 0.37037  |
    | 3005      | 0.21519   | 0.212025  | 0.253165  | … | 0.272152  | 0.231013  | 0.224684  | 0.246835 |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 11/20 dataset performed better with less scrubbing
    for aCompCor processing, 9/20 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 8/20 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 8/20 dataset performed worst with aCompCor compared to gsr
    #### plot fc categories per denoising method ####

    looking at value effect in variable in mouse
    | variable         | No       | Non-specific | Specific | Spurious |
    |------------------|----------|--------------|----------|----------|
    | s1.cat.aCompCor1 | 2.473671 | 1.387868     | 2.63695  | 4.612621 |
    | s1.cat.aCompCor2 | 2.498163 | 1.379704     | 2.612458 | 4.620785 |
    | s1.cat.aCompCor3 | 3.020655 | 1.322557     | 2.661442 | 4.106458 |
    | s1.cat.gsr1      | 3.012491 | 1.020491     | 2.432852 | 4.645277 |
    | s1.cat.gsr2      | 3.143114 | 0.938852     | 2.481835 | 4.54731  |
    | s1.cat.gsr3      | 3.281901 | 0.947016     | 2.669606 | 4.212589 |
    | s1.cat.wmcsf1    | 3.004327 | 1.257245     | 2.383868 | 4.465671 |
    | s1.cat.wmcsf2    | 3.151278 | 1.102131     | 2.457343 | 4.400359 |
    | s1.cat.wmcsf3    | 3.192097 | 1.151114     | 2.783901 | 3.983999 |
    the effect of value on variable in mouse is q =  74.75 with p-value = 0.0, dof = 24

    doing statistical analysis
    first with wmcsf3 processed data
    looking at fd.mean effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 5.813098  | 6.254599  | 5.813098  | 7.211185  |
    | other    | 19.205298 | 18.763797 | 19.131714 | 17.807211 |
    the effect of fd.mean on cat in mouse is q =  3.72 with p-value = 0.29326, dof = 3
    looking at fd.max effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 6.033848  | 7.137601  | 6.033848  | 5.886681  |
    | other    | 18.984547 | 17.880795 | 18.910964 | 19.131714 |
    the effect of fd.max on cat in mouse is q =  2.9 with p-value = 0.40764, dof = 3
    looking at s1.gsrcov.l.wmcsf3 effect in cat in mouse
    | cat      | lowest    | low       | high     | highest   |
    |----------|-----------|-----------|----------|-----------|
    | Specific | 7.60709   | 6.20384   | 5.243722 | 6.05613   |
    | other    | 17.429838 | 18.759232 | 19.71935 | 18.980798 |
    the effect of s1.gsrcov.l.wmcsf3 on cat in mouse is q =  8.21 with p-value = 0.04195, dof = 3
    looking at s1.tsnr.l effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 7.579102  | 4.562178  | 4.562178  | 8.388521  |
    | other    | 17.439294 | 20.456218 | 20.382634 | 16.629875 |
    the effect of s1.tsnr.l on cat in mouse is q =  34.62 with p-value = 0.0, dof = 3
    looking at habituation.min effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 14.695077 | 10.360029 |
    | other    | 51.285819 | 23.659074 |
    the effect of habituation.min on cat in mouse is q =  10.46 with p-value = 0.00122, dof = 1
    looking at habituation.days effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 14.915503 | 10.139603 |
    | other    | 51.506245 | 23.438648 |
    the effect of habituation.days on cat in mouse is q =  9.28 with p-value = 0.00232, dof = 1
    looking at short.habituation effect in cat in mouse
    | cat      | long      | short     |
    |----------|-----------|-----------|
    | Specific | 10.139603 | 14.915503 |
    | other    | 23.438648 | 51.506245 |
    the effect of short.habituation on cat in mouse is q =  9.28 with p-value = 0.00232, dof = 1
    looking at main.experimenter.gender effect in cat in mouse
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 13.445996 | 11.609111 |
    | other    | 33.137399 | 41.807494 |
    the effect of main.experimenter.gender on cat in mouse is q =  8.8 with p-value = 0.00302, dof = 1
    looking at rodent.sex effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 3.923139 | 21.457166 |
    | other    | 5.684548 | 68.935148 |
    the effect of rodent.sex on cat in mouse is q =  15.85 with p-value = 7e-05, dof = 1
    looking at head-plate effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 3.600294 | 21.454813 |
    | other    | 4.408523 | 70.53637  |
    the effect of head-plate on cat in mouse is q =  23.85 with p-value = 0.0, dof = 1
    looking at body.restrained effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 8.890522  | 16.164585 |
    | other    | 21.528288 | 53.416605 |
    the effect of body.restrained on cat in mouse is q =  5.2 with p-value = 0.02258, dof = 1
    looking at anesthesia.before.acquisition effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 3.379868 | 21.675239 |
    | other    | 8.082292 | 66.862601 |
    the effect of anesthesia.before.acquisition on cat in mouse is q =  1.59 with p-value = 0.20786, dof = 1
    looking at MRI.field.strength effect in cat in mouse
    | cat      | 7.0      | 9.4       | 11.7     | 15.2     |
    |----------|----------|-----------|----------|----------|
    | Specific | 1.028655 | 20.05878  | 3.232917 | 0.734754 |
    | other    | 1.469508 | 66.495224 | 4.188097 | 2.792065 |
    the effect of MRI.field.strength on cat in mouse is q =  25.81 with p-value = 1e-05, dof = 3
    looking at fMRI.sequence effect in cat in mouse
    | cat      | GE-EPI    | SE-EPI    |
    |----------|-----------|-----------|
    | Specific | 10.286554 | 14.768553 |
    | other    | 23.218222 | 51.726672 |
    the effect of fMRI.sequence on cat in mouse is q =  11.2 with p-value = 0.00082, dof = 1
    #### NOW DOING rat ####
    summary of the data that we collected
    we processed 9 datasets
    totalling 467 runs
    from 161 animals
    the smallest dataset had 5 runs
    the largest dataset had 290 runs
    we could processed 237/467 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain              |
    |-----------|-----------|--------------|----------------|---------------------|
    | 2001      | 16        | 16           | 16             | Wistar              |
    | 2002a     | 7         | 7            | 7              | Sprague-Dawley      |
    | 2002b     | 5         | 5            | 5              | Wistar              |
    | 2002c     | 5         | 5            | 5              | Wistar              |
    | 2003      | 97        | 8            | 97             | Sprague-Dawley      |
    | 2004      | 19        | 5            | 19             | Long-Evans          |
    | 2005      | 5         | 3            | 5              | Wistar              |
    | 2006      | 23        | 23           | 23             | Black hooded Lister |
    | 4001      | 290       | 89           | 60             | Long-Evans          |

    information about sex ratio
    the datasets contained 445 male runs and 22 female runs
    that corresponds to 4.71% females 
    information about animal handling
    2/9 datasets used headplates
    9/9 datasets used body restraining
    8/9 datasets used anesthesia before acquisition
    6/9 datasets were collected by men, 3 by women
    | rodent.ds | headplate | restrained | anesthesia | exp.gender | habituation.day | habituation.min |
    |           |           |            |            |            | s               |                 |
    |-----------|-----------|------------|------------|------------|-----------------|-----------------|
    | 2001      | n         | y          | y          | f          | 8               | 110             |
    | 2002a     | n         | y          | y          | f          | 4               | 120             |
    | 2002b     | n         | y          | y          | m          | 4               | 120             |
    | 2002c     | n         | y          | y          | m          | 8               | 300             |
    | 2003      | y         | y          | y          | f          | 9               | 500             |
    | 2004      | y         | y          | n          | m          | 9               | 430             |
    | 2005      | n         | y          | y          | m          | 9               | 360             |
    | 2006      | n         | y          | y          | m          | 8               | 250             |
    | 4001      | n         | y          | y          | m          | 7               | 330             |
    information about the scanner and sequence
    lowest field strength was 7.0 T
    highest field strength was 9.4 T
    | rodent.ds | field_strength | sequence | TE |
    |-----------|----------------|----------|----|
    | 2001      | 7.0            | GE-EPI   | 17 |
    | 2002a     | 7.0            | GE-EPI   | 18 |
    | 2002b     | 7.0            | SE-EPI   | 45 |
    | 2002c     | 7.0            | SE-EPI   | 45 |
    | 2003      | 7.0            | GE-EPI   | 15 |
    | 2004      | 9.4            | GE-EPI   | 15 |
    | 2005      | 7.0            | GE-EPI   | 17 |
    | 2006      | 9.4            | GE-EPI   | 13 |
    | 4001      | 7.0            | GE-EPI   | 15 |
    #### MOTION ANALYISIS ####
    mean fd across all rat datasets
    | fd.mean  |
    |----------|
    | 0.041647 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 2001      | 0.03643  |
    | 2002a     | 0.052678 |
    | 2002b     | 0.044527 |
    | 2002c     | 0.044339 |
    | 2003      | 0.04796  |
    | 2004      | 0.016877 |
    | 2005      | 0.018197 |
    | 2006      | 0.045115 |
    | 4001      | 0.039548 |
    number of scans with max FD > 0.45 mm: 71/237
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=0.74, p=0.4587, dof=203.84
    t-test for head-plate max.fd > no head-plate mean.fd
    t=-6.46, p=0.0, dof=198.95
    #### tSNR ANALYISIS ####
    tSNR across all rat datasets
    | s1.tsnr.l |
    |-----------|
    | 28.147835 |
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
    | 2002b     | 750.0        | 144.4               | 93.0                | 6.0                 |
    | 2002c     | 750.0        | 141.8               | 101.0               | 6.0                 |
    | 2003      | 90.0         | 21.092784           | 9.185567            | 10.5                |
    | 2004      | 400.0        | 80.105263           | 79.157895           | null                |
    | 2005      | 1000.0       | 181.2               | 102.2               | 3.0                 |
    | 2006      | 750.0        | 137.173913          | 79.130435           | 11.3                |
    | 4001      | 600.0        | 103.85              | 79.716667           | 13.151515           |
    mean percentage of dropped frames across all datasets
    | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |--------------------------|--------------------------|--------------------------|
    | 19.976231                | 11.674121                | 3.384946                 |
    mean percentage of dropped frames per dataset
    | rodent.ds | dropped.frames.gsr1.perc | dropped.frames.gsr2.perc | dropped.frames.gsr3.perc |
    |-----------|--------------------------|--------------------------|--------------------------|
    | 2001      | 8.802083                 | 3.489583                 | 1.166667                 |
    | 2002a     | 28.447619                | 17.809524                | 0.377778                 |
    | 2002b     | 19.253333                | 12.4                     | 0.8                      |
    | 2002c     | 18.906667                | 13.466667                | 0.8                      |
    | 2003      | 23.436426                | 10.206186                | 11.666667                |
    | 2004      | 20.026316                | 19.789474                | null                     |
    | 2005      | 18.12                    | 10.22                    | 0.3                      |
    | 2006      | 18.289855                | 10.550725                | 1.506667                 |
    | 4001      | 17.308333                | 13.286111                | 2.191919                 |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Specific     | 0.400844 |
    | Spurious     | 0.227848 |
    | Non-specific | 0.21097  |
    | No           | 0.151899 |
    | null         | 0.008439 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Specific     | 0.400844 |
    | Spurious     | 0.227848 |
    | Non-specific | 0.21097  |
    | No           | 0.151899 |
    | null         | 0.008439 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Specific     | 0.417722 |
    | Spurious     | 0.240506 |
    | Non-specific | 0.202532 |
    | No           | 0.130802 |
    | null         | 0.008439 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Non-specific  | 0.329114 |
    | Specific      | 0.295359 |
    | Spurious      | 0.240506 |
    | No            | 0.130802 |
    | null          | 0.004219 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Specific      | 0.329114 |
    | Non-specific  | 0.2827   |
    | Spurious      | 0.248945 |
    | No            | 0.130802 |
    | null          | 0.008439 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Specific      | 0.35865  |
    | Non-specific  | 0.295359 |
    | Spurious      | 0.223629 |
    | No            | 0.113924 |
    | null          | 0.008439 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Specific         | 0.341772 |
    | Non-specific     | 0.28692  |
    | Spurious         | 0.219409 |
    | No               | 0.14346  |
    | null             | 0.008439 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Specific         | 0.350211 |
    | Non-specific     | 0.299578 |
    | Spurious         | 0.21097  |
    | No               | 0.135021 |
    | null             | 0.004219 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Specific         | 0.367089 |
    | Non-specific     | 0.28692  |
    | Spurious         | 0.189873 |
    | No               | 0.147679 |
    | null             | 0.008439 |
    for gsr processing, heavy scrubbing lead to  95/235 Specific scans
    for gsr processing, light scrubbing lead to  99/235 Specific scans
    for aCompCor processing, heavy scrubbing lead to  81/235 Specific scans
    for aCompCor processing, light scrubbing lead to  87/235 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 2001      | 0.625     | 0.625     | 0.625     | … | 0.4375    | 0.625     | 0.5       | 0.5625   |
    | 2002a     | 0.714286  | 0.571429  | 0.285714  | … | 0.142857  | 0.285714  | 0.285714  | 0.285714 |
    | 2002b     | 0.0       | 0.2       | 0.4       | … | 0.4       | 0.2       | 0.2       | 0.2      |
    | 2002c     | 0.6       | 0.6       | 0.2       | … | 0.2       | 0.4       | 0.8       | 0.6      |
    | 2003      | 0.402062  | 0.391753  | 0.463918  | … | 0.391753  | 0.329897  | 0.298969  | 0.350515 |
    | 2004      | 0.210526  | 0.210526  | 0.210526  | … | 0.263158  | 0.105263  | 0.157895  | 0.210526 |
    | 2005      | 0.0       | 0.0       | 0.0       | … | 0.0       | 0.0       | 0.0       | 0.2      |
    | 2006      | 0.565217  | 0.565217  | 0.521739  | … | 0.434783  | 0.434783  | 0.521739  | 0.391304 |
    | 4001      | 0.35      | 0.366667  | 0.383333  | … | 0.35      | 0.366667  | 0.4       | 0.4      |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 3/9 dataset performed better with less scrubbing
    for aCompCor processing, 5/9 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 3/9 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 4/9 dataset performed worst with aCompCor compared to gsr
    #### plot fc categories per denoising method ####

    looking at value effect in variable in mouse
    | variable         | No       | Non-specific | Specific | Spurious |
    |------------------|----------|--------------|----------|----------|
    | s1.cat.aCompCor1 | 1.606046 | 3.212093     | 3.826169 | 2.456306 |
    | s1.cat.aCompCor2 | 1.511573 | 3.353803     | 3.920642 | 2.361833 |
    | s1.cat.aCompCor3 | 1.653283 | 3.212093     | 4.109589 | 2.12565  |
    | s1.cat.gsr1      | 1.70052  | 2.361833     | 4.487482 | 2.550779 |
    | s1.cat.gsr2      | 1.70052  | 2.361833     | 4.487482 | 2.550779 |
    | s1.cat.gsr3      | 1.464336 | 2.267359     | 4.676429 | 2.692489 |
    | s1.cat.wmcsf1    | 1.464336 | 3.684459     | 3.306566 | 2.692489 |
    | s1.cat.wmcsf2    | 1.464336 | 3.164856     | 3.684459 | 2.786963 |
    | s1.cat.wmcsf3    | 1.27539  | 3.306566     | 4.015116 | 2.503543 |
    the effect of value on variable in mouse is q =  28.05 with p-value = 0.25774, dof = 24
    doing statistical analysis
    first with gsr3 processed data
    looking at fd.mean effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 8.860759  | 13.080169 | 10.548523 | 9.2827    |
    | other    | 16.455696 | 11.814346 | 14.345992 | 15.611814 |
    the effect of fd.mean on cat in rat is q =  4.44 with p-value = 0.21757, dof = 3
    looking at fd.max effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 12.236287 | 11.814346 | 10.548523 | 7.172996  |
    | other    | 13.080169 | 13.080169 | 14.345992 | 17.721519 |
    the effect of fd.max on cat in rat is q =  5.93 with p-value = 0.11516, dof = 3
    looking at s1.gsrcov.l.gsr3 effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 9.704641  | 8.860759  | 10.126582 | 13.080169 |
    | other    | 15.611814 | 16.033755 | 14.767932 | 11.814346 |
    the effect of s1.gsrcov.l.gsr3 on cat in rat is q =  4.06 with p-value = 0.25499, dof = 3
    looking at s1.tsnr.l effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 8.438819  | 10.970464 | 10.970464 | 11.392405 |
    | other    | 16.877637 | 13.924051 | 13.924051 | 13.50211  |
    the effect of s1.tsnr.l on cat in rat is q =  2.4 with p-value = 0.49389, dof = 3
    looking at habituation.min effect in cat in rat
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 21.097046 | 20.675105 |
    | other    | 29.957806 | 28.270042 |
    the effect of habituation.min on cat in rat is q =  0.0 with p-value = 0.99069, dof = 1
    looking at habituation.days effect in cat in rat
    | cat      | low       |
    |----------|-----------|
    | Specific | 41.772152 |
    | other    | 58.227848 |
    the effect of habituation.days on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at short.habituation effect in cat in rat
    | cat      | long      | short    |
    |----------|-----------|----------|
    | Specific | 40.084388 | 1.687764 |
    | other    | 54.852321 | 3.375527 |
    the effect of short.habituation on cat in rat is q =  0.09 with p-value = 0.7581, dof = 1
    looking at main.experimenter.gender effect in cat in rat
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 24.050633 | 17.721519 |
    | other    | 26.582278 | 31.64557  |
    the effect of main.experimenter.gender on cat in rat is q =  2.82 with p-value = 0.09315, dof = 1
    looking at rodent.sex effect in cat in rat
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 3.375527 | 38.396624 |
    | other    | 5.907173 | 52.320675 |
    the effect of rodent.sex on cat in rat is q =  0.1 with p-value = 0.75419, dof = 1
    looking at head-plate effect in cat in rat
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 21.097046 | 20.675105 |
    | other    | 29.957806 | 28.270042 |
    the effect of head-plate on cat in rat is q =  0.0 with p-value = 0.99069, dof = 1
    looking at body.restrained effect in cat in rat
    | cat      | y         |
    |----------|-----------|
    | Specific | 41.772152 |
    | other    | 58.227848 |
    the effect of body.restrained on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at anesthesia.before.acquisition effect in cat in rat
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 1.687764 | 40.084388 |
    | other    | 6.329114 | 51.898734 |
    the effect of anesthesia.before.acquisition on cat in rat is q =  2.78 with p-value = 0.09554, dof = 1
    looking at MRI.field.strength effect in cat in rat
    | cat      | 7.0       | 9.4       |
    |----------|-----------|-----------|
    | Specific | 35.021097 | 6.751055  |
    | other    | 47.257384 | 10.970464 |
    the effect of MRI.field.strength on cat in rat is q =  0.13 with p-value = 0.71869, dof = 1
    looking at fMRI.sequence effect in cat in rat
    | cat      | GE-EPI    | SE-EPI   |
    |----------|-----------|----------|
    | Specific | 40.506329 | 1.265823 |
    | other    | 55.274262 | 2.953586 |
    the effect of fMRI.sequence on cat in rat is q =  0.2 with p-value = 0.65727, dof = 1
