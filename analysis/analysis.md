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


  analysis_sub = ['gsr3', 'aCompCor3']
  for analysis in analysis_sub:
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
    totalling 1586 runs
    from 198 animals
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 1368/1586 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain      |
    |-----------|-----------|--------------|----------------|-------------|
    | 1001      | 4         | 4            | 4              | C57BL/6     |
    | 1002      | 8         | 8            | 8              | C57BL/6     |
    | 1003      | 12        | 12           | 12             | C57BL/6     |
    | 1004      | 21        | 5            | 21             | C57BL/6     |
    | 1005      | 207       | 13           | 176            | C57BL/6     |
    | 1006      | 20        | 10           | 20             | C57BL/6     |
    | 1007      | 14        | 7            | 14             | C57BL/6     |
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
    the datasets contained 1320 male runs and 154 female runs
    that corresponds to 10.45% females 
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
    | 0.033747 |
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
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=1.38, p=0.1711, dof=123.06
    t-test for head-plate max.fd > no head-plate mean.fd
    t=1.69, p=0.09299, dof=129.67
    #### tSNR ANALYISIS ####
    tSNR across all mouse datasets
    | s1.tsnr.l |
    |-----------|
    | 12.683403 |
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
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Spurious     | 0.429059 |
    | No           | 0.272169 |
    | Specific     | 0.213506 |
    | Non-specific | 0.085266 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Spurious     | 0.420873 |
    | No           | 0.283083 |
    | Specific     | 0.217599 |
    | Non-specific | 0.078445 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Spurious     | 0.403138 |
    | No           | 0.294679 |
    | Specific     | 0.223056 |
    | Non-specific | 0.079127 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.414052 |
    | No            | 0.271487 |
    | Specific      | 0.209413 |
    | Non-specific  | 0.105048 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.408595 |
    | No            | 0.283765 |
    | Specific      | 0.215553 |
    | Non-specific  | 0.092087 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Spurious      | 0.373806 |
    | No            | 0.287176 |
    | Specific      | 0.232606 |
    | Non-specific  | 0.106412 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.395634 |
    | No               | 0.257844 |
    | Specific         | 0.230559 |
    | Non-specific     | 0.115962 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.396317 |
    | No               | 0.259891 |
    | Specific         | 0.228513 |
    | Non-specific     | 0.11528  |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Spurious         | 0.343111 |
    | No               | 0.313779 |
    | Specific         | 0.232606 |
    | Non-specific     | 0.110505 |
    for gsr processing, heavy scrubbing lead to  313/1466 Specific scans
    for gsr processing, light scrubbing lead to  327/1466 Specific scans
    for aCompCor processing, heavy scrubbing lead to  338/1466 Specific scans
    for aCompCor processing, light scrubbing lead to  341/1466 Specific scans
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
    | s1.cat.aCompCor1 | 2.864939 | 1.288464     | 2.561771 | 4.395938 |
    | s1.cat.aCompCor2 | 2.887676 | 1.280885     | 2.539033 | 4.403517 |
    | s1.cat.aCompCor3 | 3.486433 | 1.227831     | 2.584508 | 3.812339 |
    | s1.cat.gsr1      | 3.024102 | 0.9474       | 2.37229  | 4.767318 |
    | s1.cat.gsr2      | 3.145369 | 0.871608     | 2.417766 | 4.676368 |
    | s1.cat.gsr3      | 3.274216 | 0.879188     | 2.478399 | 4.479309 |
    | s1.cat.wmcsf1    | 3.016523 | 1.167197     | 2.326815 | 4.600576 |
    | s1.cat.wmcsf2    | 3.152948 | 1.023192     | 2.395028 | 4.539942 |
    | s1.cat.wmcsf3    | 3.190844 | 1.182356     | 2.584508 | 4.153403 |
    the effect of value on variable in mouse is q =  64.39 with p-value = 1e-05, dof = 24
    doing statistical analysis
    first with gsr3 processed data
    looking at fd.mean effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 4.576503  | 5.464481  | 5.532787  | 6.762295  |
    | other    | 20.423497 | 19.535519 | 19.467213 | 18.237705 |
    the effect of fd.mean on cat in mouse is q =  8.17 with p-value = 0.04262, dof = 3
    looking at fd.max effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 3.620219  | 7.308743  | 5.737705  | 5.669399  |
    | other    | 21.379781 | 17.691257 | 19.262295 | 19.330601 |
    the effect of fd.max on cat in mouse is q =  23.16 with p-value = 4e-05, dof = 3
    looking at s1.gsrcov.l.gsr3 effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 5.821918  | 4.657534  | 6.164384  | 5.684932  |
    | other    | 19.794521 | 19.726027 | 18.835616 | 19.315068 |
    the effect of s1.gsrcov.l.gsr3 on cat in mouse is q =  3.35 with p-value = 0.3408, dof = 3
    looking at s1.tsnr.l effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 6.830601  | 4.713115  | 4.371585  | 6.420765  |
    | other    | 18.169399 | 20.286885 | 20.628415 | 18.579235 |
    the effect of s1.tsnr.l on cat in mouse is q =  15.13 with p-value = 0.00171, dof = 3
    looking at habituation.min effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 13.369714 | 8.93588   |
    | other    | 47.885402 | 29.809004 |
    the effect of habituation.min on cat in mouse is q =  0.24 with p-value = 0.62419, dof = 1
    looking at habituation.days effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 13.506139 | 8.799454  |
    | other    | 48.158254 | 29.536153 |
    the effect of habituation.days on cat in mouse is q =  0.16 with p-value = 0.68511, dof = 1
    looking at short.habituation effect in cat in mouse
    | cat      | long      | short     |
    |----------|-----------|-----------|
    | Specific | 8.799454  | 13.506139 |
    | other    | 29.536153 | 48.158254 |
    the effect of short.habituation on cat in mouse is q =  0.16 with p-value = 0.68511, dof = 1
    looking at main.experimenter.gender effect in cat in mouse
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 11.869031 | 10.436562 |
    | other    | 31.377899 | 46.316508 |
    the effect of main.experimenter.gender on cat in mouse is q =  16.51 with p-value = 5e-05, dof = 1
    looking at rodent.sex effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 3.914328 | 18.020679 |
    | other    | 8.271787 | 69.793205 |
    the effect of rodent.sex on cat in mouse is q =  10.72 with p-value = 0.00106, dof = 1
    looking at head-plate effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 3.61528  | 18.690314 |
    | other    | 3.819918 | 73.874488 |
    the effect of head-plate on cat in mouse is q =  45.44 with p-value = 0.0, dof = 1
    looking at body.restrained effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 7.162347  | 15.143247 |
    | other    | 21.077763 | 56.616644 |
    the effect of body.restrained on cat in mouse is q =  2.87 with p-value = 0.09027, dof = 1
    looking at anesthesia.before.acquisition effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 3.13779   | 19.167804 |
    | other    | 14.665757 | 63.028649 |
    the effect of anesthesia.before.acquisition on cat in mouse is q =  3.69 with p-value = 0.05464, dof = 1
    looking at MRI.field.strength effect in cat in mouse
    | cat      | 7.0      | 9.4       | 11.7     | 15.2     |
    |----------|----------|-----------|----------|----------|
    | Specific | 0.886767 | 17.530696 | 2.933151 | 0.95498  |
    | other    | 1.432469 | 69.986357 | 3.956344 | 2.319236 |
    the effect of MRI.field.strength on cat in mouse is q =  34.05 with p-value = 0.0, dof = 3
    looking at fMRI.sequence effect in cat in mouse
    | cat      | GE-EPI    | SE-EPI    |
    |----------|-----------|-----------|
    | Specific | 9.41337   | 12.892224 |
    | other    | 28.854025 | 48.840382 |
    the effect of fMRI.sequence on cat in mouse is q =  2.55 with p-value = 0.11045, dof = 1
    doing statistical analysis
    first with aCompCor3 processed data
    looking at fd.mean effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 4.644809  | 6.284153  | 5.737705  | 6.625683  |
    | other    | 20.355191 | 18.715847 | 19.262295 | 18.374317 |
    the effect of fd.mean on cat in mouse is q =  7.38 with p-value = 0.06066, dof = 3
    looking at fd.max effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 4.986339  | 6.147541  | 6.010929  | 6.147541  |
    | other    | 20.013661 | 18.852459 | 18.989071 | 18.852459 |
    the effect of fd.max on cat in mouse is q =  3.1 with p-value = 0.37639, dof = 3
    looking at s1.gsrcov.l.aCompCor3 effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 5.129959  | 6.908345  | 5.745554  | 5.471956  |
    | other    | 19.904241 | 18.057456 | 19.220246 | 19.562244 |
    the effect of s1.gsrcov.l.aCompCor3 on cat in mouse is q =  5.97 with p-value = 0.11325, dof = 3
    looking at s1.tsnr.l effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 6.420765  | 5.737705  | 4.030055  | 7.103825  |
    | other    | 18.579235 | 19.262295 | 20.969945 | 17.896175 |
    the effect of s1.tsnr.l on cat in mouse is q =  17.11 with p-value = 0.00067, dof = 3
    looking at habituation.min effect in cat in mouse
    | cat      | low       | high     |
    |----------|-----------|----------|
    | Specific | 14.461119 | 8.799454 |
    | other    | 46.793997 | 29.94543 |
    the effect of habituation.min on cat in mouse is q =  0.11 with p-value = 0.73954, dof = 1
    looking at habituation.days effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 14.665757 | 8.594816  |
    | other    | 46.998636 | 29.740791 |
    the effect of habituation.days on cat in mouse is q =  0.29 with p-value = 0.59119, dof = 1
    looking at short.habituation effect in cat in mouse
    | cat      | long      | short     |
    |----------|-----------|-----------|
    | Specific | 8.594816  | 14.665757 |
    | other    | 29.740791 | 46.998636 |
    the effect of short.habituation on cat in mouse is q =  0.29 with p-value = 0.59119, dof = 1
    looking at main.experimenter.gender effect in cat in mouse
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 10.436562 | 12.824011 |
    | other    | 32.810368 | 43.929059 |
    the effect of main.experimenter.gender on cat in mouse is q =  0.39 with p-value = 0.53041, dof = 1
    looking at rodent.sex effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 3.766617 | 19.57164  |
    | other    | 8.419498 | 68.242245 |
    the effect of rodent.sex on cat in mouse is q =  5.55 with p-value = 0.01851, dof = 1
    looking at head-plate effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 2.18281  | 21.077763 |
    | other    | 5.252387 | 71.48704  |
    the effect of head-plate on cat in mouse is q =  2.1 with p-value = 0.14755, dof = 1
    looking at body.restrained effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 7.503411  | 15.757162 |
    | other    | 20.736698 | 56.002729 |
    the effect of body.restrained on cat in mouse is q =  3.29 with p-value = 0.06986, dof = 1
    looking at anesthesia.before.acquisition effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 3.001364  | 20.259209 |
    | other    | 14.802183 | 61.937244 |
    the effect of anesthesia.before.acquisition on cat in mouse is q =  6.86 with p-value = 0.00881, dof = 1
    looking at MRI.field.strength effect in cat in mouse
    | cat      | 7.0      | 9.4       | 11.7     | 15.2     |
    |----------|----------|-----------|----------|----------|
    | Specific | 0.613915 | 19.372442 | 2.319236 | 0.95498  |
    | other    | 1.705321 | 68.144611 | 4.570259 | 2.319236 |
    the effect of MRI.field.strength on cat in mouse is q =  8.17 with p-value = 0.04268, dof = 3
    looking at fMRI.sequence effect in cat in mouse
    | cat      | GE-EPI    | SE-EPI    |
    |----------|-----------|-----------|
    | Specific | 9.345157  | 13.915416 |
    | other    | 28.922237 | 47.81719  |
    the effect of fMRI.sequence on cat in mouse is q =  0.58 with p-value = 0.44477, dof = 1
    #### NOW DOING rat ####
    summary of the data that we collected
    we processed 9 datasets
    totalling 468 runs
    from 161 animals
    the smallest dataset had 5 runs
    the largest dataset had 291 runs
    we could processed 237/468 runs.
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
    | 4001      | 291       | 89           | 60             | Long-Evans          |
    information about sex ratio
    the datasets contained 446 male runs and 22 female runs
    that corresponds to 4.7% females 
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
    | 0.041656 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 2001      | 0.03643  |
    | 2002a     | 0.054875 |
    | 2002b     | 0.044527 |
    | 2002c     | 0.044339 |
    | 2003      | 0.04796  |
    | 2004      | 0.016877 |
    | 2005      | 0.018197 |
    | 2006      | 0.045115 |
    | 4001      | 0.039548 |
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=0.74, p=0.46141, dof=204.53
    t-test for head-plate max.fd > no head-plate mean.fd
    t=-6.45, p=0.0, dof=196.61
    #### tSNR ANALYISIS ####
    tSNR across all rat datasets
    | s1.tsnr.l |
    |-----------|
    | 28.134013 |
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
    | 2002a     | 1500.0       | 467.666667          | 289.5               | 5.666667            |
    | 2002b     | 750.0        | 144.4               | 93.0                | 6.0                 |
    | 2002c     | 750.0        | 141.8               | 101.0               | 6.0                 |
    | 2003      | 90.0         | 21.092784           | 9.185567            | 10.5                |
    | 2004      | 400.0        | 80.105263           | 79.157895           | null                |
    | 2005      | 1000.0       | 181.2               | 102.2               | 3.0                 |
    | 2006      | 750.0        | 137.173913          | 79.130435           | 11.3                |
    | 4001      | 600.0        | 103.85              | 79.716667           | 13.151515           |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Specific     | 0.396624 |
    | Spurious     | 0.232068 |
    | Non-specific | 0.21097  |
    | No           | 0.151899 |
    | null         | 0.008439 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Specific     | 0.396624 |
    | Spurious     | 0.232068 |
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
    | Non-specific  | 0.324895 |
    | Specific      | 0.295359 |
    | Spurious      | 0.244726 |
    | No            | 0.130802 |
    | null          | 0.004219 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Specific      | 0.329114 |
    | Non-specific  | 0.278481 |
    | Spurious      | 0.253165 |
    | No            | 0.130802 |
    | null          | 0.008439 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Specific      | 0.35865  |
    | Non-specific  | 0.291139 |
    | Spurious      | 0.227848 |
    | No            | 0.113924 |
    | null          | 0.008439 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Specific         | 0.341772 |
    | Non-specific     | 0.2827   |
    | Spurious         | 0.223629 |
    | No               | 0.14346  |
    | null             | 0.008439 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Specific         | 0.345992 |
    | Non-specific     | 0.299578 |
    | Spurious         | 0.21519  |
    | No               | 0.135021 |
    | null             | 0.004219 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Specific         | 0.362869 |
    | Non-specific     | 0.28692  |
    | Spurious         | 0.194093 |
    | No               | 0.147679 |
    | null             | 0.008439 |
    for gsr processing, heavy scrubbing lead to  94/235 Specific scans
    for gsr processing, light scrubbing lead to  99/235 Specific scans
    for aCompCor processing, heavy scrubbing lead to  81/235 Specific scans
    for aCompCor processing, light scrubbing lead to  86/235 Specific scans
    s1-s1 specificity analysis per dataset (only showing specific values)
    | rodent.ds | s1.gsr1.S | s1.gsr2.S | s1.gsr3.S | … | s1.wmcsf3 | s1.aCompC | s1.aCompC | s1.aComp |
    |           | pecific   | pecific   | pecific   |   | .Specific | or1.Speci | or2.Speci | Cor3.Spe |
    |           |           |           |           |   |           | fic       | fic       | cific    |
    |-----------|-----------|-----------|-----------|---|-----------|-----------|-----------|----------|
    | 2001      | 0.625     | 0.625     | 0.625     | … | 0.4375    | 0.625     | 0.5       | 0.5625   |
    | 2002a     | 0.571429  | 0.428571  | 0.285714  | … | 0.142857  | 0.285714  | 0.142857  | 0.142857 |
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
    for comparing gsr to aCompCor processing, 5/9 dataset performed worst with aCompCor compared to gsr
    #### plot fc categories per denoising method ####
    looking at value effect in variable in mouse
    | variable         | No       | Non-specific | Specific | Spurious |
    |------------------|----------|--------------|----------|----------|
    | s1.cat.aCompCor1 | 1.606046 | 3.164856     | 3.826169 | 2.503543 |
    | s1.cat.aCompCor2 | 1.511573 | 3.353803     | 3.873406 | 2.409069 |
    | s1.cat.aCompCor3 | 1.653283 | 3.212093     | 4.062352 | 2.172886 |
    | s1.cat.gsr1      | 1.70052  | 2.361833     | 4.440246 | 2.598016 |
    | s1.cat.gsr2      | 1.70052  | 2.361833     | 4.440246 | 2.598016 |
    | s1.cat.gsr3      | 1.464336 | 2.267359     | 4.676429 | 2.692489 |
    | s1.cat.wmcsf1    | 1.464336 | 3.637222     | 3.306566 | 2.739726 |
    | s1.cat.wmcsf2    | 1.464336 | 3.117619     | 3.684459 | 2.834199 |
    | s1.cat.wmcsf3    | 1.27539  | 3.259329     | 4.015116 | 2.550779 |
    the effect of value on variable in mouse is q =  26.78 with p-value = 0.31488, dof = 24
    doing statistical analysis
    first with gsr3 processed data
    looking at fd.mean effect in cat in rat
    | cat      | lowest    | low       | high     | highest   |
    |----------|-----------|-----------|----------|-----------|
    | Specific | 8.474576  | 13.559322 | 10.59322 | 9.322034  |
    | other    | 16.525424 | 11.440678 | 14.40678 | 15.677966 |
    the effect of fd.mean on cat in rat is q =  5.76 with p-value = 0.12392, dof = 3
    looking at fd.max effect in cat in rat
    | cat      | lowest    | low       | high      | highest  |
    |----------|-----------|-----------|-----------|----------|
    | Specific | 11.864407 | 11.864407 | 11.016949 | 7.20339  |
    | other    | 13.135593 | 13.135593 | 13.983051 | 17.79661 |
    the effect of fd.max on cat in rat is q =  5.76 with p-value = 0.12392, dof = 3
    looking at s1.gsrcov.l.gsr3 effect in cat in rat
    | cat      | lowest    | low       | high     | highest   |
    |----------|-----------|-----------|----------|-----------|
    | Specific | 9.745763  | 8.474576  | 10.59322 | 13.135593 |
    | other    | 15.254237 | 16.525424 | 14.40678 | 11.864407 |
    the effect of s1.gsrcov.l.gsr3 on cat in rat is q =  4.51 with p-value = 0.2117, dof = 3
    looking at s1.tsnr.l effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 8.050847  | 11.016949 | 11.440678 | 11.440678 |
    | other    | 16.949153 | 13.983051 | 13.559322 | 13.559322 |
    the effect of s1.tsnr.l on cat in rat is q =  3.11 with p-value = 0.37428, dof = 3
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
    doing statistical analysis
    first with aCompCor3 processed data
    looking at fd.mean effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 8.474576  | 11.864407 | 7.627119  | 8.474576  |
    | other    | 16.525424 | 13.135593 | 17.372881 | 16.525424 |
    the effect of fd.mean on cat in rat is q =  4.32 with p-value = 0.22916, dof = 3
    looking at fd.max effect in cat in rat
    | cat      | lowest    | low       | high      | highest  |
    |----------|-----------|-----------|-----------|----------|
    | Specific | 11.864407 | 9.745763  | 7.627119  | 7.20339  |
    | other    | 13.135593 | 15.254237 | 17.372881 | 17.79661 |
    the effect of fd.max on cat in rat is q =  5.63 with p-value = 0.1308, dof = 3
    looking at s1.gsrcov.l.aCompCor3 effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 9.322034  | 8.474576  | 11.016949 | 7.627119  |
    | other    | 15.677966 | 16.525424 | 13.983051 | 17.372881 |
    the effect of s1.gsrcov.l.aCompCor3 on cat in rat is q =  2.56 with p-value = 0.46433, dof = 3
    looking at s1.tsnr.l effect in cat in rat
    | cat      | lowest    | low       | high      | highest  |
    |----------|-----------|-----------|-----------|----------|
    | Specific | 7.627119  | 9.745763  | 8.474576  | 10.59322 |
    | other    | 17.372881 | 15.254237 | 16.525424 | 14.40678 |
    the effect of s1.tsnr.l on cat in rat is q =  2.12 with p-value = 0.54744, dof = 3
    looking at habituation.min effect in cat in rat
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 20.253165 | 16.033755 |
    | other    | 30.801688 | 32.911392 |
    the effect of habituation.min on cat in rat is q =  0.94 with p-value = 0.33157, dof = 1
    looking at habituation.days effect in cat in rat
    | cat      | low      |
    |----------|----------|
    | Specific | 36.28692 |
    | other    | 63.71308 |
    the effect of habituation.days on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at short.habituation effect in cat in rat
    | cat      | long      | short    |
    |----------|-----------|----------|
    | Specific | 35.443038 | 0.843882 |
    | other    | 59.493671 | 4.219409 |
    the effect of short.habituation on cat in rat is q =  1.31 with p-value = 0.25318, dof = 1
    looking at main.experimenter.gender effect in cat in rat
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 18.565401 | 17.721519 |
    | other    | 32.067511 | 31.64557  |
    the effect of main.experimenter.gender on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
    looking at rodent.sex effect in cat in rat
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 4.219409 | 32.067511 |
    | other    | 5.063291 | 58.649789 |
    the effect of rodent.sex on cat in rat is q =  0.5 with p-value = 0.48009, dof = 1
    looking at head-plate effect in cat in rat
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 20.253165 | 16.033755 |
    | other    | 30.801688 | 32.911392 |
    the effect of head-plate on cat in rat is q =  0.94 with p-value = 0.33157, dof = 1
    looking at body.restrained effect in cat in rat
    | cat      | y        |
    |----------|----------|
    | Specific | 36.28692 |
    | other    | 63.71308 |
    the effect of body.restrained on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at anesthesia.before.acquisition effect in cat in rat
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 1.687764 | 34.599156 |
    | other    | 6.329114 | 57.383966 |
    the effect of anesthesia.before.acquisition on cat in rat is q =  1.42 with p-value = 0.23356, dof = 1
    looking at MRI.field.strength effect in cat in rat
    | cat      | 7.0       | 9.4       |
    |----------|-----------|-----------|
    | Specific | 30.801688 | 5.485232  |
    | other    | 51.476793 | 12.236287 |
    the effect of MRI.field.strength on cat in rat is q =  0.38 with p-value = 0.53805, dof = 1
    looking at fMRI.sequence effect in cat in rat
    | cat      | GE-EPI    | SE-EPI   |
    |----------|-----------|----------|
    | Specific | 34.599156 | 1.687764 |
    | other    | 61.181435 | 2.531646 |
    the effect of fMRI.sequence on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
