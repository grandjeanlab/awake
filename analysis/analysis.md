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
  df = pl.read_csv("../assets/tables/"+rodent+"_metadata_processed.tsv", separator="\t")
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

  ax = sns.barplot(x=df_summary["rodent.ds"], y=df_summary["total_included"], hue=df_summary["strain"], dodge=False)
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  plt.savefig("../assets/plot/"+rodent+"_run_included_per_dataset.svg")

  print("information about sex ratio")
  print("the datasets contained "+str(df_summary["male"].sum())+ " male runs and " +str(df_summary["female"].sum())+ " female runs")
  print("that corresponds to "+str(round(100*df_summary["female"].sum()/(df_summary["female"].sum()+df_summary["male"].sum()),2))+"% females ")
  print("information about animal handling")
  print(str((df_summary["headplate"] == 'y').sum()) + " datasets used headplates")
  print(str((df_summary["restrained"] == 'y').sum()) + " datasets used body restraining")
  print(str((df_summary["anesthesia"].is_in(['Isoflurane', 'Sevoflurane'])).sum()) + " datasets used anesthesia before acquisition")
  print(str((df_summary["exp.gender"] == 'm').sum()) + " datasets were collected by men, " + str((df_summary["exp.gender"] == 'f').sum())+ " by women")

  print(df_summary.select("rodent.ds", "headplate", "restrained", "anesthesia","exp.gender", "habituation.days","habituation.min"))

  fig, ax1 = plt.subplots()
  ax2 = ax1.twinx()
  ax = sns.barplot(x=df_summary["rodent.ds"], y=df_summary["habituation.min"], hue=df_summary["headplate"], dodge=False, ax=ax1)
  ax = sns.swarmplot(x=df_summary["rodent.ds"], y=df_summary["habituation.days"], color='black', marker='o', label='habituation days', legend=False, ax=ax2)
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  plt.savefig("../assets/plot/"+rodent+"_habituation_per_dataset.svg")

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

  ax = makeswarmplot('framewise displacement per dataset', df["rodent.ds"], df["fd.mean"], hue=df['head-plate'])
  ax.figure.savefig("../assets/plot/"+rodent+"_fd_per_dataset.svg")

  df_to_plot = df.select('head-plate','fd.mean').rename({'head-plate':'value','fd.mean':'cont_variable'})
  ax = makeviolinplot(df_to_plot)
  ax.figure.savefig("../assets/plot/"+rodent+"_fd_headplate_violin.svg")
  
  t = ttest(df.filter(pl.col('head-plate')=='y')['fd.mean'], df.filter(pl.col('head-plate')=='n')['fd.mean'])
  print('t-test for head-plate mean.fd > no head-plate mean.fd') 
  print(f't={round(t["T"].item(),2)}, p={round(t["p-val"].item(),5)}, dof={round(t["dof"].item(),2)}')



#let's extract some infomation about tsnr and summarize it per dataset
  print("#### tSNR ANALYISIS ####")
  print("tSNR across all "+rodent+" datasets")
  print(df.select("s1.tsnr.l").mean())

  ax = makeswarmplot('S1 tSNR per dataset', df["rodent.ds"], df["s1.tsnr.l"], hue=df['MRI.field.strength'].cast(pl.String))
  ax.figure.savefig("../assets/plot/"+rodent+"_tsnr_per_dataset.svg")

  df_to_plot = df.select('MRI.field.strength','s1.tsnr.l').rename({'MRI.field.strength':'value','s1.tsnr.l':'cont_variable'}).sort(by='value').with_columns(pl.col('value').cast(pl.String))
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
  ax = makestackplot(df_melted)
  plt.savefig("../assets/plot/"+rodent+"_specificity.svg")

  analysis_sub = ['gsr3', 'aCompCor3']
  for analysis in analysis_sub:
    ax = makespecificityplot(rodent+': S1 specificity with ' + analysis, df['s1.specific.'+analysis], df['s1.unspecific.'+analysis])
    plt.savefig("../assets/plot/specific/"+rodent+"_s1_specificity_"+analysis+".svg")

    ax = makespecificityplot(rodent+': Thal FC specificity with ' + analysis, df['thal.specific.'+analysis], df['s1.unspecific.'+analysis])
    plt.savefig("../assets/plot/specific/"+rodent+"_thal_specificity_"+analysis+".svg")

    cat1='s1.cat.'+analysis
    cat2 = 'fd.mean'
    df_to_plot = split_continuous(df, cat1, cat2)
    df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
    ax = makestackplot(df_to_plot)
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_fd.svg")

    ax = makeviolinplot(df_to_plot)
    ax.set_xlabel('Connectivity category')
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_fd_violin.svg")

    cat2='s1.gsrcov.l.'+analysis
    df_to_plot = split_continuous(df, cat1, cat2)
    df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
    ax = makestackplot(df_to_plot)
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_gsrcov.svg")

    ax = makeviolinplot(df_to_plot)
    ax.set_xlabel('Connectivity category')
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_gsrcov_violin.svg")

    cat2='s1.tsnr.l'
    df_to_plot = split_continuous(df, cat1, cat2)
    df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
    ax = makestackplot(df_to_plot)
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_tsnr.svg")

    ax = makeviolinplot(df_to_plot)
    ax.set_xlabel('Connectivity category')
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_tsnr_violin.svg")

    cat2='habituation.min'
    df_to_plot = split_continuous(df, cat1, cat2)
    df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
    ax = makeviolinplot(df_to_plot)
    ax.set_xlabel('Connectivity category')
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_habituation_min.svg")

    cat2='habituation.days'
    df_to_plot = split_continuous(df, cat1, cat2)
    df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
    ax = makeviolinplot(df_to_plot)
    ax.set_xlabel('Connectivity category')
    plt.savefig("../assets/plot/"+rodent+"_"+analysis+"_habituation_days.svg")

    df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
    cat2='main.experimenter.gender'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_expgender.svg')

    cat2='rodent.sex'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_sex.svg')

    cat2='head-plate'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_headplate.svg')

    cat2='body.restrained'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_restrained.svg')

    cat2='short.habituation'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_shorthabituation.svg')

    cat2='anesthesia.before.acquisition'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_anesthesia.svg')

    cat2='MRI.field.strength'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_MRI.svg')

    cat2='fMRI.sequence'
    df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
    ax = sns.barplot(x=df_to_plot['cat'],y=df_to_plot['len'],hue=df_to_plot[cat2])
    plt.savefig('../assets/plot/'+rodent+'_'+analysis+'_fMRI.svg')

    print('doing statistical analysis')
    print('first with '+analysis+' processed data')
    cat2 = 'fd.mean'
    chi2_continuous(df, 'cat', cat2, rodent)
    cat2 = 's1.gsrcov.l.'+analysis
    chi2_continuous(df, 'cat', cat2, rodent)
    cat2 = 's1.tsnr.l'
    chi2_continuous(df, 'cat', cat2, rodent)
    cat2 = 'habituation.min'
    chi2_continuous(df, 'cat', cat2, rodent, 3, ['low','intermediate','high'])
    cat2 = 'habituation.days'
    chi2_continuous(df, 'cat', cat2, rodent, 3, ['low','intermediate','high'])

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
    we processed 19 datasets
    totalling 1344 runs
    from 177 animals
    the smallest dataset had 4 runs
    the largest dataset had 479 runs
    we could processed 1157/1344 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain      |
    |-----------|-----------|--------------|----------------|-------------|
    | 1001      | 4         | 4            | 4              | C57BL/6     |
    | 1002      | 8         | 8            | 8              | C57BL/6     |
    | 1003      | 12        | 12           | 12             | C57BL/6     |
    | 1004      | 21        | 5            | 21             | C57BL/6     |
    | 1005a     | 13        | 4            | 13             | C57BL/6     |
    | 1006      | 20        | 10           | 20             | C57BL/6     |
    | 1007      | 14        | 7            | 14             | C57BL/6     |
    | 1008      | 36        | 9            | 34             | C57BL/6     |
    | 1009      | 107       | 5            | 34             | 129S2/SvPas |
    | 1011      | 51        | 17           | 48             | C57BL/6     |
    | 1012      | 26        | 26           | 21             | null        |
    | 1013      | 12        | 3            | 12             | C57BL/6     |
    | 1014      | 4         | 4            | 4              | null        |
    | 1015      | 6         | 6            | 6              | C57BL/6J    |
    | 3001      | 112       | 6            | 112            | C57BL/6     |
    | 3002      | 10        | 10           | 10             | C57BL/6     |
    | 3003      | 479       | 19           | 414            | F1 C6/129P  |
    | 3004      | 54        | 9            | 54             | C57BL/6     |
    | 3005      | 355       | 13           | 316            | C57BL/6     |

    information about sex ratio
    the datasets contained 1082 male runs and 150 female runs
    that corresponds to 12.18% females 
    information about animal handling
    16 datasets used headplates
    13 datasets used body restraining
    0 datasets used anesthesia before acquisition
    11 datasets were collected by men, 8 by women
    | rodent.ds | headplate | restrained | anesthesia | exp.gender | habituation.day | habituation.min |
    |           |           |            |            |            | s               |                 |
    |-----------|-----------|------------|------------|------------|-----------------|-----------------|
    | 1001      | y         | n          | n          | m          | 7               | 315             |
    | 1002      | y         | n          | n          | m          | 0               | 0               |
    | 1003      | y         | y          | n          | m          | 19              | 695             |
    | 1004      | y         | y          | n          | f          | 10              | 233             |
    | 1005a     | y         | y          | y          | f          | 7               | 1320            |
    | 1006      | y         | n          | y          | m          | 7               | 260             |
    | 1007      | y         | y          | n          | m          | 9               | 515             |
    | 1008      | n         | y          | n          | m          | 8               | 420             |
    | 1009      | y         | y          | y          | m          | 10              | 50              |
    | 1011      | y         | y          | y          | m          | 0               | 0               |
    | 1012      | n         | y          | y          | f          | 5               | 150             |
    | 1013      | y         | n          | n          | f          | 5               | 100             |
    | 1014      | y         | y          | y          | f          | 7               | 180             |
    | 1015      | y         | n          | y          | m          | 4               | 300             |
    | 3001      | y         | y          | y          | m          | 4               | 42              |
    | 3002      | y         | y          | n          | f          | 16              | 432             |
    | 3003      | y         | y          | y          | m          | 4               | 42              |
    | 3004      | n         | y          | y          | f          | 5               | 165             |
    | 3005      | y         | n          | y          | f          | 4               | 42              |
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
    | 1015      | 11.7           | GE-EPI   | 14.0   |
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
    | 1015      | null     |
    | 3001      | 0.035982 |
    | 3002      | 0.022321 |
    | 3003      | 0.040797 |
    | 3004      | 0.03515  |
    | 3005      | 0.046262 |
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=3.22, p=0.00165, dof=123.93
    #### tSNR ANALYISIS ####
    tSNR across all mouse datasets
    | s1.tsnr.l |
    |-----------|
    | 11.651712 |
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
    | 1015      | null         | null                | null                | null                |
    | 3001      | 200.0        | 38.142857           | 28.821429           | 4.404255            |
    | 3002      | 1920.0       | 400.2               | 388.0               | 3.0                 |
    | 3003      | 200.0        | 55.074879           | 41.342995           | 4.953608            |
    | 3004      | 180.0        | 51.296296           | 39.150943           | 7.75                |
    | 3005      | 200.0        | 52.234177           | 38.547468           | 4.965986            |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Spurious     | 0.457371 |
    | No           | 0.238247 |
    | Specific     | 0.211952 |
    | Non-specific | 0.09243  |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Spurious     | 0.447809 |
    | No           | 0.250996 |
    | Specific     | 0.21753  |
    | Non-specific | 0.083665 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Spurious     | 0.426295 |
    | No           | 0.271713 |
    | Specific     | 0.223108 |
    | Non-specific | 0.078884 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.446215 |
    | No            | 0.23745  |
    | Specific      | 0.203187 |
    | Non-specific  | 0.113147 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.439044 |
    | No            | 0.253386 |
    | Specific      | 0.209562 |
    | Non-specific  | 0.098008 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Spurious      | 0.4      |
    | No            | 0.264542 |
    | Specific      | 0.224701 |
    | Non-specific  | 0.110757 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.41992  |
    | No               | 0.231873 |
    | Specific         | 0.231076 |
    | Non-specific     | 0.117131 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.416733 |
    | No               | 0.23506  |
    | Specific         | 0.230279 |
    | Non-specific     | 0.117928 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Spurious         | 0.36255  |
    | No               | 0.290837 |
    | Specific         | 0.232669 |
    | Non-specific     | 0.113944 |
    for gsr processing, heavy scrubbing lead to  266/1255 Specific scans
    for gsr processing, light scrubbing lead to  280/1255 Specific scans
    for aCompCor processing, heavy scrubbing lead to  290/1255 Specific scans
    for aCompCor processing, light scrubbing lead to  292/1255 Specific scans
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
    | 1015      | 0.0       | 0.0       | 0.0       | … | 0.0       | 0.0       | 0.0       | 0.0      |
    | 3001      | 0.303571  | 0.3125    | 0.267857  | … | 0.214286  | 0.303571  | 0.294643  | 0.223214 |
    | 3002      | 0.5       | 0.5       | 0.7       | … | 0.7       | 0.1       | 0.0       | 0.1      |
    | 3003      | 0.149758  | 0.154589  | 0.161836  | … | 0.178744  | 0.207729  | 0.198068  | 0.21256  |
    | 3004      | 0.592593  | 0.62963   | 0.611111  | … | 0.555556  | 0.407407  | 0.388889  | 0.37037  |
    | 3005      | 0.21519   | 0.212025  | 0.253165  | … | 0.272152  | 0.231013  | 0.224684  | 0.246835 |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 11/19 dataset performed better with less scrubbing
    for aCompCor processing, 7/19 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 7/19 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 8/19 dataset performed worst with aCompCor compared to gsr
    #### plot fc categories per denoising method ####
    doing statistical analysis
    first with gsr3 processed data
    looking at fd.mean effect in cat in mouse
    | cat      | lowest    | low   | high      | highest   |
    |----------|-----------|-------|-----------|-----------|
    | Specific | 4.166667  | 6.25  | 5.128205  | 6.891026  |
    | other    | 20.833333 | 18.75 | 19.871795 | 18.108974 |
    the effect of fd.mean on cat in mouse is q =  12.52 with p-value = 0.00579, dof = 3
    looking at s1.gsrcov.l.gsr3 effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 5.707395  | 4.581994  | 6.993569  | 5.144695  |
    | other    | 19.292605 | 20.418006 | 18.006431 | 19.855305 |
    the effect of s1.gsrcov.l.gsr3 on cat in mouse is q =  9.14 with p-value = 0.02744, dof = 3
    looking at s1.tsnr.l effect in cat in mouse
    | cat      | lowest   | low       | high      | highest   |
    |----------|----------|-----------|-----------|-----------|
    | Specific | 7.13141  | 4.246795  | 4.567308  | 6.490385  |
    | other    | 17.86859 | 20.753205 | 20.432692 | 18.509615 |
    the effect of s1.tsnr.l on cat in mouse is q =  17.31 with p-value = 0.00061, dof = 3
    looking at habituation.min effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 15.61753  | 6.693227  |
    | other    | 55.936255 | 21.752988 |
    the effect of habituation.min on cat in mouse is q =  0.33 with p-value = 0.5628, dof = 1
    looking at habituation.days effect in cat in mouse
    | cat      | low       | high     |
    |----------|-----------|----------|
    | Specific | 15.61753  | 6.693227 |
    | other    | 56.414343 | 21.2749  |
    the effect of habituation.days on cat in mouse is q =  0.61 with p-value = 0.43311, dof = 1
    looking at short.habituation effect in cat in mouse
    | cat      | long     | short     |
    |----------|----------|-----------|
    | Specific | 6.693227 | 15.61753  |
    | other    | 21.2749  | 56.414343 |
    the effect of short.habituation on cat in mouse is q =  0.61 with p-value = 0.43311, dof = 1
    looking at main.experimenter.gender effect in cat in mouse
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 10.836653 | 11.474104 |
    | other    | 25.099602 | 52.589641 |
    the effect of main.experimenter.gender on cat in mouse is q =  24.29 with p-value = 0.0, dof = 1
    looking at rodent.sex effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 4.549431 | 17.322835 |
    | other    | 9.536308 | 68.591426 |
    the effect of rodent.sex on cat in mouse is q =  11.22 with p-value = 0.00081, dof = 1
    looking at head-plate effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 4.223108 | 18.087649 |
    | other    | 4.462151 | 73.227092 |
    the effect of head-plate on cat in mouse is q =  46.03 with p-value = 0.0, dof = 1
    looking at body.restrained effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 7.49004   | 14.820717 |
    | other    | 21.673307 | 56.015936 |
    the effect of body.restrained on cat in mouse is q =  3.12 with p-value = 0.07729, dof = 1
    looking at anesthesia.before.acquisition effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 2.948207  | 19.36255  |
    | other    | 14.023904 | 63.665339 |
    the effect of anesthesia.before.acquisition on cat in mouse is q =  3.28 with p-value = 0.07027, dof = 1
    looking at MRI.field.strength effect in cat in mouse
    | cat      | 7.0      | 9.4       | 11.7     | 15.2     |
    |----------|----------|-----------|----------|----------|
    | Specific | 1.035857 | 16.89243  | 3.266932 | 1.115538 |
    | other    | 1.673307 | 68.525896 | 4.780876 | 2.709163 |
    the effect of MRI.field.strength on cat in mouse is q =  29.73 with p-value = 0.0, dof = 3
    looking at fMRI.sequence effect in cat in mouse
    | cat      | GE-EPI   | SE-EPI    |
    |----------|----------|-----------|
    | Specific | 8.207171 | 14.103586 |
    | other    | 24.38247 | 53.306773 |
    the effect of fMRI.sequence on cat in mouse is q =  2.65 with p-value = 0.10369, dof = 1
    doing statistical analysis
    first with aCompCor3 processed data
    looking at fd.mean effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 4.887821  | 6.730769  | 5.208333  | 6.570513  |
    | other    | 20.112179 | 18.269231 | 19.791667 | 18.429487 |
    the effect of fd.mean on cat in mouse is q =  7.33 with p-value = 0.06204, dof = 3
    looking at s1.gsrcov.l.aCompCor3 effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 5.778491  | 6.340289  | 5.617978  | 5.617978  |
    | other    | 19.261637 | 18.619583 | 19.341894 | 19.422151 |
    the effect of s1.gsrcov.l.aCompCor3 on cat in mouse is q =  1.01 with p-value = 0.79803, dof = 3
    looking at s1.tsnr.l effect in cat in mouse
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 6.410256  | 5.448718  | 4.887821  | 6.650641  |
    | other    | 18.589744 | 19.551282 | 20.112179 | 18.349359 |
    the effect of s1.tsnr.l on cat in mouse is q =  5.69 with p-value = 0.12789, dof = 3
    looking at habituation.min effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 16.89243  | 6.374502  |
    | other    | 54.661355 | 22.071713 |
    the effect of habituation.min on cat in mouse is q =  0.14 with p-value = 0.70431, dof = 1
    looking at habituation.days effect in cat in mouse
    | cat      | low       | high      |
    |----------|-----------|-----------|
    | Specific | 16.89243  | 6.374502  |
    | other    | 55.139442 | 21.593625 |
    the effect of habituation.days on cat in mouse is q =  0.03 with p-value = 0.86211, dof = 1
    looking at short.habituation effect in cat in mouse
    | cat      | long      | short     |
    |----------|-----------|-----------|
    | Specific | 6.374502  | 16.89243  |
    | other    | 21.593625 | 55.139442 |
    the effect of short.habituation on cat in mouse is q =  0.03 with p-value = 0.86211, dof = 1
    looking at main.experimenter.gender effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 9.083665 | 14.183267 |
    | other    | 26.85259 | 49.880478 |
    the effect of main.experimenter.gender on cat in mouse is q =  1.42 with p-value = 0.23299, dof = 1
    looking at rodent.sex effect in cat in mouse
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 4.461942 | 18.897638 |
    | other    | 9.623797 | 67.016623 |
    the effect of rodent.sex on cat in mouse is q =  6.71 with p-value = 0.00958, dof = 1
    looking at head-plate effect in cat in mouse
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 2.549801 | 20.717131 |
    | other    | 6.135458 | 70.59761  |
    the effect of head-plate on cat in mouse is q =  2.12 with p-value = 0.1453, dof = 1
    looking at body.restrained effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 7.729084  | 15.537849 |
    | other    | 21.434263 | 55.298805 |
    the effect of body.restrained on cat in mouse is q =  2.78 with p-value = 0.09547, dof = 1
    looking at anesthesia.before.acquisition effect in cat in mouse
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 2.709163  | 20.557769 |
    | other    | 14.262948 | 62.47012  |
    the effect of anesthesia.before.acquisition on cat in mouse is q =  7.18 with p-value = 0.00736, dof = 1
    looking at MRI.field.strength effect in cat in mouse
    | cat      | 7.0      | 9.4       | 11.7     | 15.2     |
    |----------|----------|-----------|----------|----------|
    | Specific | 0.717131 | 18.964143 | 2.47012  | 1.115538 |
    | other    | 1.992032 | 66.454183 | 5.577689 | 2.709163 |
    the effect of MRI.field.strength on cat in mouse is q =  4.93 with p-value = 0.1768, dof = 3
    looking at fMRI.sequence effect in cat in mouse
    | cat      | GE-EPI    | SE-EPI    |
    |----------|-----------|-----------|
    | Specific | 7.968127  | 15.298805 |
    | other    | 24.621514 | 52.111554 |
    the effect of fMRI.sequence on cat in mouse is q =  0.38 with p-value = 0.53635, dof = 1
    #### NOW DOING rat ####
    summary of the data that we collected
    we processed 9 datasets
    totalling 468 runs
    from 161 animals
    the smallest dataset had 5 runs
    the largest dataset had 291 runs
    we could processed 235/468 runs.
    below is a summary of the data included per dataset
    | rodent.ds | total_run | total_animal | total_included | strain         |
    |-----------|-----------|--------------|----------------|----------------|
    | 2001      | 16        | 16           | 16             | Wistar         |
    | 2002a     | 7         | 7            | 7              | Sprague-Dawley |
    | 2002b     | 5         | 5            | 5              | Wistar         |
    | 2002c     | 5         | 5            | 5              | Wistar         |
    | 2003      | 97        | 8            | 97             | Sprague-Dawley |
    | 2004      | 19        | 5            | 19             | null           |
    | 2005      | 5         | 3            | 3              | null           |
    | 2006      | 23        | 23           | 23             | null           |
    | 4001      | 291       | 89           | 60             | null           |
    information about sex ratio
    the datasets contained 423 male runs and 22 female runs
    that corresponds to 4.94% females 
    information about animal handling
    2 datasets used headplates
    8 datasets used body restraining
    0 datasets used anesthesia before acquisition
    5 datasets were collected by men, 3 by women
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
    | 2006      | null      | null       | null       | null       | null            | null            |
    | 4001      | n         | y          | y          | m          | 7               | 330             |
    information about the scanner and sequence
    lowest field strength was 7T
    highest field strength was 7T
    | rodent.ds | field_strength | sequence | TE   |
    |-----------|----------------|----------|------|
    | 2001      | 7              | GE-EPI   | 17   |
    | 2002a     | 7              | GE-EPI   | 18   |
    | 2002b     | 7              | SE-EPI   | 45   |
    | 2002c     | 7              | SE-EPI   | 45   |
    | 2003      | 7              | GE-EPI   | 15   |
    | 2004      | null           | GE-EPI   | 15   |
    | 2005      | 7              | GE-EPI   | 17   |
    | 2006      | null           | null     | null |
    | 4001      | 7              | GE-EPI   | 15   |
    #### MOTION ANALYISIS ####
    mean fd across all rat datasets
    | fd.mean  |
    |----------|
    | 0.041561 |
    mean fd per dataset
    | rodent.ds | fd.mean  |
    |-----------|----------|
    | 2001      | 0.03643  |
    | 2002a     | 0.054875 |
    | 2002b     | 0.044527 |
    | 2002c     | 0.044339 |
    | 2003      | 0.04796  |
    | 2004      | 0.016877 |
    | 2005      | 0.022399 |
    | 2006      | null     |
    | 4001      | 0.039548 |
    t-test for head-plate mean.fd > no head-plate mean.fd
    t=0.88, p=0.37768, dof=199.39
    #### tSNR ANALYISIS ####
    tSNR across all rat datasets
    | s1.tsnr.l |
    |-----------|
    | 28.752543 |
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
    | 2005      | 2400.0       | 285.666667          | 179.666667          | 3.0                 |
    | 2006      | null         | null                | null                | null                |
    | 4001      | 600.0        | 103.85              | 79.716667           | 13.151515           |
    #### FC specifiticy analysis ####
     
    overall FC specificity for gsr1
    | s1.cat.gsr1  | count    |
    |--------------|----------|
    | Specific     | 0.348936 |
    | Spurious     | 0.314894 |
    | Non-specific | 0.187234 |
    | No           | 0.144681 |
    | null         | 0.004255 |
     
    overall FC specificity for gsr2
    | s1.cat.gsr2  | count    |
    |--------------|----------|
    | Specific     | 0.348936 |
    | Spurious     | 0.319149 |
    | Non-specific | 0.182979 |
    | No           | 0.144681 |
    | null         | 0.004255 |
     
    overall FC specificity for gsr3
    | s1.cat.gsr3  | count    |
    |--------------|----------|
    | Specific     | 0.370213 |
    | Spurious     | 0.331915 |
    | Non-specific | 0.178723 |
    | No           | 0.114894 |
    | null         | 0.004255 |
     
    overall FC specificity for wmcsf1
    | s1.cat.wmcsf1 | count    |
    |---------------|----------|
    | Spurious      | 0.33617  |
    | Non-specific  | 0.280851 |
    | Specific      | 0.26383  |
    | No            | 0.119149 |
     
    overall FC specificity for wmcsf2
    | s1.cat.wmcsf2 | count    |
    |---------------|----------|
    | Spurious      | 0.344681 |
    | Specific      | 0.285106 |
    | Non-specific  | 0.238298 |
    | No            | 0.12766  |
    | null          | 0.004255 |
     
    overall FC specificity for wmcsf3
    | s1.cat.wmcsf3 | count    |
    |---------------|----------|
    | Specific      | 0.319149 |
    | Spurious      | 0.319149 |
    | Non-specific  | 0.255319 |
    | No            | 0.102128 |
    | null          | 0.004255 |
     
    overall FC specificity for aCompCor1
    | s1.cat.aCompCor1 | count    |
    |------------------|----------|
    | Spurious         | 0.314894 |
    | Specific         | 0.302128 |
    | Non-specific     | 0.255319 |
    | No               | 0.123404 |
    | null             | 0.004255 |
     
    overall FC specificity for aCompCor2
    | s1.cat.aCompCor2 | count    |
    |------------------|----------|
    | Spurious         | 0.302128 |
    | Specific         | 0.297872 |
    | Non-specific     | 0.280851 |
    | No               | 0.114894 |
    | null             | 0.004255 |
     
    overall FC specificity for aCompCor3
    | s1.cat.aCompCor3 | count    |
    |------------------|----------|
    | Specific         | 0.323404 |
    | Spurious         | 0.276596 |
    | Non-specific     | 0.268085 |
    | No               | 0.12766  |
    | null             | 0.004255 |
    for gsr processing, heavy scrubbing lead to  82/234 Specific scans
    for gsr processing, light scrubbing lead to  87/234 Specific scans
    for aCompCor processing, heavy scrubbing lead to  71/234 Specific scans
    for aCompCor processing, light scrubbing lead to  76/234 Specific scans
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
    | 2005      | 0.333333  | 0.333333  | 0.0       | … | 0.0       | 0.0       | 0.0       | 0.0      |
    | 2006      | 0.0       | 0.0       | 0.0       | … | 0.0       | 0.0       | 0.0       | 0.0      |
    | 4001      | 0.35      | 0.366667  | 0.383333  | … | 0.35      | 0.366667  | 0.4       | 0.4      |
    it seems that less scrubbing lead to better outcomes
    for gsr processing, 3/9 dataset performed better with less scrubbing
    for aCompCor processing, 4/9 dataset performed better with less scrubbing
    for comparing gsr to aCompCor processing, 2/9 dataset performed better with aCompCor compared to gsr
    for comparing gsr to aCompCor processing, 4/9 dataset performed worst with aCompCor compared to gsr
    #### plot fc categories per denoising method ####
    doing statistical analysis
    first with gsr3 processed data
    looking at fd.mean effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 8.056872  | 12.796209 | 10.42654  | 9.952607  |
    | other    | 17.061611 | 12.322275 | 14.218009 | 15.165877 |
    the effect of fd.mean on cat in rat is q =  3.98 with p-value = 0.26383, dof = 3
    looking at s1.gsrcov.l.gsr3 effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 10.42654  | 7.109005  | 10.900474 | 12.796209 |
    | other    | 14.691943 | 18.009479 | 13.744076 | 12.322275 |
    the effect of s1.gsrcov.l.gsr3 on cat in rat is q =  5.91 with p-value = 0.11585, dof = 3
    looking at s1.tsnr.l effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 7.582938  | 10.42654  | 12.322275 | 10.900474 |
    | other    | 17.535545 | 14.691943 | 12.322275 | 14.218009 |
    the effect of s1.tsnr.l on cat in rat is q =  4.42 with p-value = 0.21941, dof = 3
    looking at habituation.min effect in cat in rat
    | cat      | low       | intermediate |
    |----------|-----------|--------------|
    | Specific | 17.924528 | 23.113208    |
    | other    | 25.943396 | 33.018868    |
    the effect of habituation.min on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
    looking at habituation.days effect in cat in rat
    | cat      | low       | intermediate |
    |----------|-----------|--------------|
    | Specific | 12.735849 | 28.301887    |
    | other    | 21.226415 | 37.735849    |
    the effect of habituation.days on cat in rat is q =  0.36 with p-value = 0.54614, dof = 1
    looking at short.habituation effect in cat in rat
    | cat      | long      | short    |
    |----------|-----------|----------|
    | Specific | 35.319149 | 1.702128 |
    | other    | 59.574468 | 3.404255 |
    the effect of short.habituation on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
    looking at main.experimenter.gender effect in cat in rat
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 26.886792 | 14.150943 |
    | other    | 29.716981 | 29.245283 |
    the effect of main.experimenter.gender on cat in rat is q =  4.18 with p-value = 0.04098, dof = 1
    looking at rodent.sex effect in cat in rat
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 3.773585 | 37.264151 |
    | other    | 5.660377 | 53.301887 |
    the effect of rodent.sex on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
    looking at head-plate effect in cat in rat
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 17.924528 | 23.113208 |
    | other    | 27.358491 | 31.603774 |
    the effect of head-plate on cat in rat is q =  0.06 with p-value = 0.80152, dof = 1
    looking at body.restrained effect in cat in rat
    | cat      | y         |
    |----------|-----------|
    | Specific | 41.037736 |
    | other    | 58.962264 |
    the effect of body.restrained on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at anesthesia.before.acquisition effect in cat in rat
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 1.886792 | 39.150943 |
    | other    | 7.075472 | 51.886792 |
    the effect of anesthesia.before.acquisition on cat in rat is q =  2.6 with p-value = 0.10703, dof = 1
    looking at MRI.field.strength effect in cat in rat
    | cat      | 7         |
    |----------|-----------|
    | Specific | 43.005181 |
    | other    | 56.994819 |
    the effect of MRI.field.strength on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at fMRI.sequence effect in cat in rat
    | cat      | GE-EPI    | SE-EPI   |
    |----------|-----------|----------|
    | Specific | 39.622642 | 1.415094 |
    | other    | 55.660377 | 3.301887 |
    the effect of fMRI.sequence on cat in rat is q =  0.16 with p-value = 0.6909, dof = 1
    doing statistical analysis
    first with aCompCor3 processed data
    looking at fd.mean effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 7.582938  | 11.374408 | 8.056872  | 9.004739  |
    | other    | 17.535545 | 13.744076 | 16.587678 | 16.113744 |
    the effect of fd.mean on cat in rat is q =  3.01 with p-value = 0.39072, dof = 3
    looking at s1.gsrcov.l.aCompCor3 effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 9.004739  | 8.530806  | 10.42654  | 8.056872  |
    | other    | 16.113744 | 16.587678 | 14.218009 | 17.061611 |
    the effect of s1.gsrcov.l.aCompCor3 on cat in rat is q =  1.35 with p-value = 0.71778, dof = 3
    looking at s1.tsnr.l effect in cat in rat
    | cat      | lowest    | low       | high      | highest   |
    |----------|-----------|-----------|-----------|-----------|
    | Specific | 7.109005  | 10.900474 | 8.056872  | 9.952607  |
    | other    | 18.009479 | 14.218009 | 16.587678 | 15.165877 |
    the effect of s1.tsnr.l on cat in rat is q =  3.17 with p-value = 0.3662, dof = 3
    looking at habituation.min effect in cat in rat
    | cat      | low       | intermediate |
    |----------|-----------|--------------|
    | Specific | 17.924528 | 17.924528    |
    | other    | 25.943396 | 38.207547    |
    the effect of habituation.min on cat in rat is q =  1.44 with p-value = 0.22986, dof = 1
    looking at habituation.days effect in cat in rat
    | cat      | low       | intermediate |
    |----------|-----------|--------------|
    | Specific | 12.264151 | 23.584906    |
    | other    | 21.698113 | 42.45283     |
    the effect of habituation.days on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
    looking at short.habituation effect in cat in rat
    | cat      | long      | short    |
    |----------|-----------|----------|
    | Specific | 31.489362 | 0.851064 |
    | other    | 63.404255 | 4.255319 |
    the effect of short.habituation on cat in rat is q =  0.77 with p-value = 0.38169, dof = 1
    looking at main.experimenter.gender effect in cat in rat
    | cat      | f         | m         |
    |----------|-----------|-----------|
    | Specific | 20.754717 | 15.09434  |
    | other    | 35.849057 | 28.301887 |
    the effect of main.experimenter.gender on cat in rat is q =  0.02 with p-value = 0.88943, dof = 1
    looking at rodent.sex effect in cat in rat
    | cat      | f        | m         |
    |----------|----------|-----------|
    | Specific | 4.245283 | 31.603774 |
    | other    | 5.188679 | 58.962264 |
    the effect of rodent.sex on cat in rat is q =  0.42 with p-value = 0.51457, dof = 1
    looking at head-plate effect in cat in rat
    | cat      | n         | y         |
    |----------|-----------|-----------|
    | Specific | 17.924528 | 17.924528 |
    | other    | 27.358491 | 36.792453 |
    the effect of head-plate on cat in rat is q =  0.79 with p-value = 0.37477, dof = 1
    looking at body.restrained effect in cat in rat
    | cat      | y         |
    |----------|-----------|
    | Specific | 35.849057 |
    | other    | 64.150943 |
    the effect of body.restrained on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at anesthesia.before.acquisition effect in cat in rat
    | cat      | n        | y         |
    |----------|----------|-----------|
    | Specific | 1.886792 | 33.962264 |
    | other    | 7.075472 | 57.075472 |
    the effect of anesthesia.before.acquisition on cat in rat is q =  1.34 with p-value = 0.24651, dof = 1
    looking at MRI.field.strength effect in cat in rat
    | cat      | 7         |
    |----------|-----------|
    | Specific | 37.305699 |
    | other    | 62.694301 |
    the effect of MRI.field.strength on cat in rat is q =  0.0 with p-value = 1.0, dof = 0
    looking at fMRI.sequence effect in cat in rat
    | cat      | GE-EPI    | SE-EPI   |
    |----------|-----------|----------|
    | Specific | 33.962264 | 1.886792 |
    | other    | 61.320755 | 2.830189 |
    the effect of fMRI.sequence on cat in rat is q =  0.0 with p-value = 1.0, dof = 1
