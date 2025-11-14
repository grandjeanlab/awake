# figure plots
Joanes Grandjean
2025-11-10

``` python
%autoindent 
```

``` python
import polars as pl
#import plotting modules
import matplotlib.pyplot as plt
import seaborn as sns
from met_brewer import met_brew, is_colorblind_friendly
import matplotlib
#set matplotlib rendering to TkAgg
matplotlib.use('TkAgg')


def makeviolinplot(df_to_plot):
  ax = sns.violinplot(x=df_to_plot['value'], y=df_to_plot['cont_variable'], hue=df_to_plot['value'],  inner=None)
  ax = sns.stripplot(x=df_to_plot['value'], y=df_to_plot['cont_variable'],  dodge=True, color='black', alpha=0.3)
  return ax

def makestackplot(df_melted):
  ax = sns.displot(x=df_melted['variable'], hue=df_melted['value'], multiple='stack', hue_order=cat_index, stat='probability', legend = None)
  ax.ax.invert_yaxis()
  ax.ax.tick_params(axis='x', rotation=90)
  return ax

def makespecificityplot(x,y, hue=None):
  ax = sns.jointplot(x=x, y=y, hue=hue, kind="kde")
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

pl.Config(
    tbl_formatting="MARKDOWN",
    tbl_hide_column_data_types=True,
    tbl_hide_dataframe_shape=True,
    tbl_rows=100,
)


analysis_list = [ "gsr1", "gsr2", "gsr3","wmcsf1", "wmcsf2", "wmcsf3", "aCompCor1", "aCompCor2", "aCompCor3" ]  

cat_index = ["Specific", "Non-specific", "Spurious", "No"]

#set color theme for all the plots
plt.rcParams['font.size'] = 6
plt.rcParams['axes.spines.right'] = False
plt.rcParams['axes.spines.top'] = False


met_colors_name = "Hiroshige"

met_colors = met_brew(name=met_colors_name)
print(is_colorblind_friendly(met_colors_name))

sns.set_theme(style="white", palette=met_colors)
sns_saveparms = "bbox_inches='tight', dppi=600"
```

    Palette 'Hiroshige' has '10' discrete colors
    Palette 'Hiroshige' is colorblind friendly.
    True

``` python
def figure1(df, df_summary, rodent, save_plot):
  #figure 1ab
  if rodent == 'mouse':
    pannel_label = 'a'
  else:
    pannel_label = 'b'
  plt.figure(figsize=(10,10))
  ax = sns.barplot(x=df_summary["rodent.ds"], y=df_summary["habituation.min"], hue=df_summary["headplate"], dodge=False)
  for i in df_summary.with_row_index().iter_rows():
    ax.text(i[0], i[19], i[18], color='black', ha='center', va='bottom')
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  ax.set(xlabel='Dataset ID', ylabel='Habituation time (min)')
  ax.legend_.set_title('Head-plate')
  if save_plot:
    ax.figure.savefig("../figure/Figure1"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()
  #figure 1ce
  if rodent == 'mouse':
    pannel_label = 'c'
  else:
    pannel_label = 'e'
  plt.figure(figsize=(5,10))
  ax = sns.stripplot(x=df["rodent.ds"], y=df["fd.mean"], hue=df['head-plate'])
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  ax.set(xlabel='Dataset ID', ylabel='Framewise displacement (mm)')
  ax.legend_.set_title('Head-plate')
  if save_plot:
    ax.figure.savefig("../figure/Figure1"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()
  #figure 1df
  if rodent == 'mouse':
    pannel_label = 'd'
  else:
    pannel_label = 'f'
  df_to_plot = df.select('head-plate','fd.mean').rename({'head-plate':'value','fd.mean':'cont_variable'})
  plt.figure(figsize=(5,10))
  ax = makeviolinplot(df_to_plot)
  ax.set(xlabel='Head-plate', ylabel='Framewise displacement (mm)')
  if save_plot:
    ax.figure.savefig("../figure/Figure1"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()
  #figure 1gi
  if rodent == 'mouse':
    pannel_label = 'g'
  else:
    pannel_label = 'i'
  plt.figure(figsize=(5,10))
  ax = sns.scatterplot(x=df["habituation.min"], y=df["fd.mean"], hue=df['head-plate'])
  ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
  ax.set(xlabel='Habituation (min)', ylabel='Framewise displacement (mm)')
  ax.legend_.set_title('Head-plate')
  if save_plot:
    ax.figure.savefig("../figure/Figure1"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()
  #figure 1hj
  if rodent == 'mouse':
    pannel_label = 'h'
  else:
    pannel_label = 'j'
  df_to_plot = df.select('short.habituation','fd.mean').rename({'short.habituation':'value','fd.mean':'cont_variable'})
  plt.figure(figsize=(5,10))
  ax = makeviolinplot(df_to_plot)
  ax.set(xlabel='Habituation', ylabel='Framewise displacement (mm)')
  if save_plot:
    ax.figure.savefig("../figure/Figure1"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()

def figure2(df, rodent, save_plot):
#figure 2ab
  if rodent == 'mouse':
    pannel_label = 'a'
  else:
    pannel_label = 'b'
  df_melted = df.melt(id_vars="rodent.ds", value_vars=['s1.cat.' + x for x in analysis_list])
  plt.figure(figsize=(10,5))
  ax = makestackplot(df_melted)
  ax.set(xlabel='Denoising model', ylabel='Runs (%)')
  if save_plot:
    ax.figure.savefig("../figure/Figure2"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()
#figure 2cd
  if rodent == 'mouse':
    analysis = 'wmcsf3'
    pannel_label = 'c'
  else:
    analysis = 'gsr3'
    pannel_label = 'd'
  plt.figure(figsize=(10,5))
  ax = makespecificityplot(df['s1.specific.'+analysis], df['s1.unspecific.'+analysis])
  if save_plot:
    ax.figure.savefig("../figure/Figure2"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
    plt.show()

def figure3(df, rodent, save_plot):
#figure3ab
  if rodent == 'mouse':
      analysis = 'wmcsf3'
      pannel_label = 'a'
  else:
      analysis = 'gsr3'
      pannel_label = 'b'
  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
  cat2 = 'fd.mean'
  df_to_plot = split_continuous(df, cat1, cat2)
  df_to_plot = df_to_plot.rename({cat1:'value', 'quartiles':'variable', cat2:'cont_variable'})
  plt.figure(figsize=(10,10))
  ax = makestackplot(df_to_plot)
  ax.set(xlabel='Framewise displacement (bin)', ylabel='Runs (%)')
  if save_plot:
      ax.figure.savefig("../figure/Figure3"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
      plt.show()
#figure3ce
  if rodent == 'mouse':
      analysis = 'wmcsf3'
      pannel_label = 'c'
  else:
      analysis = 'gsr3'
      pannel_label = 'e'
  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
  cat2='short.habituation'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(5,10))
  ax = sns.barplot(x=df_to_plot[cat2],y=df_to_plot['len'],hue=df_to_plot['cat'])
  ax.set(xlabel='Habituation', ylabel='Runs (%)')
  if save_plot:
      ax.figure.savefig("../figure/Figure3"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
      plt.show()
#figure3df
  if rodent == 'mouse':
      analysis = 'wmcsf3'
      pannel_label = 'd'
  else:
      analysis = 'gsr3'
      pannel_label = 'f'
  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
  cat2='head-plate'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(5,10))
  ax = sns.barplot(x=df_to_plot[cat2],y=df_to_plot['len'],hue=df_to_plot['cat'])
  ax.set(xlabel='Head-plate', ylabel='Runs (%)')
  if save_plot:
      ax.figure.savefig("../figure/Figure3"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
      plt.show()
#figure3gi
  if rodent == 'mouse':
      analysis = 'wmcsf3'
      pannel_label = 'g'
  else:
      analysis = 'gsr3'
      pannel_label = 'i'
  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
  cat2='main.experimenter.gender'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(5,10))
  ax = sns.barplot(x=df_to_plot[cat2],y=df_to_plot['len'],hue=df_to_plot['cat'])
  ax.set(xlabel='Experimenter gender', ylabel='Runs (%)')
  if save_plot:
      ax.figure.savefig("../figure/Figure3"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
      plt.show()
#figure3hj
  if rodent == 'mouse':
      analysis = 'wmcsf3'
      pannel_label = 'h'
  else:
      analysis = 'gsr3'
      pannel_label = 'j'
  cat1='s1.cat.'+analysis
  df = df.with_columns(pl.when(pl.col(cat1)=='Specific').then(pl.lit('Specific')).otherwise(pl.lit('other')).alias('cat'))
  cat2='fMRI.sequence'
  df_to_plot = df.select('cat',cat2).drop_nulls().group_by(['cat',cat2]).agg(pl.len()).sort(by=cat2).with_columns((pl.col('len') / pl.col('len').sum().over('cat')).alias('rel_count')).sort(by='cat')
  plt.figure(figsize=(5,10))
  ax = sns.barplot(x=df_to_plot[cat2],y=df_to_plot['len'],hue=df_to_plot['cat'])
  ax.set(xlabel='fMRI sequence', ylabel='Runs (%)')
  if save_plot:
      ax.figure.savefig("../figure/Figure3"+pannel_label+".svg", bbox_inches='tight',dpi=600)
  else:
      plt.show()
```

``` python
save_plot = True

rodent_list = ['mouse', 'rat']
for rodent in rodent_list:
  print("#### NOW DOING " + rodent + " ####")
  df = pl.read_csv("../assets/tables/"+rodent+"_metadata_processed.tsv", separator="\t").sort(by='rodent.ds')
  df_summary=pl.read_csv("../assets/tables/"+rodent+"_summary_processed.tsv", separator="\t").sort(by='rodent.ds')
  figure1(df, df_summary, rodent, save_plot)
  figure2(df, rodent, save_plot)
  figure3(df, rodent, save_plot)
```

    #### NOW DOING mouse ####

    #### NOW DOING rat ####
