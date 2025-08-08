# summary of the datasets
Joanes Grandjean
Invalid Date

We accumulated datasets from several laboratories. Here, we carry out a
descriptive summary of the dataset and its processing.

``` python
#only run this in interactive mode to turn off the autoindent in ipython
%autoindent 
```

``` python
import pandas as pd
def get_rodent_summary(df):
  #get a summary of the total number of runs, excluded runs and number of animals
  df_summary = df.groupby("rodent.ds").size().reset_index(name='total_run')
  df_summary["total_included"] = df.groupby("rodent.ds")["exclude"].apply(lambda x: x.isna().sum()).values
  df_summary["included_percentage"] = (df_summary["total_included"] / df_summary["total_run"]) * 100
  df_summary["total_mouse"] = df.groupby("rodent.ds")["rodent.sub"].nunique().values
  #get a summary of the mice used 
  df_summary["strain"] = df.groupby("rodent.ds")["rodent.strain"].unique().values
  df_summary["male"] =  df.groupby("rodent.ds")["rodent.sex"].apply(lambda x: (x == 'm').sum()).values
  df_summary["female"] = df.groupby("rodent.ds")["rodent.sex"].apply(lambda x: (x == 'f').sum()).values
  df_summary["unknownsex"] = df.groupby("rodent.ds")["rodent.sex"].apply(lambda x: x.isna().sum()).values
  #get a summary of the animal preparation
  df_summary["headplate"] = df.groupby("rodent.ds")["head-plate"].unique().values
  df_summary["restrained"] = df.groupby("rodent.ds")["body.restrained"].unique().values
  df_summary["anesthesia"] = df.groupby("rodent.ds")["anesthesia.before.acquisition"].unique().values
  df_summary["exp.gender"] = df.groupby("rodent.ds")["main.experimenter.gender"].unique().values
  #get a summary of the MRI parameters
  df_summary["field_strength"] = df.groupby("rodent.ds")["MRI.field.strength"].unique().values
  df_summary["sequence"] = df.groupby("rodent.ds")["fMRI.sequence"].unique().values
  df_summary["TR"] = df.groupby("rodent.ds")["MRI.TR"].unique().values
  df_summary["TE"] = df.groupby("rodent.ds")["MRI.TE"].unique().values
  # get aidaqc summary
  df_summary["aidaqc.tsnr.mean"] = df.groupby("rodent.ds")["aidaqc.tsnr"].mean().values
  df_summary["aidaqc.tsnr.std"] = df.groupby("rodent.ds")["aidaqc.tsnr"].std().values
  df_summary["aidaqc.dist.mean"] = df.groupby("rodent.ds")["aidaqc.dist"].mean().values
  df_summary["aidaqc.dist.std"] = df.groupby("rodent.ds")["aidaqc.dist"].std().values
  return df_summary


for rodent in [ "mouse", "rat" ]:
  df = pd.read_csv("../assets/tables/"+rodent+"_metadata.tsv", sep="\t")
  df["scan"] = "sub-0" + df["rodent.sub"].astype("str") + "_ses-" + df["rodent.session"].astype("str") + "_run-" + df["rodent.run"].astype("str")
  aidaqc = pd.read_csv("../assets/tables/"+rodent+"_caculated_features_func.csv", sep=",")
  aidaqc = aidaqc.rename(columns={"tSNR (Averaged Brain ROI)": "aidaqc.tsnr", "Displacement factor (std of Mutual information)": "aidaqc.dist"})
  aidaqc["scan"] = aidaqc["FileAddress"].apply(lambda x : x.split("func/")[1].split("_task")[0])
  df = df.set_index("scan").join(aidaqc.set_index("scan"))
  df_summary = get_rodent_summary(df)
  #saving the table for later use
  df_summary.to_csv("../assets/tables/"+rodent+"_summary.tsv", sep="\t", index=False)
  print("summary data for "+rodent)
  print(df_summary[["total_run", "total_included", "included_percentage", "total_mouse", "strain"]])
```

    summary data for mouse
        total_run  total_included  included_percentage  total_mouse         strain
    0           4               4           100.000000            4      [C57BL/6]
    1           8               8           100.000000            8      [C57BL/6]
    2          12              12           100.000000           12      [C57BL/6]
    3          21              21           100.000000            5      [C57BL/6]
    4          13              13           100.000000            4      [C57BL/6]
    5          10               0             0.000000           10      [C57BL/6]
    6          20              20           100.000000           10      [C57BL/6]
    7          14              14           100.000000            7      [C57BL/6]
    8          36              34            94.444444            9      [C57BL/6]
    9         107              34            31.775701            5  [129S2/SvPas]
    10         10               0             0.000000           10      [C57BL/6]
    11         51              48            94.117647           17      [C57BL/6]
    12        112             112           100.000000            6      [C57BL/6]
    13         10              10           100.000000           10      [C57BL/6]
    14        479             414            86.430063           19   [B6129PF/J1]
    summary data for rat
       total_run  total_included  included_percentage  total_mouse  \
    0         16              16           100.000000           16   
    1          7               7           100.000000            7   
    2         10              10           100.000000           10   
    3         97              97           100.000000            8   
    4        291              60            20.618557           89   

                 strain  
    0          [Wistar]  
    1  [Sprague-Dawley]  
    2          [Wistar]  
    3  [Sprague-Dawley]  
    4             [nan]  
