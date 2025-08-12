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
import polars as pl

def get_rodent_summary(df):
    # Get a summary of the total number of runs, excluded runs and number of animals
    df_summary = df.group_by("rodent.ds").agg([
        pl.count().alias("total_run"),
        (pl.col("exclude").is_null().sum()).alias("total_included"),
        pl.col("rodent.sub").n_unique().alias("total_mouse"),
        pl.col("rodent.strain").unique().first().alias("strain"),
        (pl.col("rodent.sex") == 'm').sum().alias("male"),
        (pl.col("rodent.sex") == 'f').sum().alias("female"),
        pl.col("rodent.sex").is_null().sum().alias("unknownsex"),
        pl.col("head-plate").unique().first().alias("headplate"),
        pl.col("body.restrained").unique().first().alias("restrained"),
        pl.col("anesthesia.before.acquisition").unique().first().alias("anesthesia"),
        pl.col("main.experimenter.gender").unique().first().alias("exp.gender"),
        pl.col("MRI.field.strength").unique().first().alias("field_strength"),
        pl.col("fMRI.sequence").unique().first().alias("sequence"),
        pl.col("MRI.TE").unique().first().alias("TE"),
    ])
    df_summary = df_summary.with_columns([
        (pl.col("total_included") / pl.col("total_run") * 100).alias("included_percentage")
    ])
    return df_summary


    

for rodent in [ "mouse", "rat" ]:
# Read the CSV file
  df = pl.read_csv("../assets/tables/"+rodent+"_metadata.tsv", separator="\t", ignore_errors=True)
# Create a new column "scan"
  df = df.with_columns([
      pl.concat_str([
          pl.lit("sub-0"),
          pl.col("rodent.sub").cast(pl.Utf8),
          pl.lit("_ses-"),
          pl.col("rodent.session").cast(pl.Utf8),
          pl.lit("_run-"),
          pl.col("rodent.run").cast(pl.Utf8)
      ]).alias("scan")
  ])
# Get the rodent summary
  df_summary = get_rodent_summary(df)
# Save the table for later use
  df_summary.write_csv("../assets/tables/"+rodent+"_summary.tsv", separator="\t")
  print("now doing "+rodent)
  print(df_summary["rodent.ds", "total_run", "total_included", "total_mouse", "male", "female", "exp.gender"])
```

    now doing mouse
    shape: (18, 7)
    ┌───────────┬───────────┬────────────────┬─────────────┬──────┬────────┬────────────┐
    │ rodent.ds ┆ total_run ┆ total_included ┆ total_mouse ┆ male ┆ female ┆ exp.gender │
    │ ---       ┆ ---       ┆ ---            ┆ ---         ┆ ---  ┆ ---    ┆ ---        │
    │ str       ┆ u32       ┆ u32            ┆ u32         ┆ u32  ┆ u32    ┆ str        │
    ╞═══════════╪═══════════╪════════════════╪═════════════╪══════╪════════╪════════════╡
    │ 3002      ┆ 10        ┆ 10             ┆ 10          ┆ 10   ┆ 0      ┆ null       │
    │ 1005a     ┆ 13        ┆ 13             ┆ 4           ┆ 13   ┆ 0      ┆ f          │
    │ 1012      ┆ 26        ┆ 26             ┆ 26          ┆ 14   ┆ 12     ┆ f          │
    │ 3001      ┆ 112       ┆ 112            ┆ 6           ┆ 0    ┆ 0      ┆ null       │
    │ 1005b     ┆ 10        ┆ 0              ┆ 10          ┆ 10   ┆ 0      ┆ f          │
    │ …         ┆ …         ┆ …              ┆ …           ┆ …    ┆ …      ┆ …          │
    │ 1010      ┆ 10        ┆ 0              ┆ 10          ┆ 10   ┆ 0      ┆ null       │
    │ 1003      ┆ 12        ┆ 12             ┆ 12          ┆ 7    ┆ 5      ┆ m          │
    │ 1001      ┆ 4         ┆ 4              ┆ 4           ┆ 4    ┆ 0      ┆ m          │
    │ 3004      ┆ 54        ┆ 54             ┆ 9           ┆ 0    ┆ 54     ┆ f          │
    │ 1002      ┆ 8         ┆ 8              ┆ 8           ┆ 8    ┆ 0      ┆ m          │
    └───────────┴───────────┴────────────────┴─────────────┴──────┴────────┴────────────┘
    now doing rat
    shape: (7, 7)
    ┌───────────┬───────────┬────────────────┬─────────────┬──────┬────────┬────────────┐
    │ rodent.ds ┆ total_run ┆ total_included ┆ total_mouse ┆ male ┆ female ┆ exp.gender │
    │ ---       ┆ ---       ┆ ---            ┆ ---         ┆ ---  ┆ ---    ┆ ---        │
    │ str       ┆ u32       ┆ u32            ┆ u32         ┆ u32  ┆ u32    ┆ str        │
    ╞═══════════╪═══════════╪════════════════╪═════════════╪══════╪════════╪════════════╡
    │ 2002a     ┆ 7         ┆ 7              ┆ 7           ┆ 7    ┆ 0      ┆ f          │
    │ 2001      ┆ 16        ┆ 16             ┆ 16          ┆ 8    ┆ 8      ┆ f          │
    │ 2005      ┆ 5         ┆ 5              ┆ 3           ┆ 0    ┆ 5      ┆ m          │
    │ 2004      ┆ 19        ┆ 19             ┆ 5           ┆ 10   ┆ 9      ┆ m          │
    │ 2002b     ┆ 10        ┆ 10             ┆ 10          ┆ 10   ┆ 0      ┆ m          │
    │ 2003      ┆ 97        ┆ 97             ┆ 8           ┆ 97   ┆ 0      ┆ f          │
    │ 4001      ┆ 291       ┆ 60             ┆ 89          ┆ 291  ┆ 0      ┆ null       │
    └───────────┴───────────┴────────────────┴─────────────┴──────┴────────┴────────────┘
