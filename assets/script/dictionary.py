import polars as pl
import os
import glob
import warnings
from sklearn.exceptions import ConvergenceWarning
from nilearn.decomposition import CanICA, DictLearning


#get the TMPDIR environment variable for caching
tmpdir=os.getenv("TMPDIR")

rodent = 'mouse'
codedir = '/home/traaffneu/joagra/code/awake'
datadir = '/project/4180000.36/awake/complete_output_mouse/'
preprocess = 'gsr3'

mask = os.path.join(codedir, 'assets', 'template',rodent,'mask04.nii.gz')

df = pl.read_csv(os.path.join(codedir,'assets','tables',rodent+"_metadata.tsv"), separator="\t", ignore_errors=True)

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


#list all files that have the pattern 'bold_cleaned.nii.gz' in joined path of datadir and preprocess
file_list = glob.glob(os.path.join(datadir, preprocess,'confound_correction_datasink/cleaned_timeseries/_split_name*/*'))

#list in df the column 'scan' where 'exclude' is empty
scan_fliter = df.filter(pl.col("exclude").is_null()).select("scan")

#filter file_list to only include files that have the scan name in scan_fliter
file_list = [f for f in file_list if any(scan in f for scan in scan_fliter["scan"].to_list())]

print('starting dictionary learning with {} files'.format(len(file_list)))

dict_learning = DictLearning(
    n_components=20,
    memory=tmpdir,
    memory_level=1,
    verbose=1,
    mask=mask,
    random_state=0,
    n_epochs=3,
    smoothing_fwhm=None,
    standardize="zscore_sample",
    n_jobs=-1,
)

dict_learning.fit(file_list)

print('outputing dictionary to file')
dict_learning.components_img_.to_filename(os.path.join(codedir,'assets','template',rodent,'dict_'+preprocess+'.nii.gz'))


print('starting canica with {} files'.format(len(file_list)))

canica = CanICA(
        n_components=20,
        memory=tmpdir,
        memory_level=1,
        verbose=1,
        mask=mask,
        smoothing_fwhm=None,
        standardize="zscore_sample",
        n_jobs=-1,
)

canica.fit(file_list)
print('outputing canica to file')
canica.components_img_.to_filename(os.path.join(codedir,'assets','template',rodent,'canica_'+preprocess+'.nii.gz'))

print('done')
