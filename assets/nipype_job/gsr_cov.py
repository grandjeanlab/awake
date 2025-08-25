from os.path import join
import sys

from nipype import Node, Workflow, IdentityInterface, Function
from nipype.interfaces.io import DataFinder, DataSink
from nipype.interfaces.nilearn import SignalExtraction


def keep_filename(out_paths):
    return out_paths.split('/')[-1].split('.')[0] + '.tsv'

rodent = sys.argv[1] if len(sys.argv) > 1 else 'mouse'
analysis = sys.argv[2] if len(sys.argv) > 2 else 'gsr1'

data_dir = '/project/4180000.36/awake/complete_output_'+rodent
wf_dir = '/project/4180000.36/awake/analysis_'+rodent+'/workflow'
ds_dir = '/project/4180000.36/awake/analysis_'+rodent+'/datasink'
label_img = '/home/traaffneu/joagra/code/awake/assets/template/'+rodent+'/s1bf.nii.gz'


#get the list of nifti images in the folder of interest
df_node = Node(DataFinder(), name='df_node')
df_node.inputs.root_paths = join(data_dir, analysis, 'data_diagnosis_datasink/GS_cov_nii/')
df_node.inputs.match_regex = r'^(?P<basedir>.+)/(?P<filename>[^/]+)\.nii\.gz$'
df_list = df_node.run('MultiProc').outputs

#iterate over list of nifti images 
infosource = Node(IdentityInterface(fields=['out_paths']),name="infosource")
infosource.iterables = [('out_paths', df_list.out_paths)]

#singal extrraction node
se_node = Node(SignalExtraction(), name="se_node")
se_node.inputs.label_files = label_img  
se_node.inputs.class_labels = ['s1_left', 's1_right']

#stupid filename node because i could not make it work otherwise
filename_node = Node(Function(input_names=["out_paths"],
                       output_names=["filename"],
                       function=keep_filename),
              name='filename_node')

#datasink node
ds_node = Node(DataSink(), name="ds_node")
ds_node.inputs.base_directory = join(ds_dir, 'gsr_cov', analysis)
ds_node.inputs.parameterization = False

#define the workflow
wf = Workflow(name='workflow_gsr_cov'+analysis)
wf.base_dir = wf_dir

wf.connect(infosource, 'out_paths', se_node, 'in_file')
wf.connect(infosource, 'out_paths', filename_node, 'out_paths')
wf.connect(filename_node, 'filename', se_node, 'out_file')
wf.connect(se_node, 'out_file', ds_node, 'out_file')
wf.run(plugin='MultiProc')

