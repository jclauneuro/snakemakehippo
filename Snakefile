from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'


#load participants.tsv file, and strip off sub- from participant_id column
df = pd.read_table(config['participants_tsv'])
subjects = df.participant_id.to_list() 
subjects = [ s.strip('sub-') for s in subjects ]


#get rois and hemis from config
rois = config['rois']
hemis = config['hemis']



wildcard_constraints:
    subject="[a-zA-Z0-9]+"



rule all:
    input:
        final = expand('working/sub-{template}_{roi}_{hemi}_fuzzy.nii.gz',template=config['template'],roi=rois,hemi=hemis)
    group: 'average_roi'

rule import_template_roi:
    input: join(config['template_seg_dir'],config['template_seg_nii'])
    output: 'working/template_masks/sub-{template}_hemi-{hemi}_desc-{roi}_mask.nii.gz'
    log: 'logs/import_template_roi/{template}_{roi}_{hemi}.log'
    group: 'template_rois_to_subject'
    shell: 'cp -v {input} {output} &> {log}'

rule transform_to_subject:
    input: 
        roi = rules.import_template_roi.output,
        affine = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(**wildcards))),
        invwarp = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_invwarp_nii'].format(**wildcards))),
        ref = join(config['subject_dir'],'sub-{subject}/anat/',config['subject_T1_ref_nii'])
    output: 'working/sub-{subject}/masks/sub-{subject}_space-{template}_hemi-{hemi}_desc-{roi}_mask.nii.gz'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{roi}_{hemi}.log'
    group: 'template_rois_to_subject'
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.roi} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'

# rule transform from subject back to template!
rule transform_subject_back_to_template:
    input:
        roi = rules.transform_to_subject.output,
        affine = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(**wildcards))),
        warp = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_warp_nii'].format(**wildcards))),
        ref = join(config['ants_template_dir'],config['ants_ref_nii']), 
    output: 'working/sub-{subject}/masks/sub-{template}_space-sub-{subject}_hemi-{hemi}_desc-{roi}_mask.nii.gz',
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_from_subject_back_to_template/{template}_sub-{subject}_{roi}_{hemi}.log'
    group: 'template_rois_to_template'
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.roi} -o {output} -r {input.ref} -t {input.warp} -t [{input.affine},0] &> {log}'

# final rule is just to average all the subjects together to create probabilistic atlas for ROI
rule average_roi:
    input: rois = expand('working/sub-{subject}/masks/sub-{template}_space-sub-{subject}_hemi-{hemi}_desc-{roi}_mask.nii.gz',subject=subjects,allow_missing=True)
    output: 'working/sub-{template}_{roi}_{hemi}_fuzzy.nii.gz'
    log: 'logs/average_roi/{template}_{roi}_{hemi}_average_roi.log'
    group: 'average_roi'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    shell:
        'AverageImages 3 {output} 0 {input.rois}'

