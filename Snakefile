from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'

prepdwi_dir = config['prepdwi_dir']
seg_dir = config['seg_dir']
bids_dir = config['bids_dir']
singularity_dir = config['singularity_dir']

seed_res = config['probtrack']['seed_res']
nsamples = config['probtrack']['nsamples']

participants_tsv = join(bids_dir,'participants.tsv')
subjects_table = pd.read_table(participants_tsv)

#get list of subjects based on seg_dir
subjects = os.listdir(seg_dir)
subjects = [ s.strip('sub-') for s in subjects ]
subjects = sorted(subjects)

template = 'snsx32'
ants_template_dir = '/project/6007967/jclau/snsx32/templates/snsx32_v0.2/snsx32_v0.2_i09'
labels = ['ZI']
hemis = ['L','R']


target_csv = 'cfg/HarvardOxford_noremap.csv'
ntargets = 49 #hardcoded right now..
tmaps = ['{tmap:04}'.format(tmap=t) for t in range(ntargets) ]
print(tmaps)

#just test with one subject:
#subjects = subjects[0]

# note: diffparc preprocessing may not be parallel-safe -- could try using shadow rules to enforce isolated temporary directories
#       for now, when you first run a subject, keep labels and hemis limited to 1 instance..  

wildcard_constraints:
    subject="[a-zA-Z0-9]+"

#this specifies that these rules always run locally (instead of submitting as jobs)
#localrules: all, gen_parc_cfg, copy_seed_to_diffparc

rule all:
    input: 
        connmap = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.mat',subject=subjects,label=labels,hemi=hemis,seed_res=seed_res,nsamples=nsamples),
        connmap_npz = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.npz',subject=subjects,label=labels,hemi=hemis,seed_res=seed_res,nsamples=nsamples,template=template)

rule gen_parc_cfg:
    input: 'cfg/parcellate_HarvardOxford.cfg'
    output: 'generated_cfg/parcellate_HarvardOxford_{label}_{hemi}_res-{seed_res}_nsamples-{nsamples}.cfg'
    params:
        parcellation_name = '{label}_{hemi}_nsamples-{nsamples}_res-{seed_res}',
        target_labels_txt = 'cfg/HarvardOxford_noremap.csv',
        target_mapping_txt = 'cfg/HarvardOxford_noremap.csv',
        bids_tags = 'seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}',
        seed_file = 'labels/t1/zona_warped_seg/{label}_{hemi}.nii.gz'
    shell:
        'cp -v {input} {output} &&'
        'echo parcellation_name={params.parcellation_name} >> {output} &&'
        'echo target_mapping_txt=`realpath {params.target_mapping_txt}` >> {output} &&'
        'echo target_labels_txt=`realpath {params.target_labels_txt}` >> {output} &&'
        'echo seed_file={params.seed_file} >> {output} &&'
        'echo bids_tags={params.bids_tags} >> {output}'

rule copy_seed_to_diffparc:
    input: join(seg_dir,'sub-{subject}/binary/sub-{subject}_exprater_{label}_{hemi}_bin.nii.gz')
    output: 'diffparc/work/sub-{subject}/labels/t1/zona_warped_seg/{label}_{hemi}.nii.gz'
    shell: 'cp -v {input} {output}'

rule run_diffparc:
    input: 
        bids_dir = bids_dir,
        in_prepdwi_dir = prepdwi_dir,
        parcellate_type = 'generated_cfg/parcellate_HarvardOxford_{label}_{hemi}_res-{seed_res}_nsamples-{nsamples}.cfg',
        seed = 'diffparc/work/sub-{subject}/labels/t1/zona_warped_seg/{label}_{hemi}.nii.gz'
    output:
        connmap = 'diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.mat',
        connmap_4d = 'diffparc/work/sub-{subject}/bedpost.{label}_{hemi}_nsamples-{nsamples}_res-{seed_res}/connMap.4d.nii.gz'
    params:
        app = f'{singularity_dir}/khanlab_diffparc-sumo_latest.sif',
        out_dir = 'diffparc',
        analysis_level = 'participant',
        participant_label = '{subject}',
        other_opts = '--skip_postproc --seed_res {seed_res} --nsamples {nsamples}',
    log: 'logs/run_diffparc/{subject}_{label}_{hemi}_res-{seed_res}_nsamples-{nsamples}.out'
    threads: 8
    resources: 
        mem_mb = 16000,
        time = 60*6 #6 hrs
    shell:
        'singularity run {params.app} {input.bids_dir} {params.out_dir} {params.analysis_level} --participant_label {params.participant_label} --in_prepdwi_dir {input.in_prepdwi_dir} --parcellate_type {input.parcellate_type} {params.other_opts} &> {log}'


rule split_connmap:
    input: 
        connmap_4d = 'diffparc/work/sub-{subject}/bedpost.{label}_{hemi}_nsamples-{nsamples}_res-{seed_res}/connMap.4d.nii.gz',
    params:
        connmap_3d_dir = 'diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap'
    output:
        connmap_3d = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap/conn_{tmap}.nii.gz',tmap=tmaps,allow_missing=True)
    envmodules: 'fsl'
    shell:
        'mkdir -p {params.connmap_3d_dir} && fslsplit {input.connmap_4d} {params.connmap_3d_dir}/conn_ -t'



rule transform_conn_to_template:
    input:
        connmap_3d = 'diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap/conn_{tmap}.nii.gz',
        affine =  lambda wildcards: glob(join(ants_template_dir,'snsx32_v0.2_i09sub-{subject}_acq-MP2RAGE_run-01_T1w*GenericAffine.mat'.format(subject=wildcards.subject))),
        warp =  lambda wildcards: glob(join(ants_template_dir,'snsx32_v0.2_i09sub-{subject}_acq-MP2RAGE_run-01_T1w*[0-9]Warp.nii.gz'.format(subject=wildcards.subject))),
        ref = join(ants_template_dir,'snsx32_v0.2_i09template0.nii.gz')
    output:
        connmap_3d = 'diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap/conn_{tmap}.nii.gz'
    envmodules: 'ants'
    shell:
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.connmap_3d} -o {output.connmap_3d} -r {input.ref} -t {input.warp} -t {input.affine}'


rule save_connmap_template_npz:
    input:
        connmap_3d = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap/conn_{tmap}.nii.gz',tmap=tmaps,allow_missing=True),
        mask = join(config['template_seg_dir'],'sub-SNSX32NLin2020Asym_hemi-{hemi}_desc-{label}_mask.nii.gz')
    output:
        connmap_npz = 'diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.npz'
    script: 'scripts/save_connmap_template_npz.py'


            

#this is slow:
rule merge_connmap:
    input:
        connmap_3d = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap/conn_{tmap}.nii.gz',tmap=tmaps,allow_missing=True)
    output:
        connmap_4d = 'diffparc/sub-{subject}/anat/sub-{subject}_space-{template}_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap_4d.nii.gz'
    envmodules: 'fsl'
    shell:
        'fslmerge -t {output} {input}'
    


