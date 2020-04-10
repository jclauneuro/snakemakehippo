from os.path import join
from glob import glob
import pandas as pd

configfile: 'config.yml'


#currently this isn't used:
participants_tsv = join(config['bids_dir'],'participants.tsv')
subjects_table = pd.read_table(participants_tsv)

#get list of subjects based on seed_seg_dir
subjects = os.listdir(config['seed_seg_dir'])
subjects = [ s.strip('sub-') for s in subjects ]
subjects = sorted(subjects)


#get list of ROIs
f = open(config['targets_txt'],'r')
targets = [line.strip() for line in f.readlines()]
f.close()

#get seeds and hemis from config
seeds = config['seeds']
hemis = config['hemis']



#for saving data to dropbox (for colab notebook)
if config['enable_dropbox'] == True:
    from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
    DBox = DropboxRemoteProvider(oauth2_access_token=config['dropbox_token'])

#just test with one subject:
if config['test_single_subj'] == True:
    subjects = subjects[0]


wildcard_constraints:
    subject="[a-zA-Z0-9]+"

#this specifies that these rules always run locally (instead of submitting as jobs)
#localrules: all, gen_parc_cfg, copy_seed_to_diffparc



rule all:
    input: 
        connmap_group_npz = expand('diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz',seed=seeds,hemi=hemis,template=config['template'])



rule import_targets:
    input: 
        lh = join(config['targets_seg_dir'],config['targets_seg_lh']),
        rh = join(config['targets_seg_dir'],config['targets_seg_rh'])
    output: 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    envmodules: 'fsl'
    singularity: config['singularity_neuroglia']
    log: 'logs/import_targets_hcp_mmp_sym/sub-{subject}.log'
    shell:
        'fslmaths {input.lh} -max {input.rh} {output} &> {log}'

rule import_template_seed:
    input: join(config['template_seg_dir'],config['template_seg_nii'])
    output: 'diffparc/template_masks/sub-{template}_hemi-{hemi}_desc-{seed}_mask.nii.gz'
    log: 'logs/import_template_seed/{template}_{seed}_{hemi}.log'
    shell: 'cp -v {input} {output} &> {log}'


rule transform_to_subject:
    input: 
        seed = 'diffparc/template_masks/sub-{template}_hemi-{hemi}_desc-{seed}_mask.nii.gz',
        affine = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(**wildcards))),
        invwarp = lambda wildcards: glob(join(config['ants_template_dir'],config['ants_invwarp_nii'].format(**wildcards))),
        ref = 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    output: 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}.nii.gz'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_to_subject/{template}_sub-{subject}_{seed}_{hemi}.log'
    shell:
        'antsApplyTransforms -d 3 --interpolation NearestNeighbor -i {input.seed} -o {output} -r {input.ref} -t [{input.affine},1] -t {input.invwarp} &> {log}'
    
    
rule resample_targets:
    input: 
        dwi = join(config['prepdwi_dir'],'bedpost','sub-{subject}','mean_S0samples.nii.gz'),
        targets = 'diffparc/sub-{subject}/masks/lh_rh_targets_native.nii.gz'
    params:
        seed_resolution = config['probtrack']['seed_resolution']
    output:
        mask = 'diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz',
        mask_res = 'diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz',
        targets_res = 'diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_targets/sub-{subject}.log'
    group: 'split_targets'
    shell:
        'fslmaths {input.dwi} -bin {output.mask} &&'
        'mri_convert {output.mask} -vs {params.seed_resolution} {params.seed_resolution} {params.seed_resolution} {output.mask_res} -rt nearest &&'
        'reg_resample -flo {input.targets} -res {output.targets_res} -ref {output.mask_res} -NN 0  &> {log}'

rule resample_seed:
    input: 
        seed = 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}.nii.gz',
        mask_res = 'diffparc/sub-{subject}/masks/brain_mask_dwi_resampled.nii.gz'
    output:
        seed_res = 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}_resampled.nii.gz',
    singularity: config['singularity_neuroglia']
    log: 'logs/resample_seed/{template}_sub-{subject}_{seed}_{hemi}.log'
    shell:
        'reg_resample -flo {input.seed} -res {output.seed_res} -ref {input.mask_res} -NN 0 &> {log}'

    
    

rule split_targets:
    input: 
        targets = 'diffparc/sub-{subject}/masks/lh_rh_targets_dwi.nii.gz',
    params:
        target_num = lambda wildcards: targets.index('{target}'.format(**wildcards)) + 1
    output:
        target_seg = 'diffparc/sub-{subject}/targets/{target}.nii.gz'
    singularity: config['singularity_neuroglia']
    log: 'logs/split_targets/sub-{subject}/{target}.log'
    group: 'split_targets'
    shell:
        'fslmaths {input} -thr {params.target_num} -uthr {params.target_num} {output} &> {log}'

rule gen_targets_txt:
    input:
        target_seg = expand('diffparc/sub-{subject}/targets/{target}.nii.gz',target=targets,allow_missing=True)
    output:
        target_txt = 'diffparc/sub-{subject}/target_images.txt'
    log: 'logs/get_targets_txt/sub-{subject}.log'
    run:
        f = open(output.target_txt,'w')
        for s in input.target_seg:
            f.write(f'{s}\n')
        f.close()


rule run_probtrack:
    input:
        seed_res = 'diffparc/sub-{subject}/masks/seed_from-{template}_{seed}_{hemi}_resampled.nii.gz',
        target_txt = 'diffparc/sub-{subject}/target_images.txt',
        mask = 'diffparc/sub-{subject}/masks/brain_mask_dwi.nii.gz'
    params:
        bedpost_merged = join(config['prepdwi_dir'],'bedpost','sub-{subject}','merged'),
        probtrack_opts = config['probtrack']['opts'],
        probtrack_dir = 'diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}'
    output:
        target_seg = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}.nii.gz',target=targets,allow_missing=True)
    singularity: config['singularity_neuroglia']
    threads: 4
    resources: 
        mem_mb = 16000, #16GB 
        time = 120 #2 hrs
    log: 'logs/run_probtrack/{template}_sub-{subject}_{seed}_{hemi}.log'
    shell:
        'probtrackx2 --samples={params.bedpost_merged}  --mask={input.mask} --seed={input.seed_res} ' 
        '--targetmasks={input.target_txt} --seedref={input.seed_res} --nsamples={config[''probtrack''][''nsamples'']} ' 
        '--dir={params.probtrack_dir} {params.probtrack_opts} -V 2  &> {log}'


rule transform_conn_to_template:
    input:
        connmap_3d = 'diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}.nii.gz',
        affine =  lambda wildcards: glob(join(config['ants_template_dir'],config['ants_affine_mat'].format(subject=wildcards.subject))),
        warp =  lambda wildcards: glob(join(config['ants_template_dir'],config['ants_warp_nii'].format(subject=wildcards.subject))),
        ref = join(config['ants_template_dir'],config['ants_ref_nii'])
    output:
        connmap_3d = 'diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}_space-{template}.nii.gz'
    envmodules: 'ants'
    singularity: config['singularity_neuroglia']
    log: 'logs/transform_conn_to_template/sub-{subject}_{seed}_{hemi}_{template}/{target}.log'
    group: 'post_track'
    shell:
        'antsApplyTransforms -d 3 --interpolation Linear -i {input.connmap_3d} -o {output.connmap_3d} -r {input.ref} -t {input.warp} -t {input.affine} &> {log}'


rule save_connmap_template_npz:
    input:
        connmap_3d = expand('diffparc/sub-{subject}/probtrack_{template}_{seed}_{hemi}/seeds_to_{target}_space-{template}.nii.gz',target=targets,allow_missing=True),
        mask = 'diffparc/template_masks/sub-{template}_hemi-{hemi}_desc-{seed}_mask.nii.gz'
    output:
        connmap_npz = 'diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz'
    log: 'logs/save_connmap_to_template_npz/sub-{subject}_{seed}_{hemi}_{template}.log'
    group: 'post_track'
    script: 'scripts/save_connmap_template_npz.py'

rule gather_connmap_group:
    input:
        connmap_npz = expand('diffparc/sub-{subject}/connmap/sub-{subject}_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz',subject=subjects,allow_missing=True)
    output:
        connmap_group_npz = 'diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz'
    log: 'logs/gather_connmap_group/{seed}_{hemi}_{template}.log'
    run:
        import numpy as np
        
        #load first file to get shape
        data = np.load(input['connmap_npz'][0])
        affine = data['affine']
        mask = data['mask']
        conn_shape = data['conn'].shape
        nsubjects = len(input['connmap_npz'])
        conn_group = np.zeros([nsubjects,conn_shape[0],conn_shape[1]])
        
        for i,npz in enumerate(input['connmap_npz']):
            data = np.load(npz)
            conn_group[i,:,:] = data['conn']
            
        #save conn_group, mask and affine
        np.savez(output['connmap_group_npz'], conn_group=conn_group,mask=mask,affine=affine)
 
                   
localrules: save_to_dropbox

rule save_to_dropbox:
    input: 
        connmap_group_npz = 'diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz' if config['enable_dropbox']==True else ''
    output: 
        connmap_group_npz_dropbox = DBox.remote('diffparc_hcp/remote/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz') if config['enable_dropbox']==True else ''
    shell: 'cp {input} {output}'
