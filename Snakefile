from os.path import join

configfile: 'config.yml'

prepdwi_dir = config['prepdwi_dir']
seg_dir = config['seg_dir']
bids_dir = config['bids_dir']
singularity_dir = config['singularity_dir']

seed_res = config['probtrack']['seed_res']
nsamples = config['probtrack']['nsamples']


subjects = ['C001']
labels = ['ZI']
hemis = ['L']

wildcard_constraints:
    subject="[a-zA-Z0-9]+"

#this specifies that these rules always run locally (instead of submitting as jobs)
localrules: all, gen_parc_cfg, copy_seed_to_diffparc

rule all:
    input: 
        connmap = expand('diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.mat',subject=subjects,label=labels,hemi=hemis,seed_res=seed_res,nsamples=nsamples)

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
        connmap = 'diffparc/sub-{subject}/anat/sub-{subject}_space-T1w_res-{seed_res}_seed-{label}_hemi-{hemi}_targets-cortical_nsamples-{nsamples}_connMap.mat'
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
        'singularity run {params.app} {input.bids_dir} {params.out_dir} {params.analysis_level} --participant_label {params.participant_label} --parcellate_type {input.parcellate_type} {params.other_opts} &> {log}'



