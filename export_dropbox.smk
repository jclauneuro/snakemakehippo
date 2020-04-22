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
from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
DBox = DropboxRemoteProvider(oauth2_access_token=config['dropbox_token'])


wildcard_constraints:
    subject="[a-zA-Z0-9]+"


localrules: save_connmap,save_clusters

rule all:
    input: 
        clusters = expand('remote/diffparc/clustering/group_space-{template}_seed-{seed}_hemi-{hemi}_method-spectralcosine_k-{k}_cluslabels.nii.gz',seed=seeds,hemi=hemis,template=config['template'],k=range(2,config['max_k']+1)),
        connmaps = expand('remote/diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz',seed=seeds,hemi=hemis,template=config['template'])

rule save_connmap:
    input: 'diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz' 
    output: DBox.remote('remote/diffparc/connmap/group_space-{template}_seed-{seed}_hemi-{hemi}_connMap.npz') 
    shell: 'cp {input} {output}'

rule save_clusters:
    input: 'diffparc/clustering/group_space-{template}_seed-{seed}_hemi-{hemi}_method-spectralcosine_k-{k}_cluslabels.nii.gz'
    output: DBox.remote('remote/diffparc/clustering/group_space-{template}_seed-{seed}_hemi-{hemi}_method-spectralcosine_k-{k}_cluslabels.nii.gz')
    shell: 'cp {input} {output}'


