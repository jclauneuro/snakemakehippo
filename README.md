# Snakemake workflow: Zona Incerta Diffusion Parcellation

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for performing connectivity-based segmentation of the zona-incerta using probabilistic tractography.
It requires pre-processed DWI and bedpost (e.g. from prepdwi), and makes use of transforms from ANTS buildtemplate on the SNSX32 dataset.

## Authors

* Ali Khan (@akhanf)

## Usage

### Running on Graham

#### Step 1: Install neuroglia-helpers

`snakemake` and `snakemake_slurm` are wrappers in [neuroglia-helpers](http://github.com/khanlab/neuroglia-helpers), make sure you use the `graham_khanlab` cfg file.

#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`.

You can set the `test_single_subj` variable to True to reduce computing time when testing.

#### Step 3: 

Test your configuration by performing a dry-run via

    snakemake --use-singularity -n

Execute the workflow locally using an interactive job:
    
    regularInteractive -g  

    snakemake --use-singularity --cores 8 --resources gpus=1

To execute the workflow for all subjects, submitting a job for each rule group, use:

    snakemake_slurm

Alternatively, you can use the `snakemake_batch` wrapper to submit in N batches (e.g. 32), and submit all the batches at once, using the `gather_connmap_group` rule to split batches:

    snakemake_batch gather_connmap_group 32


To export files to dropbox, use:
    snakemake -s export_dropbox.smk


See the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/executable.html) for further details.

# Step 4: Investigate results

After successful execution, you can create a self-contained interactive HTML report with all results via:

    snakemake --report report.html

This report can, e.g., be forwarded to your collaborators.

### Advanced

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a personal or lab account.
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for the concrete project/run on your machine.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).
4. Modify the config, and any necessary sheets (and probably the workflow) as needed.
5. Commit any changes and push the project-branch to your fork on github.
6. Run the analysis.
7. Optional: Merge back any valuable and generalizable changes to the [upstream repo](https://github.com/snakemake-workflows/zona-diffparc) via a [**pull request**](https://help.github.com/en/articles/creating-a-pull-request). This would be **greatly appreciated**.
8. Optional: Push results (plots/tables) to the remote branch on your fork.
9. Optional: Create a self-contained workflow archive for publication along with the paper (snakemake --archive).
10. Optional: Delete the local clone/workdir to free space.


## Testing

No test cases yet

