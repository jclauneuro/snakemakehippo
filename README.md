# Snakemake workflow for SNSX32 hippocampal template creation.

[![Snakemake](https://img.shields.io/badge/snakemake-≥3.12.0-brightgreen.svg)](https://snakemake.bitbucket.io)

This is a Snakemake workflow for propagating hippocampal segmentations from Jordan Dekraker's AutoTop: https://github.com/jordandekraker/Hippocampal_AutoTop
It requires pre-processed T2SPACE data, and makes use of transforms from ANTS buildtemplate on the SNSX32 dataset.

## Authors

* Jonathan C. Lau (@jclauneuro)

## Usage


### Step 1: Install Snakemake

Install Snakemake using [conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html):

    conda create -c bioconda -c conda-forge -n snakemake snakemake

For installation details, see the [instructions in the Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).


#### Step 2: Configure workflow

Configure the workflow according to your needs via editing the file `config.yaml`, and editing the `participants.tsv`

#### Step 3: Dry-run

Test your configuration by performing a dry-run via

    snakemake -np

#### Step 4: Execution on graham (compute canada)

There are a few different ways to execute the workflow:
  1. Execute the workflow locally using an interactive job
  2. Execute the workflow using the `cc-slurm` profile

##### Interactive Job

Execute the workflow locally using an interactive job:

    salloc --time=3:00:00 --gres=gpu:t4:1 --cpus-per-task=8 --ntasks=1 --mem=32000 --account=YOUR_CC_ACCT srun snakemake --use-singularity --cores 8 --resources gpus=1 mem_mb=32000 

##### Use the cc-slurm profile

The cc-slurm profile sets up default options for running on compute canada systems. More info in the README here: https://github.com/khanlab/cc-slurm

If you haven't used it before, deploy the cc-slurm profile using:

    cookiecutter gh:khanlab/cc-slurm -o ~/.config/snakemake -f    

Note: you must have cookiecutter installed (e.g. `pip install cookiecutter`)

Then to execute the workflow for all subjects, submitting a job for each rule group, use:

    snakemake --profile cc-slurm


##### Export to Dropbox

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

