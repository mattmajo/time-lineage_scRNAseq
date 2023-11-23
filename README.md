# scRNA-seq analysis of blood induction and T-cell differentiation from PSCs

## Getting started
A preprocessed dataset is available at `./data_preprocessed/processed-2.h5ad` in Anndata format. Can be loaded via `scanpy.read('path/to/file')`.

### Git reposistory structure
If you're new to Git, see the ["Start using Git from the command line"](https://docs.gitlab.com/ee/gitlab-basics/start-using-git.html) documentation. This [Git cheatsheet](https://www.atlassian.com/git/tutorials/atlassian-git-cheatsheet) is also really useful.

To facilitate collaboration, easy sharing of results, and version control, each person should create a branch that they can update as their work progresses. First clone the repository to your workstation then create a branch with an easy identifier for the name (like your intials). For example, if my name was Andrew Hagner I would navigate to the repository and create the branch using:
```
git checkout -b AH
```
And then continue with my work. When I have some updates I want to add to the repository I use the command `git status` to see what has been changed. Then, I run:
```
git add new_or_updated_file.py
git commit -m "A quick description of the update"
git push origin AH
```
It's good practice to push small updates regularly rather than large ones infreqently. It makes it much easier to revert to a previous version without losing too much work when something catastrophic happens in your code.

### Conda environments and Jupyter Lab
This assumes that Anaconda is installed on your computer. If not, [install it from here](https://docs.anaconda.com/anaconda/install/). Details for how to manage conda environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

The file `environment.yml` contains the channels and dependencies needed to create a basic conda environment to start working with the data. To use it, open your terminal and navigate to the directory containing the file and type the command:
```
conda env create -f environment.yml
```
The environment will be named 'scanpy' and can be activated using `conda activate scanpy`. To access the environment from within Jupyter Lab/Notebook, activate it then type:
```
ipython kernel install --user --name=scanpy
```
Likewise, other environments can be added to Jupyter but ipykernel must be installed before running that last command (use `conda install ipykernel`).