# SpaceHack - method modules

### Implementing a new dataset module

1. Create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the dataset module you want to implement. There are currently around 30 methods to implement, but if you come up with a new idea, please **create** a new issue, add the appropriate **tags**, and **assign** the task to yourself.
 2. Add **metadata** to our metadata [spreadsheet on the methods tab](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=0). Please fill in as much as you can as metadata is helpful! If you feel the need, please also add new columns or add additional notes.
 3. Now you are ready to create a new git **[branch](https://learngitbranching.js.org/)**. Try to give your new branch an intuitive prefix such as `method...`. You can create a new branch in several ways: (i) [create a branch directly from the issue board](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-a-branch-for-an-issue) and then `git checkout` that branch, or (ii) via the command line:
```
# clone the template repository
git clone https://github.com/SpatialHackathon/SpaceHack2023.git
# create and switch to a new branch for your e.g. method "SpaGCN"
git branch method_SpaGCN_naveedishaque # try to make the branch name unique!
git checkout method_SpaGCN_naveedishaque
# link the branch to the issue via the issue board: https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue
```
 4. Modify the files, filenames, and code in `template/`, referring to the examples in the `method` (for the SpaGCN and BayesSpace methods).
 5. Test. Before you make a pull request make sure that you (... do something...?). You should try to run your method on the default dataset (LIBD Visium DLPFC) and run the results through one of the default metrics (ARI or V). We are currently working on validators and automatic testing scripts... but this is tricky. Reach out to Niklas Muller-Botticher when you are ready to test!
 6. Create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=cli)
 7. Code review (by whom?) and merge your contributed module into the GitHub main branch!

For examples have a look [here for a method in Python](spaGCN/) or [here for a method in R](BayesSpace/).

### Method module layout and interface

Please define a conda recipe as a `yml` file. This should list all the dependencies and also explicitly define the versions.

Input:
 - Coordinates
 - Counts
 - …

Output:
 - Predicted labels
 - Embedding (optional)
