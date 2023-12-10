# SpaceHack - templates for modules

This directory contains the templates that you can use to implement a new dataset, method or metric.

For further instructions have a look in the corresponding directory;

- [data](/data)
- [method](/method)
- [metric](/metric)

### How to contribute a module

Module contribution will be managed via GitHub. The steps to contribute a module are:
 1. Create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the module you want to implement. There are currently 90 issues to claim, but if you come up with a new idea, please **create** a new issue, add the appropriate **tags**, and **assign** the task to yourself.
 2. Add **metadata** to our metadata [spreadsheet](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit). Please fill in as much as you can as metadata is helpful! If you feel the need, please also add new columns or add additional notes. The metadata should be added to the appropriate tabs:
    - [datasets](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=1453488771)
    - [computational methods](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=0)
    - [evaluation metrics](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=4776337)
    - [simulations and technical evaluation](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=640974611)
 3. Now you are ready to create a new git **[branch](https://learngitbranching.js.org/)**. Try to give your new branch an intuitive prefix such as `data_...`, `method_...` or `metric_..`. You can create a new branch in several ways: (i) [create a branch directly from the issue board](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-a-branch-for-an-issue) and then `git checkout` that branch, or (ii) via the command line:
```
# clone the template repository
git clone https://github.com/SpatialHackathon/SpaceHack2023.git
# create and switch to a new branch for your e.g. method "X"
git branch method_x_naveedishaque # try to make the branch name unique!
git checkout method_x_naveedishaque
# link the branch to the issue via the issue board: https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue
```
 4. Modify the files, filenames, and code in `template/`, referring to the examples in the `data`, `method`, or `metric` subfolder. If your method requires a specific type or preprocessing, please reach out to the organisers!
 5. Test. We are currently working on validators and automatic testing scripts... but this is tricky. Reach out to Niklas Muller-Botticher when you are ready to test!
 6. Create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=cli)
 7. Code review (by whom?) and merge your contributed module into the GitHub main branch!
