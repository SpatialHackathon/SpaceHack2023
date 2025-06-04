# Contribution guide

First off, thank you for considering contributing to the SpaceHack 2023 project!

Following these guidelines helps to communicate that you respect the time of the developers managing and developing this open source project. In return, they should reciprocate that respect in addressing your issue, assessing changes, and helping you finalize your pull requests.

SpaceHack 2023 is an open source project and we love to receive contributions from our community — you! There are many ways to contribute, from writing tutorials or blog posts, improving the documentation, submitting bug reports and feature requests or writing modules (for data, methods, metrics) that can be incorporated into SpaceHack 2023 itself.

# Getting started

## Creating new modules

 1. First create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the module you want to implement and add the appropriate **tags**, and **assign** the task to yourself.
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

Easy!

# How to report a bug

When filing an issue, make sure to answer these five questions:

1. What version of Go are you using (go version)?
2. What operating system and processor architecture are you using?
3. What did you do?
4. What did you expect to see?
5. What did you see instead?

### Feature requests

If you find yourself wishing for a feature or modules that don't exist in SpaceHack 2023, you are probably not alone. There are bound to be others out there with similar needs. Many of the features that SpaceHack 2023 has today have been added because our users saw the need. Open an issue on our issues list on GitHub that describes the feature you would like to see, why you need it, and how it should work.


# Code review process

The core team looks at Pull Requests on an ad-hoc basis. 
After feedback has been given we expect responses within four weeks. After four weeks we may close the pull request if it isn't showing any activity.

# Community

You can chat with the core team on [Slack](https://spacehack20.slack.com).














