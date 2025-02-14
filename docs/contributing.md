# Creating a new module

Our workflow is set up to allow everyone to contribute "modules" in their preferred programming language (.. as long as that is either R or Python). A module can either be a dataset, a computational method, or an evaluation metric.
![Workflow](img/workflow.svg)

This repository contains some templates and examples of how to implement your module so that it interfaces seamlessly with other modules in the workflow. For example, if you want to implement a new method, you do not need to worry about input data or evaluation metrics as long as you follow the template for reading input and writing output - if you correctly adhere to the input and output guidelines, you should be able to interface with our default data modules and default evaluation metrics modules. The default modules are:

 - data: LIBD Visium DLPFC dataset (4 samples, each with 3 replicates)
 - methods: BayesSpace and SpaGCN
 - evaluation metrics: ARI and V

## How to contribute a module

Module contribution will be managed via GitHub. The steps to contribute a module are:

1. Fork (or if you are part of the SpaceHack community branch) the latest version of the [SpaceHack repository](https://github.com/SpatialHackathon/SpaceHack2023)

2. Make a copy of a [template]({{ repo_branch_url }}/templates/) depending on whether you are implementing a data, method, or metric module. You can have a look at existing modules if you are unsure what to do.

3. Modify the files, filenames, and code in your copied template and move it to the correct module directory.

4. Test. Before you make a pull request make sure that you can run the default modules. If it is a dataset, try to run e.g. SpaGCN or BayesSpace (or any other implemented method). If you implemented a method, try to run it on a dataset, and then test your output with an existing method.

5. Create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=cli).

6. Wait for code review and your module being merged into the workflow.

Easy!

## License

Currently, we have adopted the "MIT No Attribution" (MIT-0) License. 
Ideally, contributed modules should not have a license of their own, 
and therefore adopt the repository license.
More on MIT-0 [here](https://github.com/aws/mit-0).
