# SpaceHack - dataset modules

### Implementing a new dataset module

1. Create or claim a **GitHub issue** from the [SpaceHack issue board.](https://github.com/SpatialHackathon/SpaceHack2023/issues) that describes the dataset module you want to implement. There are currently around 20 datasets to implement, but if you come up with a new idea, please **create** a new issue, add the appropriate **tags**, and **assign** the task to yourself.
 2. Add **metadata** to our metadata [spreadsheet on the dataset tab](https://docs.google.com/spreadsheets/d/1QCeAF4yQG4bhZSGPQwwVBj_XF7ADY_2mK5xivAIfHsc/edit#gid=1453488771). Please fill in as much as you can as metadata is helpful! If you feel the need, please also add new columns or add additional notes.
 3. Now you are ready to create a new git **[branch](https://learngitbranching.js.org/)**. Try to give your new branch an intuitive prefix such as `data_...`. You can create a new branch in several ways: (i) [create a branch directly from the issue board](https://docs.github.com/en/issues/tracking-your-work-with-issues/creating-a-branch-for-an-issue) and then `git checkout` that branch, or (ii) via the command line:
```
# clone the template repository
git clone https://github.com/SpatialHackathon/SpaceHack2023.git
# create and switch to a new branch for your e.g. data "X"
git branch data_x_naveedishaque # try to make the branch name unique!
git checkout data_x_naveedishaque
# link the branch to the issue via the issue board: https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue
```
 4. Modify the files, filenames, and code in `template/`, referring to the examples in the `data` (for the LIBD Visium DLPFC dataset).
 5. Test. Before you make a pull request make sure that you (... do something...?). You should try to run your dataset through one of the two implemented methods (BayesSpace or SpaGCN). We are currently working on validators and automatic testing scripts... but this is tricky. Reach out to Niklas Muller-Botticher when you are ready to test!
 6. Create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request?tool=cli)
 7. Code review (by whom?) and merge your contributed module into the GitHub main branch!

For an example have a look [here](libd_dlpfc/).

### Dataset module layout and interface

Please define a conda recipe as a `yml` file for all the dependencies required for your script. This should list all the dependencies and also explicitly define the versions.

Input:
 - Output directory

Output:
 - Data organized by samples
 - Metadata

File structure:
```
dataset
├── sample_1
│   ├── coordinates.tsv
│   ├── counts.mtx
│   ├── features.tsv
│   ├── observations.tsv
│   ├── labels.tsv (optional)
│   └── H_E.{tiff,json} (optional)
├── sample_2
│   ├── …
│   └── …
├── experiment.json
└── samples.tsv
```
#### File structure (keep the headers the same - these are important for interfacing!)
head `coordinates.tsv`
```
     x    y
AAAA 1234 9876
ATAC 1357 9753
CAAG 3579 7531
```
`counts.mtx`. This should be in MatrixMarket format.

head `features.tsv`/`observations.tsv`
```
     row col selected
AAAA 1   1   true
ATAC 2   3   true
CAAG 5   2   false
```
The column `selected` is used for subsetting but is optional. `row` and `col` is needed in `observations.tsv` for bead-array based methods (Visium/ST).


head `labels.tsv` (annotations of the ground truth domain cluster annotations)
```
     label   label_confidence
AAAA Domain1 True
ATAC Domain1 True
CAAG Domain2 False
```
The column `label_confidence` is optional and used to indicate those cells 
and/or labels that are ground truth, if not all labels are high enough
confidence to be considered ground truth.


`image.tiff`. Images can be added in any format as appropriate (does not have to be tiff). If an image is available, please also add a json with relevant metadata (e.g. scale, but this might evolve during the hackathon)

`experiment.json`. Currently only technology (e.g. Visium, ST, MERSCOPE, MERFISH, Stereo-seq, Slide-seq, Xenium, STARmap, STARmap+, osmFISH, seqFISH) but more fields might be added.

`samples.tsv`. Sample directory and all relevant metadata, e.g. patient, replicate, slice, … and if applicable #clusters
