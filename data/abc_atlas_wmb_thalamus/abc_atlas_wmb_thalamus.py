#!/usr/bin/env python

# Author_and_contribution: Niklas Mueller-Boetticher; created template
# Author_and_contribution: Thomas Chartrand; co-wrote initial versions of 
#                                            thalamus subsetting functions
# Author_and_contribution: Meghan Turner; wrote code in file, co-wrote 
#                                         initial versions of thalamus 
#                                         subsetting functions

import argparse
from pathlib import Path
import json

import numpy as np
import pandas as pd
import anndata as ad
import scipy
# additional packages are imported before the functions in which they are used

'''
This file Contains code for loading ABC Atlas WMB TH dataset for SpaceHack2.0,
including functions to:
- Subset the ABC Whole Mouse Brain (WMB) Atlas to a thalamus (TH) subset
- Download the ABC Atlas WMB data from the AWS S3 bucket
- Add confidence flag for CCFv3/ARA annotations
- Convert dataset to SpaceHack format
and executable code to write the dataset to disk in the SpaceHack2.0 directory
structure.
'''

# versions for downloading ABC Atlas WMB data
CURRENT_VERSION = '20230830'
BRAIN_LABEL = 'C57BL6J-638850'
CCF_VERSION = '20230630'





''' 
--------------------------------------------------------------------------------
THALAMUS SUBSET FUNCTIONS ------------------------------------------------------
--------------------------------------------------------------------------------
Functions for loading a TH + ZI (thalamus + zona incerta) subset of the 
ABC (Allen Brain Cell) Atlas WMB (Whole Mouse Brain) MERFISH dataset.
- load_th_subset_adata: load a TH+ZI subset of the ABC Atlas MERFISH dataset
- filter_adata_by_class: remove non-neuronal & non-TH neuronal classes
- label_thalamus_spatial_subset: label cells that are in the TH+ZI subset
- sectionwise_dilation: dilate a stack of 2D binary masks by a specified radius
- cleanup_mask_regions: remove too-small mask regions that are likely mistakes
- label_thalamus_masked_cells: label cells that are inside the TH+ZI CCF mask
'''

import scipy.ndimage as ndi
import nibabel

def load_th_subset_adata(version=CURRENT_VERSION, specimen=BRAIN_LABEL, 
                         counts_transform='raw'):
    '''Load a TH + ZI (thalamus + zona incerta) subset of the ABC (Allen 
    Brain Cell) Atlas WMB (Whole Mouse Brain) MERFISH dataset as an AnnData
    object.
    
    Parameters
    ----------
    version : str, default=CURRENT_VERSION
        which release version of the ABC Atlas to load
    specimen : str, default=BRAIN_LABEL
        which specimen to load from the ABC Atlas
    counts_transform : {'log2', 'raw'}, default='raw'
        string specifying which transformation of the gene expression counts
        to load from the expression matrices
    '''
    flip_y = True  # so coronal sections display dorsal/superior side up
    field_name='TH_ZI_dataset'
    
    # we'll combine everything (counts+metadata) into an AnnData obj to 
    # ensure subsetting is applied consistently across all data
    adata = get_counts_adata(version=version, specimen=specimen, 
                             counts_transform=counts_transform)
   
    # Load cell (observation) metadata
    cells_md_df = get_cell_metadata_df(version=version, specimen=specimen,
                                       flip_y=flip_y, drop_unused=True)
    # label TH+ZI subset of cells
    cells_md_df = label_thalamus_spatial_subset(cells_md_df,
                                                flip_y=flip_y,
                                                field_name=field_name)
    # subset cell metadata to just TH+ZI labeled cells
    cells_md_df = cells_md_df[cells_md_df[field_name]].copy().drop(columns=[field_name])
    # use cell_labels to correctly subset counts to match cell metadata
    cell_labels = cells_md_df.index
    adata = adata[adata.obs_names.intersection(cell_labels)]
    # add cell metadata to obs
    adata.obs = adata.obs.join(cells_md_df[cells_md_df.columns.difference(adata.obs.columns)])
    
    # Load gene (feature) metadata
    gene_md_df = get_gene_metadata_df(version=version, specimen=specimen)
    # remove 'blank' codenames from genes list
    gene_md_df = gene_md_df[~gene_md_df.index.str.contains('Blank')]
    # use gene_symbol to correctly subset by genes
    adata.var = adata.var.reset_index().set_index('gene_symbol')
    gene_symbols = gene_md_df.index
    adata = adata[:, adata.var_names.intersection(gene_symbols)]
    # add cell metadata to obs
    adata.var = adata.var.join(gene_md_df[gene_md_df.columns.difference(adata.var.columns)])
    # reset index back to gene_idetifier and rename columns as preferred by
    # SpaceHack community
    adata.var.rename(columns={'gene_identifier':'gene_id'}, inplace=True)
    adata.var = adata.var.reset_index().set_index('gene_id')
    adata.var.rename(columns={'gene_symbol':'gene_name'}, inplace=True)

    # filter out non-neuronal cells, plus non-thalamic classes
    adata = filter_adata_by_class(adata, filter_nonneuronal=True,
                                  filter_midbrain=True)

    # curate sections to those with reasonable alignment to CCF parcellation
    adata = adata[( (adata.obs['z_reconstructed'] > 6.3)
                    & (adata.obs['z_reconstructed'] < 6.9) 
                  )
                  | ( (adata.obs['z_reconstructed'] > 7.1)
                    & (adata.obs['z_reconstructed'] < 8.1) 
                    )]
    
    return adata
    

def filter_adata_by_class(th_zi_adata, filter_nonneuronal=True,
                          filter_midbrain=True):
    ''' Filters anndata object to only include cells from specific taxonomy 
    classes. Defaults to include only neuronal cells from TH+ZI classes

    Parameters
    ----------
    th_zi_adata
        anndata object containing the ABC Atlas MERFISH dataset
    filter_nonneuronal : bool, default=True
        filters out non-neuronal classes
    filter_midbrain : bool, default=True
        filters out midbrain classes; may be useful to keep these if interested
        in analyzing midbrain-thalamus boundary in the anterior
    '''
    # hardcoded class categories
    th_zi_dataset_classes = ['12 HY GABA', '17 MH-LH Glut', '18 TH Glut']
    midbrain_classes = ['19 MB Glut', '20 MB GABA']
    nonneuronal_classes = ['30 Astro-Epen', '31 OPC-Oligo', '33 Vascular',
                           '34 Immune']

    # always keep th_zi_dataset_classes
    classes_to_keep = th_zi_dataset_classes.copy()

    # optionally include midbrain and/or nonneuronal classes
    if not filter_midbrain:
        classes_to_keep += midbrain_classes
    if not filter_nonneuronal:
        classes_to_keep += nonneuronal_classes

    th_zi_adata = th_zi_adata[th_zi_adata.obs['class'].isin(classes_to_keep)]
    return th_zi_adata

def label_thalamus_spatial_subset(cells_df, flip_y=False, 
                                  field_name='TH_ZI_dataset',
                                  distance_px=20):
    '''Labels cells that are in the thalamus spatial subset of the ABC atlas.
    
    Turns a rasterized image volume that includes all thalamus (TH) and zona
    incerta (ZI) CCF structures in a binary mask, then dilates by 200um (20px)
    to ensure inclusion of the vast majority cells in known thalamic subclasses.
    Labels cells that fall in this dilate binary mask as in the 'TH_ZI_dataset' 
    
    Parameters
    ----------
    cells_df : pandas dataframe
        dataframe of cell metadata
    flip_y : bool, default=False
        flip y-axis orientation of th_mask so coronal section is dorsal-side up.
        MUST be set to true if flip_y=True in get_cell_metadata_df()
    field_name : str, default='TH_ZI_dataset'
        name for column containing the thalamus dataset boolean flag
    distance_px : int, default=20
        dilation radius in pixels (1px = 10um)
        
    Returns
    -------
    cells_df 
        with a new boolean column specifying which cells are in the TH+ZI dataset
    '''
    # use reconstructed (in MERFISH space) coordinates from cells_df
    coords = ['x_reconstructed','y_reconstructed','z_reconstructed']
    resolutions = np.array([10e-3, 10e-3, 200e-3])
    
    # load 'resampled CCF' (rasterized, in MERFISH space) image volumes from the
    # ABC Atlas dataset (z resolution limited to merscope slices)
    ccf_img = get_ccf_labels_image()
    
    # ccf_img voxels are labelled by brain structure parcellation_index;
    # get a list of all indices that correspond to TH or ZI (sub)structures
    ccf_df = get_ccf_terms_df()
    th_zi_ind = np.hstack(
                    (ccf_df.loc[ccf_df['parcellation_term_acronym']=='TH', 
                                'parcellation_index'].unique(),
                     ccf_df.loc[ccf_df['parcellation_term_acronym']=='ZI', 
                                'parcellation_index'].unique())
                )
    
    # generate binary mask
    th_mask = np.isin(ccf_img, th_zi_ind) # takes about 5 sec
    # flip y-axis to match flipped cell y-coordinates
    if flip_y:
        th_mask = np.flip(th_mask, axis=1)
    # dilate by 200um to try to capture more TH/ZI cells
    mask_img = sectionwise_dilation(th_mask, distance_px, true_radius=False)
    # remove too-small mask regions that are likely mistaken parcellations
    mask_img = cleanup_mask_regions(mask_img, area_ratio_thresh=0.1)
    # label cells that fall within dilated TH+ZI mask; by default, 
    cells_df = label_thalamus_masked_cells(cells_df, mask_img, coords,  
                                           resolutions, field_name=field_name)
    # exclude the 1 anterior-most and 1 posterior-most thalamus sections due to
    # poor overlap between mask & thalamic cells
    cells_df[field_name] = (cells_df[field_name] 
                            & (4.81 < cells_df[coords[2]]) 
                            & (cells_df[coords[2]] < 8.39))

    return cells_df

def sectionwise_dilation(mask_img, distance_px, true_radius=False):
    '''Dilates a stack of 2D binary masks by a specified radius (in px).
    
    Parameters
    ----------
    mask_img : array_like
        stack of 2D binary mask, shape (x, y, n_sections)
    distance_px : int
        dilation radius in pixels
    true_radius : bool, default=False
        specifies the method used by ndimage's binary_dilation to dilate
          - False: dilates by 1 px per iteration for iterations=distance_px
          - True: dilates once using a structure of radius=distance_px
        both return similar results but true_radius=False is significantly faster
    '''
    dilated_mask_img = np.zeros_like(mask_img)
    
    if true_radius:
        # generate a circular structure for dilation
        coords = np.mgrid[-distance_px:distance_px+1, -distance_px:distance_px+1]
        struct = np.linalg.norm(coords, axis=0) <= distance_px
        
    for i in range(mask_img.shape[2]):
        if true_radius:
            dilated_mask_img[:,:,i] = ndi.binary_dilation(mask_img[:,:,i], 
                                                          structure=struct)
        else:
            dilated_mask_img[:,:,i] = ndi.binary_dilation(mask_img[:,:,i], 
                                                           iterations=distance_px)
    return dilated_mask_img


def cleanup_mask_regions(mask_img, area_ratio_thresh=0.1):
    ''' Removes, sectionwise, any binary mask regions whose areas are smaller
    than the specified ratio, as compared to the largest region in the mask.
    
    Parameters
    ----------
    mask_img : array_like
        stack of 2D binary mask, shape (x, y, n_sections)
    area_ratio_thresh : float, default=0.1
        threshold for this_region:largest_region area difference ratio; removes
        any regions smaller than this threshold
    '''
    new_mask_img = np.zeros_like(mask_img)
    for sec in range(mask_img.shape[2]):
        mask_2d = mask_img[:,:,sec]
        labeled_mask, n_regions = ndi.label(mask_2d)

        # calculate the area of the largest region
        largest_region = np.argmax(ndi.sum(mask_2d, labeled_mask, 
                                           range(n_regions+1)))
        largest_area = np.sum(labeled_mask==largest_region)

        # filter out regions with area ratio smaller than the specified threshold
        regions_to_keep = [label for label 
                           in range(1, n_regions+1) 
                           if ( (np.sum(labeled_mask==label) / largest_area) 
                                >= area_ratio_thresh
                              )
                          ]
        # make a new mask with only the remaining objects
        new_mask_img[:,:,sec] = np.isin(labeled_mask, regions_to_keep)

    return new_mask_img


def label_thalamus_masked_cells(cells_df, mask_img, coords, resolutions,
                                field_name='TH_ZI_dataset'):
    '''Labels cells that are inside the TH+ZI mask from the CCF parcellation.
    
    Parameters
    ----------
    cells_df : pandas dataframe
        dataframe of cell metadata
    mask_img : array_like
        stack of 2D binary masks, shape (x, y, n_sections)
    coords : list
        column names in cells_df that contain the cells xyz coordinates, 
        list of strings of length 3
    resolutions : array
        xyz resolutions used to compare coords to mask_img positions
    field_name : str, default='TH_ZI_dataset'
        name for column containing the thalamus dataset boolean flag
    drop_end_sections : bool, default=True
        drops the anterior-most section and posterior-most section, which
        contain CCF thalamus parcellation labels but have poor overlap with
        thalamic cell types

    Returns
    -------
    cells_df 
        with a new boolean column specifying which cells are in the thalamus dataset
    '''
    coords_index = np.rint(cells_df[coords].values / resolutions).astype(int)
    # tuple() makes this like calling mask_img[coords_index[:,0], coords_index[:,1], coords_index[:,2]]
    cells_df[field_name] = mask_img[tuple(coords_index.T)]
    
    return cells_df


'''
--------------------------------------------------------------------------------
ABC ATLAS DOWNLOAD FUNCTIONS ---------------------------------------------------
--------------------------------------------------------------------------------
Functions for downloading the ABC Atlas MERFISH dataset from the AWS S3 
bucket.
- get_counts_adata: load gene expression counts matrix as AnnData object
- get_cell_metadata_df: load cell (obs) metadata as DataFrame
- get_gene_metadata_df: load gene (var / feature) metadata as DataFrame
- get_ccf_labels_image: load rasterized image volumes of CCFv3 parcellation
- get_ccf_terms_df: load CCF parcellation term membership metadata as DataFrame
'''
import tempfile
import boto3
from botocore import UNSIGNED
from botocore.client import Config

# AWS S3 bucket: https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html
S3_BUCKET_NAME = 'allen-brain-cell-atlas'
S3_REGION_NAME = 'us-west-2'

def get_counts_adata(version=CURRENT_VERSION, specimen=BRAIN_LABEL, 
                     counts_transform='raw'):
    '''Load gene expression counts matrix as AnnData obj directly into 
    memory using tempfile.'''
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))
        
        # Load one set of gene expr counts, raw or log2 (takes ~2min + ~3GB)
        counts_file = f'{specimen}-{counts_transform}.h5ad'
        counts_dir = f'expression_matrices/MERFISH-{specimen}/{version}/'
        counts_temp_path = f'{temp_dir}/{counts_file}'
        s3_client.download_file(S3_BUCKET_NAME, counts_dir+counts_file, 
                                counts_temp_path)
        adata = ad.read_h5ad(counts_temp_path) #, backed='r')
    
    return adata

    
def get_cell_metadata_df(version=CURRENT_VERSION, specimen=BRAIN_LABEL,
                         drop_unused=True, flip_y=True):
    '''Load cell (obs) metadata as DataFrame from CSV directly into memory
    using tempfile.
    Parameters
    ----------
    drop_unused : bool, default=True
        don't load unused columns (color, other xyz coords, etc)
    flip_y : bool, default=True
        flip y_reconstructed coords so up is positive
        
    '''
    # Load only the columns we will use (e.g. load only "reconstructed" xyz 
    # coords as they are in the same coordinate system as the CCF image 
    # volumes; only load some color columns)
    float_columns = [
        'average_correlation_score',
        'x_reconstructed', 'y_reconstructed', 'z_reconstructed', # millimeters
        ]
    cat_columns = [
        'brain_section_label', 'cluster_alias', 
        'neurotransmitter', 'class',
        'subclass', 'supertype', 'cluster', 
        'class_color', 'subclass_color', 'supertype_color', 'cluster_color',
        'parcellation_index',  
        'parcellation_division',
        'parcellation_structure', 'parcellation_substructure',
        'parcellation_substructure_color'
        ]
    dtype = dict(cell_label='string', 
                **{x: 'float' for x in float_columns}, 
                **{x: 'category' for x in cat_columns})
    usecols = list(dtype.keys()) if drop_unused else None
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))

        # Load metadata as pandas DataFrames from CSV file
        cell_metadata_dir = f'metadata/MERFISH-{specimen}-CCF/{version}/views/'
        cell_metadata_file = 'cell_metadata_with_parcellation_annotation.csv'
        cell_md_temp_path = f'{temp_dir}/{cell_metadata_file}'
        s3_client.download_file(S3_BUCKET_NAME, 
                                cell_metadata_dir+cell_metadata_file, 
                                cell_md_temp_path)
        cell_md_df = pd.read_csv(cell_md_temp_path, dtype=dtype, 
                                 usecols=usecols, index_col='cell_label',
                                 engine='c')
    
    # round z_reconstructed coords to 10ths place to correct overprecision
    cell_md_df['z_reconstructed'] = cell_md_df['z_reconstructed'].round(1)
    
    if flip_y:
        cell_md_df['y_reconstructed'] *= -1
    
    return cell_md_df


def get_gene_metadata_df(version=CURRENT_VERSION, specimen=BRAIN_LABEL):
    '''Load gene (var / feature) metadata as DataFrame from CSV directly
    into memory using tempfile.'''
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))

        # Load metadata as pandas DataFrames from CSV file
        gene_metadata_dir = f'metadata/MERFISH-{specimen}/{version}/'
        gene_metadata_file = 'gene.csv'
        gene_md_temp_path = f'{temp_dir}/{gene_metadata_file}'
        s3_client.download_file(S3_BUCKET_NAME, 
                                gene_metadata_dir+gene_metadata_file, 
                                gene_md_temp_path)
        gene_md_df = pd.read_csv(gene_md_temp_path, index_col='gene_symbol')
    
    return gene_md_df

    
def get_ccf_labels_image(version=CCF_VERSION, specimen=BRAIN_LABEL):
    '''Loads rasterized image volumes of CCFv3 parcellation as 3D numpy array.

    Loads the "resampled CCF" labels, which have been aligned into the 
    MERFISH space/coordinates; should be paired with "reconstructed" coords.
    Voxels are labelled with assigned brain structure parcellation ID #.
    Rasterized voxels are 10 x 10 x 10 micrometers. First (x) axis is 
    axis is left-right, the second (y) axis is superior-inferior 
    (dorsal-ventral) and third (z) anterior-posterior.
    '''
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))
        # Load image volumes using nibabel
        ccf_labels_dir = f'image_volumes/MERFISH-{specimen}-CCF/{version}/'
        ccf_labels_file = 'resampled_annotation.nii.gz'
        ccf_labels_temp_path = f'{temp_dir}/{ccf_labels_file}'
        s3_client.download_file(S3_BUCKET_NAME, 
                                ccf_labels_dir+ccf_labels_file, 
                                ccf_labels_temp_path)
        img = nibabel.load(ccf_labels_temp_path)
        ccf_imdata = np.array(img.dataobj)

    return ccf_imdata

    
def get_ccf_terms_df(version=CCF_VERSION):
    '''Load CCF parcellation term membership metadata as DataFrame directly
    into memory using tempfile.'''
    with tempfile.TemporaryDirectory() as temp_dir:
        # Initialize S3 client with unsigned credentials for public bucket
        s3_client = boto3.client('s3', region_name=S3_REGION_NAME,
                                 config=Config(signature_version=UNSIGNED))
        # Load metadata as pandas DataFrames from CSV file    
        ccf_metadata_dir = f'metadata/Allen-CCF-2020/{version}/'
        ccf_terms_file = 'parcellation_to_parcellation_term_membership.csv'
        ccf_terms_temp_path = f'{temp_dir}/{ccf_terms_file}'
        s3_client.download_file(S3_BUCKET_NAME, 
                                ccf_metadata_dir+ccf_terms_file, 
                                ccf_terms_temp_path)
        ccf_terms_df = pd.read_csv(ccf_terms_temp_path)

    return ccf_terms_df


'''
----------------------------------------------------------------------------
MANUAL LABEL CURATION FUNCTION ---------------------------------------------
----------------------------------------------------------------------------
'''
def add_label_confidence_to_obs(adata):
    adata = adata.copy()
    '''Add confidence boolean flag for CCFv3/ARA substructure annotations.
    Flags only cells that have good enough alignment to their ARA structures
    such that they can be used to assess method performance.
    CCFv3 (Common Coordinate Framework); ARA (Allen Reference Atlas)

    '''
    # cell_md_df = adata.obs.copy()
    # manual mapping of sections and CCFv3 substructures pairs that have 
    # good alignment
    sec_ccf_dict = {
            6.4 : ['LGv', 'LP', 'PF', 'PO', 'RT', 'VPL', 'VPM', 'ZI-unassigned'],
            6.8 : ['LGv', 'LH', 'MH', 'PO', 'RT', 'VM','VPL', 'VPM', 
                   'ZI-unassigned'],
            8.0 : ['AD','AMd','AMv','AV','IAD','PT','RE','RT']
        }

    # only keep annotations for good section, substructure pairs
    def keep_only_good_annotations(row):
        section = row['z_reconstructed']
        annotation = row['parcellation_substructure']
        
        if section in sec_ccf_dict and annotation in sec_ccf_dict[section]:
            return True
        else:
            return False

     # Apply the function to create the new 'label_curated' column
    adata.obs['label_confidence'] = False
    confidence_bool = adata.obs.apply(keep_only_good_annotations, axis=1)
    adata.obs['label_confidence'] = confidence_bool

    return adata


'''
--------------------------------------------------------------------------------
SPACEHACK2.0 FORMAT FUNCTIONS --------------------------------------------------
--------------------------------------------------------------------------------
'''
def split_adata_into_components(adata):
    features_df = adata.var.loc[:,:]
    observations_df = adata.obs.loc[
                        :, ['brain_section_label', 'average_correlation_score', 
                            'class', 'cluster', 'cluster_alias', 'neurotransmitter',
                            'parcellation_division', 'parcellation_structure',
                            'parcellation_substructure', 'subclass', 'supertype'].copy()
                      ]
    coordinates_df = adata.obs.loc[:, ['x_reconstructed', 
                                       'y_reconstructed', 
                                       'z_reconstructed'].copy()
                                  ]
    coordinates_df.rename(columns={'x_reconstructed':'x', 
                                   'y_reconstructed':'y', 
                                   'z_reconstructed':'z'},
                          inplace=True)
    counts = adata.X.astype('int')  # CSR sparse matrix as dtype-int64
    #TODO Label df should be aligned with what was used for n_clusters
    labels_df = adata.obs[['parcellation_substructure', 'label_confidence']].copy()
    labels_df.rename(columns={'parcellation_substructure':'label'}, 
                     inplace=True)
    
    return coordinates_df, observations_df, features_df, counts, labels_df


def write_sample(path, sample, coordinates_df, observations_df, features_df,
                 counts, labels_df=None, img=None):
    '''Write a sample to disk in the SpaceHack2.0 directory structure.
    '''
    sample_path = Path(path) / sample
    Path.mkdir(sample_path, exist_ok=True)

    coordinates_df.to_csv(sample_path / 'coordinates.tsv', sep='\t', index_label='')
    features_df.to_csv(sample_path / 'features.tsv', sep='\t', index_label='')
    observations_df.to_csv(sample_path / 'observations.tsv', sep='\t', index_label='')
    scipy.io.mmwrite(sample_path / 'counts.mtx', counts)

    if labels_df is not None:
        labels_df.to_csv(sample_path / 'labels.tsv', sep='\t', index_label='')

    if img is not None:
        # TODO write to image_file
        # H_E.json must contain the scale
        pass


def generate_sample_df(adata):
    sections = sorted(adata.obs['z_reconstructed'].unique())

    # set data that's the same across sections/samples
    patient_ls = [BRAIN_LABEL]*len(sections)
    position_ls = [np.nan]*len(sections)
    replicate_ls = [np.nan]*len(sections)

    # generate data that needs to be set per sections/samples
    sample_ls = []
    directory_ls = []
    n_clusters_ls = []
    for sec in sections:
        sec_adata = adata[adata.obs['z_reconstructed']==sec]
        sec_str = str(int(sec))
        
        sample_ls.append(sec_str)
        directory_ls.append(BRAIN_LABEL+'_'+sec_str)

        # number of mapped subclasses should give an accruate starting point
        # for determining the n_clusters to input into various methods, but
        # if this fails to work, alternatives to try would be length of 
        # unique: ['parcellation_substructure', 'supertype']
        n_cl_sec = len(sec_adata.obs['subclass'].unique())
        n_clusters_ls.append(n_cl_sec)

    # put it all together in a df
    sample_df = pd.DataFrame(data={'patient':patient_ls, 
                                    'sample':sample_ls, 
                                    'position':position_ls, 
                                    'replicate':replicate_ls,
                                    'directory':directory_ls, 
                                    'n_clusters':n_clusters_ls
                                   })
    return sample_df


'''
--------------------------------------------------------------------------------
EXECUTED CODE ------------------------------------------------------------------
--------------------------------------------------------------------------------
'''
if __name__=='__main__':
    # Parse output directory from user -----------------------------------------
    parser = argparse.ArgumentParser(
        description="Load data for ABC Atlas Mouse Brain - Thalamus dataset"
    )
    parser.add_argument(
        "-o", "--out_dir", help="Output directory to write files to.", required=True
    )
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    Path.mkdir(out_dir, exist_ok=True)
    
    # Generate dataset ---------------------------------------------------------
    # load ABC Atlas WMB, thalamus subset, dataset
    adata_abc = load_th_subset_adata(counts_transform='raw') # can also load 'log2' counts

    # add label confidence for parcellation substructures 
    adata_abc = add_label_confidence_to_obs(adata_abc)

    # convert xyz coordinates from mm to um
    factor = 1000
    adata_abc.obs[['x_reconstructed', 
                   'y_reconstructed', 
                   'z_reconstructed']] = adata_abc.obs[['x_reconstructed', 
                                                        'y_reconstructed', 
                                                        'z_reconstructed']]*factor

    # no H&E images included with this dataset
    # img = None  # optional

    # Write out each z section as a separate sample -----------------------------
    sections = sorted(adata_abc.obs['z_reconstructed'].unique())
    for sec in sections:
        sec_adata = adata_abc[adata_abc.obs['z_reconstructed']==sec]
        
        # split adata into SpaceHack data structures
        (coordinates_df, 
         observations_df, 
         features_df, 
         counts, 
         labels_df) = split_adata_into_components(sec_adata)
        
        # write section to file
        sec_str = str(int(sec))
        write_sample(out_dir, BRAIN_LABEL+'_'+sec_str, coordinates_df, 
                     observations_df, features_df, counts, labels_df=labels_df)

    # Generate & write experiment metadata ---------------------------------
    # samples.tsv
    samples_df = generate_sample_df(adata_abc)
    samples_df.to_csv(out_dir / "samples.tsv", sep='\t', index_label='')

    # experiment.json
    technology = 'MERFISH'
    species = 'mouse'
    is_3D = True
    with open(out_dir / 'experiment.json', 'w') as f:
        exp_info = {'technology':technology,
                    'species':species,
                    'is_3D':is_3D}
        json.dump(exp_info, f)

    # license file
    license_text = (
        'ABC Atlas - Mouse Whole Brain by Allen Institute for Brain Science'+'\n'+
        ' '+'\n'+
        'ABC Atlas - Mouse Whole Brain'+'\n'+ 
        '(https://knowledge.brain-map.org/data/LVDBJAW8BI5YSS1QUBG/collections)'+'\n'+ 
        'MERSCOPE v1 whole brain Data Collection is licensed under a'+'\n'+ 
        'Creative Commons Attribution 4.0 International License, and'+'\n'+ 
        '10x scRNAseq whole brain Data Collection is licensed under a'+'\n'+ 
        'Creative Commons Attribution-NonCommercial 4.0 International License.'+'\n'+
        ' '+'\n'+
        'See https://alleninstitute.org/citation-policy/ for the Allen Institute Citation'+'\n'+ 
        'Policy and https://alleninstitute.org/terms-of-use/ for the Allen Institute'+'\n'+ 
        'Terms of Use.'+'\n'+
        ' '+'\n'+
        'See https://creativecommons.org/licenses/by/4.0/ and'+'\n'+
        'https://creativecommons.org/licenses/by-nc/4.0/ for a copy of each license.'
        )
    with open(out_dir / 'LICENSE.txt', "w") as f:
        f.write(license_text)