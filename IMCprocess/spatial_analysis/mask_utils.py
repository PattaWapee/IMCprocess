#from cellpose import io
import cv2
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skimage.measure import label, regionprops
from skimage import segmentation, color
import matplotlib.image as mpimg
from skimage.transform import resize
from skimage import io
import itertools



class Mask:
    def __init__(self, mask_file, img_id, mask_type):
        '''
        Class Name: Mask

        Description: This class represents a mask of an image. 
        It contains methods to load a mask, create a pandas dataframe from the mask

        Attributes:

        filename: A string representing the file name of the mask.
        img_id: A string representing the image ID of the mask.
        mask_type: A string representing the mask type (cancer_mask, tissue_mask, cell_mask).
        mask_df: A pandas dataframe representing the mask.
        Methods:

        init(self, mask_file, img_id, mask_type): 
            Constructor method that initializes the class attributes and loads the mask.
        load_mask(self): Method that loads the mask using skimage's io.imread method 
            and stores it in the pixels attribute.
        create_mask_df(self): Method that creates a pandas dataframe 
            from the mask using skimage's label and regionprops methods and 
            stores it in the mask_df attribute.
        
        '''
        self.filename = mask_file
        self.load_mask()
        self.img_id = img_id
        self.mask_id = mask_type
        self.mask_df = self.create_mask_df()

    def load_mask(self):
        self.pixels = io.imread(self.filename)

    def create_mask_df(self):
        """
        Create a pandas dataframe from the mask
        """
        self.mask_labeled = label(self.pixels)
        self.mask_regprops = regionprops(self.mask_labeled)
        mask_df = table_region(self.mask_regprops)
        return mask_df


class Img_mask:
    def __init__(self, Img_anndata, cell_mask=None, cancer_mask=None, tissue_mask=None):
        """
        Parameters:
        ___________
        cell_mask: a Mask object representing the cell mask
        cancer_mask: a Mask object representing the cancer mask
        tissue_mask: a Mask object representing the tissue mask
        """
        self.img_mask_adata = Img_anndata
        self.cell_mask = cell_mask
        self.cancer_mask = cancer_mask
        self.tissue_mask = tissue_mask
        self.analyze_cancer_mask()
        self.analyze_tissue_mask()
        self.find_stroma_cells()
        self.create_cells_table()
        
    
    def plot_outline_mask(self, masktype, output_file=None, overlay=None):
        if masktype == 'cancer':
            mask_outline = plt_outline(self.cancer_mask.pixels, line_color=(1, 0, 0), 
                                       mode='inner', overlay=overlay, 
                                       output_file=output_file)
        elif masktype == 'tissue':
            mask_outline = plt_outline(self.tissue_mask.pixels, line_color=(1, 0, 0), 
                                       mode='inner', overlay=overlay, 
                                       output_file=output_file)
        elif masktype=='cell_mask':
            mask_outline = plt_outline(self.cell_mask.pixels, line_color=(1, 0, 0), 
                                       mode='inner', overlay=overlay, 
                                       output_file=output_file)
    

    def analyze_cancer_mask(self):
        # Find cell i located in cancer mask and return table
        self.cells_in_cancer_tb, self.cells_in_cancer = create_cell_in_region_table(
            self.cell_mask, self.cancer_mask)

    def analyze_tissue_mask(self):
        # Find cell i located in tissue mask and return table
        self.cells_in_tissue_tb, self.cells_in_tissue = create_cell_in_region_table(
            self.cell_mask, self.tissue_mask)

    def find_stroma_cells(self):
        # Find cell i located in tissue mask and not in cancer mask and return table
        all_cells_in_cancer = list(itertools.chain.from_iterable(
            self.cells_in_cancer.values()))
        all_cells_in_tissue = list(itertools.chain.from_iterable(
            self.cells_in_tissue.values()))
        self.cells_in_both_cancer_tissue = list(set(all_cells_in_cancer).intersection(
        set(all_cells_in_tissue)))
        self.cells_in_stroma = list(set(all_cells_in_tissue).difference(self.cells_in_both_cancer_tissue))

    def create_cells_table(self):
        self.cells_list_df = pd.DataFrame({'cells_in_stroma': [self.cells_in_stroma],
                          'cells_in_cancer': [self.cells_in_cancer],
                          'cells_in_cancer&tissue': [self.cells_in_both_cancer_tissue],
                          'num_cells_in_stroma': [len(self.cells_in_stroma)],
                          'num_cells_in_cancer': [len(self.cells_in_cancer)],
                          'num_cells_in_cancer&tissue': [len(self.cells_in_both_cancer_tissue)]
                          })
        cells_data = self.img_mask_adata.obs.copy().reset_index()
        cancer_indices = np.array(list(itertools.chain.from_iterable(
            self.cells_list_df.loc[0,'cells_in_cancer'].values())))-1
        stroma_indices = np.array(self.cells_list_df.loc[0,'cells_in_stroma'])-1
        cells_data['location'] = ''
        cells_data['cell_label'] = np.array(cells_data.index)+1
        cells_data.loc[cancer_indices,'location'] = 'cancer'
        cells_data.loc[stroma_indices,'location'] = 'stroma'
        self.cells_data = cells_data

def plt_outline(mask, line_color=(1, 0, 0),
                mode='inner', overlay=None,
                output_file=None
                ):
    """
    Plot the outline of a mask on a black background (overlay=None) 
    or on a specific background (overlay=tiff_filename).

    Parameters:
    ___________
    mask: numpy ndarray of the mask 
    line_color: tuple of the color of the outline
    mode: mode parameter in segmentation.mark_boundaries
    overlay: tiff file to overlay the mask on (optional)
    return:
    _______
    outline_image: numpy ndarray of the outline image

    """
    outline_image = segmentation.mark_boundaries(color.gray2rgb(mask),
                                                 mask,
                                                 color=line_color,
                                                 mode=mode)
    if overlay is not None:
        # Resize the overlay image to match the size of the outline image
        overlay_image_resized = resize(overlay, mask.shape, anti_aliasing=True)
        # Convert the grayscale overlay image to RGB
        overlay_image_resized = color.gray2rgb(overlay_image_resized)
        # Blend the overlay image and the outline image
        alpha = 0.5  # Set the blending ratio
        outline_image = alpha * outline_image + (1 - alpha) * overlay_image_resized

    plt.imshow(outline_image)
    # Convert the depth of the image to uint8
    outline_image_uint8 = np.array(outline_image * 255, dtype=np.uint8)

    if output_file is not None:
        # Save the outline image to a file
        mpimg.imsave(output_file, outline_image_uint8)
    return outline_image


def plt_outline_label(mask, region_props,
                      labeled_cells,
                      line_color=(1, 0, 0),
                      mode='inner',
                      output_file=None):
    """
    Plot the mask on a black background with the label of each region

    Parameters:
    ___________
    mask: numpy ndarray of the mask 
    region_props: region_probs from skimage.measure.regionprops
    return:
    _______
    mask_image: numpy ndarray of the mask image with labels

    """
    # Generate an outline image using the labeled mask
    outline_image = segmentation.mark_boundaries(
        color.gray2rgb(mask),
        labeled_cells,
        color=line_color,
        mode=mode)

    # Iterate over each labeled region and draw the label number on the outline image
    for region in region_props:
        # Compute the centroid of the region
        cy, cx = region.centroid

        # Draw the label number on the outline image
        cv2.putText(outline_image, str(region.label), (int(cx), int(
            cy)), cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 0, 1), 2)

    plt.imshow(outline_image)
    plt.axis('off')
    plt.show()
    # Convert the depth of the image to uint8
    outline_image_uint8 = np.array(outline_image * 255, dtype=np.uint8)

    if output_file is not None:
        # Save the outline image to a file
        mpimg.imsave(output_file, outline_image_uint8)
    return outline_image


def table_region(region_props):
    """
    Create a table of the region properties of the mask.
    Parameters:
    ___________
    region_props: region_probs from skimage.measure.regionprops 
    return:
    _______
    table: pandas dataframe of the region properties
    """
    labels = [reg.label for reg in region_props]
    areas = [reg.area for reg in region_props]
    centx = [reg.centroid[1] for reg in region_props]
    centy = [reg.centroid[0] for reg in region_props]

    reg_table = pd.DataFrame({'label': labels,
                              'area': areas,
                              'centroid_x': centx,
                              'centroid_y': centy})
    return reg_table


def cell_in_region(cell_regprops, mask_regprops):
    '''
    Check if the centroid of the cell is inside the bounding box of an item in region from mask.
    Parameters:
    ___________
    cell_regprops: cell region_probs from skimage.measure.regionprops 
    mask_regprops: other type (such as cancer) region_probs from skimage.measure.regionprops 
    return:
    _______
    cell_in_cancer: dictionary of the list of cell labels inside the mask regions
    '''

    cell_in_region = {}
    [cell_in_region.setdefault(region_j.label, [])
     for region_j in mask_regprops]
    cell_outside_region = []

    for i, cell_i in enumerate(cell_regprops):
        #print(f"Checking cell {i+1}...")
        celli_centroid = cell_i.centroid
        for j, region_j in enumerate(mask_regprops):
            #print(f'Checking cell {i+1} against region {j+1}...')
            #cell_in_region.setdefault(region_j.label, [])
            if region_j.bbox[0] <= celli_centroid[0] <= region_j.bbox[2] and region_j.bbox[1] <= celli_centroid[1] <= region_j.bbox[3]:
                #print(f"Cell {i+1} is located inside an item in cancer region {j+1}.")
                cell_in_region[region_j.label].append(cell_i.label)
                break

        # The else block will be executed if the loop completes without encountering a break statement,
        # indicating that the cell is not located inside any of the cancer regions.
        else:
            cell_outside_region.append(cell_i.label)
            #print(f"Cell {i+1} is located outside all items in cancer region.")
    return cell_in_region, cell_outside_region

def create_cell_in_region_table(Cell_mask, Region_mask):
    cells_in_Reg, cells_outside_Reg = cell_in_region(Cell_mask.mask_regprops, 
                                                        Region_mask.mask_regprops)
    # add number of cells in each region to the table
    cell_in_reg_tb =Region_mask.mask_df.copy()
    cell_in_reg_tb['num_cells'] = cell_in_reg_tb['label'].map(lambda x: len(cells_in_Reg[x]))
    cell_in_reg_tb['cells_in_region'] = cell_in_reg_tb['label'].map(
        lambda x: cells_in_Reg[x])
    return cell_in_reg_tb, cells_in_Reg
    


#################################################
# These function below is for the analysis of ###
# the cell type in the stroma or cancer regions #
#################################################


def get_celltype_fraction(obs_object, col_label, cell_in_reg_df, col_cell_type):
    obs_df = obs_object.copy()
    cellid = cell_in_reg_df.loc[0, col_cell_type]
    if isinstance(cellid, dict):
        cellid_list = []  # The target list to merge all the values into

        # Iterate through the values of the dictionary
        for lst in cellid.values():
            cellid_list.extend(lst)
    else:
        cellid_list = cellid

    # add label column to obs for mapping with cell_data
    obs_df['label'] = np.array(range(1, len(obs_df)+1))
    #How many cells for each cell types in stroma?
    obs_df[obs_df['label'].isin(cellid_list)][col_label].value_counts()

    fraction_df = pd.DataFrame(obs_df[obs_df.label.isin(
        cellid_list)][col_label].value_counts(
            normalize=True)).rename(columns={col_label:'fraction_'+ col_cell_type})
    return fraction_df

def plt_fraction_df(fraction_df, output_file):
    fraction_df.plot(kind='bar')
    # Add labels to the bars
    for i in range(len(fraction_df)):
        label = round(fraction_df.iloc[:,0][i], 2)  # Round value to one decimal place
        plt.text(i, fraction_df.iloc[:,0][i] + 0.01, label, ha='center')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
