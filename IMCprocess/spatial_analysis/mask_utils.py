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



class Mask:
    def __init__(self, mask_file, img_id, mask_type):
        '''
        
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
        mask_labeled = label(self.pixels)
        mask_regprops = regionprops(mask_labeled)
        mask_df = table_region(mask_regprops)
        return mask_df



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
