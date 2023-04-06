from cellpose import io
from skimage.measure import label, regionprops
from skimage import segmentation, color



def plt_outline(mask, color=(1, 0, 0), mode='inner', overlay=None):
    """
    Plot the outline of a mask on a black background (overlay=None) 
    or on a specific background (overlay=tiff_filename).

    Parameters:
    ___________
    mask: numpy ndarray of the mask 
    color: tuple of the color of the outline
    mode: mode parameter in segmentation.mark_boundaries
    overlay: tiff file to overlay the mask on (optional)
    return:
    _______
    outline_image: numpy ndarray of the outline image
    
    """
    outline_image = segmentation.mark_boundaries(color.gray2rgb(mask), 
                                                 mask, 
                                                 color=color, 
                                                 mode=mode)
    io.show(outline_image)
    return outline_image
