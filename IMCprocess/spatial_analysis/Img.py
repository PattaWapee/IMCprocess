
import pandas as pd
import numpy as np
from skimage import io, color

class Img:
    def __init__(self, img_file):
        self.filename = img_file
        self.load_image()
        self.cancer_mask = None
        self.tissue_mask = None
        
    def load_image(self):
        self.pixels = io.imread(self.filename)
    
    def get_dimensions(self):
        return self.pixels.shape
    
    def display(self):
        io.imshow(self.pixels)
        io.show()

    def load_cancer_mask(self, filename):
        self.cancer_mask = io.imread(filename)
        
    def load_tissue_mask(self, filename):
        self.tissue_mask = io.imread(filename)

    def analyze_cancer_mask(self):
        # Analyze cancer mask and return statistics
        # ...

    def analyze_tissue_mask(self):
        # Analyze tissue mask and return statistics
        # ...
