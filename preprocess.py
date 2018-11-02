import cv2
import numpy as np
from PIL import Image
from PIL import ImageEnhance
from glob import glob
import os

def preprocess(imgs):
    rslt = []
    
    for im in imgs:
        if im is None:
            return

        Im = Image.fromarray(im)
        Im = ImageEnhance.Contrast(Im)
        Im = Im.enhance(1.5)
        im = np.array(Im)
        rslt.append(im)
    return rslt
