import SimpleITK as sitk
import numpy as np


class Scan:

    def __init__(self, path=None, array=None, ref=None, image=None):
        if path:
            self.image = sitk.ReadImage(path)

        elif image:
            self.image = image

        if type(array) == np.ndarray:
            self.array = array
            self.image = sitk.GetImageFromArray(self.array)
            self.image.CopyInformation(ref)
        else:
            self.array = sitk.GetArrayFromImage(self.image)

    def write_image(self, path):
        sitk.WriteImage(self.image, path)
