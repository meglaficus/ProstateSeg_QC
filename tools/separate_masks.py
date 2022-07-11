import SimpleITK as sitk
import os
import numpy as np
from tqdm import tqdm

""" Separates the masks if they were originaly joined. This is needed in the for the main function"""

ORIG = ''

OUT = ''

for file_name in tqdm([i for i in os.listdir(ORIG) if i.endswith(('.mhd', '.nii.gz', '.nii'))]):
    path = os.path.join(ORIG, file_name)
    img = sitk.ReadImage(path)
    array = sitk.GetArrayFromImage(img)

    whole_array = np.zeros(array.shape)
    whole_array[array > 0] = 1
    whole_img = sitk.GetImageFromArray(whole_array)
    whole_img.CopyInformation(img)

    perif_array = np.zeros(array.shape)
    perif_array[array == 1] = 1
    perif_img = sitk.GetImageFromArray(perif_array)
    perif_img.CopyInformation(img)

    central_array = np.zeros(array.shape)
    central_array[array == 2] = 1
    central_img = sitk.GetImageFromArray(central_array)
    central_img.CopyInformation(img)

    os.makedirs(os.path.join(OUT, 'whole'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'peripheral'), exist_ok=True)
    os.makedirs(os.path.join(OUT, 'central'), exist_ok=True)

    sitk.WriteImage(whole_img, os.path.join(
        OUT, f'whole/{file_name}'))
    sitk.WriteImage(perif_img, os.path.join(
        OUT, f'peripheral/{file_name}'))
    sitk.WriteImage(central_img, os.path.join(
        OUT, f'central/{file_name}'))
