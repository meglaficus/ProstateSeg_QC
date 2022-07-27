import os
import SimpleITK as sitk
from tqdm import tqdm
import numpy as np
from ..tools.filename_tools import find_seq_num, find_scan_name

""" Works to combine the two separate files for peripheral zone mask and central zone mask in the italian label-set.
    To use if originally had joined masks. Makes pz 1 and tz 2!"""


def join_masks(peripheral_dir: str, central_dir: str, out_dir: str) -> None:
    os.makedirs(out_dir, exist_ok=True)

    print('Joining masks...')
    for mask in tqdm(os.listdir(peripheral_dir)):
        if mask.endswith('.nii.gz'):
            try:
                pt_id = find_seq_num(mask)
            except:
                pt_id = find_seq_num(mask, number_of_digits=3).zfill(4)

            perif_mask_path = os.path.join(peripheral_dir, mask)
            perif_mask_img = sitk.ReadImage(perif_mask_path)
            perif_mask_arr = sitk.GetArrayFromImage(perif_mask_img)

            central_mask_name = find_scan_name(pt_id, central_dir)
            central_mask_path = os.path.join(
                central_dir, central_mask_name)
            central_mask_img = sitk.ReadImage(central_mask_path)
            central_mask_array = sitk.GetArrayFromImage(central_mask_img)

            if perif_mask_arr.shape != central_mask_array.shape:
                raise Exception(
                    f'Shape of {mask} and {central_mask_name} do not match')

            combi_array = np.zeros(central_mask_array.shape)
            combi_array[perif_mask_arr == 1] = 1
            combi_array[central_mask_array == 1] = 2

            combi_mask_img = sitk.GetImageFromArray(combi_array)
            combi_mask_img.CopyInformation(central_mask_img)

            sitk.WriteImage(combi_mask_img, os.path.join(
                out_dir, f'{mask}'))
