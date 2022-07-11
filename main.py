from copy import deepcopy

import pandas as pd
from tools.filename_tools import *
from tqdm import tqdm
import numpy as np

from tools.functions import filter_small_components, process_scan
from tools.scan_class import Scan

# Enter paths to folders containing the scans:
WHOLE_PATH = ''
PERIPHERAL_PATH = ''
CENTRAL_PATH = ''

# Change to save somewhere else:
WHOLE_OUT = 'out/whole'
PERIPHERAL_OUT = 'out/peripheral'
CENTRAL_OUT = 'out/central'


def main(to_save: bool = True, changed_only: bool = True) -> None:
    """Main function of this program. Performs the analysis of all the masks. Requires the paths to the folders containing the masks to be defined.
    Only works if all the folder paths are defined and if they contain all the masks. Currently only works with .nii.gz, .nii and .mhd files.

    Args:
        to_save (bool, optional): to save processed masks or not. Defaults to True.
        changed_only (bool, optional): if to_save whether to save all scans or just modified ones. Defaults to True.
    """

    # Create dataframe to store all changes

    scan_names = sorted(i for i in os.listdir(
        WHOLE_PATH) if i.endswith(('.nii.gz', '.nii', '.mhd')))

    df = pd.DataFrame(columns=['scan_name', 'whole_filtered', 'whole_patched', 'whole_mismatch', 'perif_filtered',
                               'perif_patched', 'central_filtered', 'central_patched', 'strays_converted'])

    for scan_name in tqdm(scan_names):
        # Finds patient id to find matching mask in other folders
        try:
            pt_id = find_seq_num(scan_name)
        except:
            pt_id = find_seq_num(scan_name, number_of_digits=3).zfill(4)

        whole_scan_path = os.path.join(WHOLE_PATH, scan_name)
        whole_scan0 = Scan(path=whole_scan_path)

        central_scan_name = find_scan_name(pt_id, CENTRAL_PATH)
        central_scan_path = os.path.join(CENTRAL_PATH, central_scan_name)
        central_scan = Scan(path=central_scan_path)

        perif_scan_name = find_scan_name(pt_id, PERIPHERAL_PATH)
        perif_scan_path = os.path.join(PERIPHERAL_PATH, perif_scan_name)
        perif_scan = Scan(path=perif_scan_path)

        # Create another 'whole' scan to make all the three masks match up
        whole_scan_array = np.zeros(central_scan.array.shape)
        whole_scan_array[(central_scan.array == 1) |
                         (perif_scan.array == 1)] = 1
        whole_scan = Scan(array=whole_scan_array, ref=central_scan.image)
        mismatch = False
        try:
            whole_test = whole_scan.array - whole_scan0.array
            if np.count_nonzero(whole_test) != 0:
                tqdm.write(
                    f'{scan_name} whole mask does not match the peripheral and central masks! {np.sum(np.nonzero(whole_test))} voxels different.')
                mismatch = True
        except Exception as e:
            tqdm.write(
                f'{e} \n Something wrong with {scan_name}, possibly the dimensions dont match between the different masks')

        # Processes whole prostate mask
        whole_scan_aug = process_scan(
            whole_scan, to_patch_holes=True, to_filter_small_components=True)

        # Processes central zone mask

        central_scan_aug = process_scan(
            central_scan, to_patch_holes=True, to_filter_small_components=True)

        # Finds small components in central zone mask that are included in the processed whole prostate mask
        central_array_base = central_scan.array
        central_array_aug0 = deepcopy(central_array_base)
        central_array_aug0[whole_scan_aug.array == 0] = 0
        filtered_central_array, converted_strays = filter_small_components(
            central_array_aug0)

        # Does initial processing of peripheral zone mask

        perif_scan_aug_raw = process_scan(
            perif_scan, to_patch_holes=True, to_filter_small_components=True)

        # Adds the small components that were in the central zone mask and in the processed whole prostate mask.
        # These are considered 'strays' and are believed to be erroneously included in the central zone mask.
        perif_strays = central_array_aug0 - filtered_central_array
        perif_array_aug = perif_scan_aug_raw.array + perif_strays
        perif_scan_aug = Scan(array=perif_array_aug, ref=perif_scan.image)

        if to_save:
            # If changed_only is True, only write the files if there were changes else writes all files
            if changed_only:
                if whole_scan_aug.filtered or whole_scan_aug.patched or mismatch:
                    whole_scan_aug.write_image(
                        os.path.join(WHOLE_OUT, scan_name))

                if perif_scan_aug_raw.filtered or perif_scan_aug_raw.patched or converted_strays:
                    perif_scan_aug.write_image(
                        os.path.join(PERIPHERAL_OUT, scan_name))

                if central_scan_aug.filtered or central_scan_aug.patched or converted_strays:
                    central_scan_aug.write_image(
                        os.path.join(CENTRAL_OUT, scan_name))

            else:
                whole_scan_aug.write_image(os.path.join(WHOLE_OUT, scan_name))
                central_scan_aug.write_image(
                    os.path.join(CENTRAL_OUT, central_scan_name))
                perif_scan_aug.write_image(
                    os.path.join(PERIPHERAL_OUT, perif_scan_name))

        # Adds row to dataframe, logging all the findings
        df.loc[pt_id] = {'scan_name': scan_name,
                         'whole_filtered': whole_scan_aug.filtered, 'whole_patched': whole_scan_aug.patched, 'whole_mismatch': mismatch,
                         'perif_filtered': perif_scan_aug_raw.filtered, 'perif_patched': perif_scan_aug_raw.patched,
                         'central_filtered': central_scan_aug.filtered, 'central_patched': central_scan_aug.patched,
                         'strays_converted': converted_strays}

    # Saves information about all changes to a csv file
    df.sort_index(inplace=True)
    df.to_csv('change_log/all_mods.csv')

    # Creates new csv files for each zone, including only masks that were changed and logs the changes
    whole_df = df[['scan_name', 'whole_filtered',
                   'whole_patched']]
    whole_df = whole_df[(whole_df['whole_filtered'] == True) | (
        whole_df['whole_patched'] == True)]
    whole_df.to_csv('change_log/whole_mods.csv')

    central_df = df[['scan_name', 'central_filtered',
                     'central_patched', 'strays_converted']]
    central_df = central_df[(central_df['central_filtered'] == True) | (
        central_df['central_patched'] == True) | (central_df['strays_converted'] == True)]
    central_df.to_csv('change_log/central_mods.csv')

    perif_df = df[['scan_name', 'perif_filtered',
                   'perif_patched', 'strays_converted']]
    perif_df = perif_df[(perif_df['perif_filtered'] == True) | (
        perif_df['perif_patched'] == True) | (perif_df['strays_converted'] == True)]
    perif_df.to_csv('change_log/perif_mods.csv')


if __name__ == '__main__':
    os.makedirs(WHOLE_OUT, exist_ok=True)
    os.makedirs(PERIPHERAL_OUT, exist_ok=True)
    os.makedirs(CENTRAL_OUT, exist_ok=True)
    os.makedirs('change_log', exist_ok=True)

    main(to_save=True)
