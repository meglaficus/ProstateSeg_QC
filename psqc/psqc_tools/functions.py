import cc3d
import SimpleITK as sitk
import numpy as np
from copy import deepcopy
from psqc_tools.scan_class import Scan

def find_sort_components(array: np.ndarray, connectivity: int = 6) -> list:
    """finds all connected components in the array and sorts them by voxel number.

    Args:
        array (np.ndarray): mask in array format
        connectivity (int, optional): Define the connectivity directions. Defaults to 26.

    Returns:
        list: list of tuples (voxel_number, array)
    """
    labels_out, N = cc3d.largest_k(
        array, k=100,
        connectivity=connectivity, delta=0,
        return_N=True,
    )

    arrays = []

    for i in range(1, N + 1):
        new_array = np.zeros(array.shape)
        new_array[labels_out == i] = 1
        total = np.sum(new_array)
        arrays.append((total, deepcopy(new_array)))

    arrays.sort(reverse=True, key=lambda x: x[0])

    return arrays


def filter_small_components(base_array: np.ndarray):
    """Finds all connected components in the array and filters out those that are smaller than 1/10 of the largest component.

    Args:
        base_array (np.ndarray): mask array

    Returns:
        tuple[np.ndarray, bool]: tuple (filtered_array, was_anything_changed)
    """
    arrays = find_sort_components(base_array)

    if len(arrays) > 1:
        biggest = arrays[0][0]

        combined_array = np.zeros(base_array.shape)

        for size, array in arrays:
            if size > biggest / 10:
                combined_array += array

        return combined_array, True

    else:
        return base_array, False


def patch_holes(base_array: np.ndarray) -> tuple[np.ndarray, bool]:
    """Patches holes in the connected components.

    Args:
        base_array (np.ndarray): mask array

    Returns:
        tuple[np.ndarray, bool]: tuple (patched_array, was_anything_changed)
    """
    reversed_array = np.zeros(base_array.shape)
    reversed_array[base_array == 0] = 1

    arrays = find_sort_components(reversed_array)

    if len(arrays) > 1:
        biggest = arrays[0][0]

        combined_array = np.zeros(base_array.shape)

        for size, array in arrays:
            if size < biggest / 100:
                combined_array += array

        patched_array = base_array + combined_array

        return patched_array, True

    else:
        return base_array, False


def process_scan(scan: Scan, to_patch_holes: bool = True, to_filter_small_components: bool = True) -> Scan:
    """The first filters out the small components and then patches the holes in the mask.

    Args:
        scan (Scan): scan object
        to_patch_holes (bool, optional): to patch small holes in mask or not. Defaults to True.
        to_filter_small_components (bool, optional): to filter out small components or not. Defaults to True.

    Returns:
        Scan: processed scan object
    """
    aug_scan = Scan
    base_array = scan.array

    whole = np.zeros(base_array.shape)
    whole[scan.array > 0] = 1

    if to_filter_small_components:
        filtered_array, was_changed = filter_small_components(whole)
    else:
        filtered_array = whole

    if to_patch_holes:
        filtered_array, was_patched = patch_holes(filtered_array)

    else:
        was_patched = False

    if was_changed or was_patched:
        aug_scan = Scan(array=filtered_array, ref=scan.image)
    else:
        aug_scan = scan

    aug_scan.filtered = was_changed
    aug_scan.patched = was_patched

    return aug_scan
