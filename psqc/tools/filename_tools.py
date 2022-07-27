import regex as re
import os


def find_seq_num(scan_name, number_of_digits=4, ignore_miss=False) -> str:
    """ Finds 4 digit id-number of scan in the filename.

    Args:
        scan_name (_type_): filename to search in
        number_of_digits (int, optional): zfill value to search for. Defaults to 4.
        ignore_miss (bool, optional): Weather to ignore error when things go wrong. Defaults to False.

    Returns:
        str: finds the sequence number of the scan (does not correct if it is not 4 digits)
    """
    if number_of_digits == 2:
        pattern = re.compile(r'\d{2,}')
    elif number_of_digits == 3:
        pattern = re.compile(r'\d{3,}')
    elif number_of_digits == 4:
        pattern = re.compile(r'\d{4,}')
    else:
        raise Exception('Invalid number of digits chosen, pick 2, 3 or 4')

    result = re.findall(pattern, scan_name)
    if len(result) == 1:
        return result[0]
    elif len(result) == 2 and result[1] == '0000':
        return result[0]
    elif len(result) == 2 and result[1] == result[0]:
        return result[0]
    elif len(result) > 1 and result[1] != '0000':
        raise Exception(
            f'Found multiple matches for scan {scan_name}')
    elif ignore_miss:
        return None
    elif len(result) == 0:
        raise Exception(f'No matches found! for {scan_name}')


def find_scan_name(pt_id: str, directory: str, step: bool = True) -> str:
    """ Searches the directory for a filename with the given id.

    Args:
        pt_id (str): pt_id to find
        directory (str): where to look
        step (bool, optional): If it finds 0 matches should it revert to searching for zfill(3) strings. Defaults to True.

    Returns:
        str: filename with the given id
    """

    list_of_names = [i for i in os.listdir(directory)
                     if i.endswith(('.nii.gz', '.nii', '.mhd'))]

    matching_scan_names = [
        i for i in list_of_names if find_seq_num(i, ignore_miss=True) == pt_id]

    if len(matching_scan_names) == 1:
        return matching_scan_names[0]

    elif step:
        matching_scan_names = [
            i for i in list_of_names if find_seq_num(i, number_of_digits=3, ignore_miss=True) == str(int(pt_id)).zfill(3)]
        if len(matching_scan_names) == 1:
            return matching_scan_names[0]
        else:
            matching_scan_names = [
                i for i in list_of_names if find_seq_num(i, number_of_digits=2, ignore_miss=True) == str(int(pt_id)).zfill(3)]
        if len(matching_scan_names) == 1:
            return matching_scan_names[0]
        else:
            raise(Exception(
                f'Error: {len(matching_scan_names)} scans found for patient {pt_id}'))

    else:
        raise(Exception(
            f'Error: {len(matching_scan_names)} scans found for patient {pt_id}'))
