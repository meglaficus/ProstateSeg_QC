# ProstateSeg_QC


ProstateSeg_QC is a quality control algorithm to fix prostate zone and lesion segmentation masks. It finds and rectifies common errors that occur during manual segmentation.

## What it does

This is an all-in-one algorithm to process the entire dataset. It is largely based on finding and analyzing connected components.

QC of zone segmentation masks:

1. Finds and removes all small connected components from all the masks.
2. Finds and patches all small holes in the masks.
3. Finds snippets that are labeled as central zone but were very clearly just small errors when marking the peripheral zone mask onto the whole prostate mask. It converts those to peripheral zone.
4. Checks if the whole mask equals the sum of the peripheral and central zone masks. If not, it replaces the whole mask. (optional)

QC of lesion segmentation masks:
1. Finds and removes all small connected components from all the masks.
2. Finds and patches all small holes in the masks.

## Examples


<div class="center"> Filtering small components <div >  <img title="Filtering small components" src="./examples/fragment_before.png" width="250" height="250" /> <img src="./examples/fragment_after.png" width="250" height="250" alt="example 1"/> 
 </div></div>
 
<div class="center"> Patching holes <div >  <img src="./examples/hole_before.png" width="250" height="250" /> <img src="./examples/hole_after.png" width="250" height="250" alt="example 2"/> 
 </div></div>

<div class="center"> Fixing snippets <div >  <img src="./examples/snippet_before.png" width="250" height="250" /> <img src="./examples/snippet_after.png" width="250" height="250" alt="example 3"/> 
</div></div>

## Requirements

The program requires these modules to be installed:

- numpy
- pandas
- SimpleITK
- connected-components-3d
- tqdm
- regex

## Installation
`pip install prostate_seg_qc`


## Usage
```python
from pqc.qc import qc_zone, qc_lesion

# To do quality control on zone masks
qc_zone(COMBINED_PATH='path/to/combined/mask/directory')

# To do quality control on lesion masks
qc_lesion(LESION_PATH='path/to/lesion/mask/directory')
```



By default the program saves only the masks that were changed, saving them into new directories that are created in the cwd. It also creates .csv files that store the information on all the errors that were found.


## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
