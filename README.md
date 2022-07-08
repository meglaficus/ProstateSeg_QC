# ProstateSeg_QC


ProstateSeg_QC is a quality control algorithm to fix prostate zone segmentation masks. It finds and rectifies common errors that occur during manual segmentation.

## Method

This is an all-in-one algorithm to process the entire dataset. It is largely based on finding and analyzing connected components.

The algorithm requires separate files for the whole prostate mask, the peripheral zone mask and the non-peripheral (central) zone mask. It finds and removes all small connected components from all the masks. It also patches all the holes in the masks. 

It also pays special attention to all the snippets that are labeled as central zone but were very clearly just small errors when marking the peripheral zone mask onto the whole prostate mask.


## Usage

First clone or fork repo. Then open main.py and enter the paths to relevant directories at top of file as show bellow. Then just run main.

```python
WHOLE_PATH = '/<example path>/whole'
PERIPHERAL_PATH = '/<example path>/peripheral'
CENTRAL_PATH = '/<example path>/central'
```
By default the program saves only the masks that were changed, saving them into new directories that are created in the cwd. It also creates .csv files that store the information on all the errors that were found.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.
