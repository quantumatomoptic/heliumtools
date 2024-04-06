# Tools

For some time now, the ROI selection functionalities that we use in all our programs rely on the *apply_ROI* function. This standardization allows to take into account various ROI formats such as lists, dictionaries with min and max, dictionaries with center and size, etc.

## ROI Format
The ROIs managed by the program must be dictionaries where each entry ("Vx", "Vy" for example) can be:

* a dictionary of the form {"min": *float*, "max":*float*},
* a dictionary of the form {"position": *float*, "size":*float*},
* a list or tuple with the minimum and maximum.

To add a new ROI format, you need to modify the *get_roi_min_max(roi, key)* function which is called by all the other functions.
## Available functions
* __apply_ROI(dataframe, ROI)__: applies the ROI to the dataframe. It is necessary that every entry of the ROI is contained in the dataframe.
* __apply_ROD(dataframe, ROI)__: same as below but with a region of disinterest (exclusion).
* __get_roi_min_max(roi, key)__: returns the min and max of a certain key of the ROI.
* __get_roi_size(roi, axis)__: returns the size of a certain key of the ROI.

