# Dataset class
The idea of this class is to store in a file all the parameters required to perform correlation codes without executing the entire notebook in a well (or not that well) defined order. The example notebook is available [here](../heliumtools/exemples/dataset.ipynb).
## Initialisation

When instanciating the class, you must gives it a valid directory or at least a directory in which it can creates a folder that will be the data location. 
All parameters of your dataset are then stored in a nice readble `yaml` file called 'properties.yml' that is readen when you instanciate the class. 
Important parameters are:
* `name` : it is *where* your data are stored. 
* `sequences`: the list of sequences that are stored in the dataset
* `filters` : filter you want to apply before loading data and metadatas
* `raw_ROI` : one has to be carrefull as we have always a thousand definition for the ROIs. Here, the roi hat is used to choose atoms in (X, Y, T) is defined as raw_ROI. The attribut ROI is kept for the ROI defined in momentum space. 

## Methods
`set(**kwargs)` is the method to set an argument. \
For example, `Dataset.set(hello = 4)` will add the property hello to the dataset with value 4. It is equivalent do `Dataset.hello = 4` BUT each times the `set` function is called, it saves in the YAML file the new parameter and checks also fi the YAML file was updated since the last use.

`add_filter(key, val)` is a method that can be used to add a filter to the filters dictionary. You can also use `set(filters={"Vz":[-90,90]})`  but this will erase previous filters.


`load(filtered = True)` returns you the data and matadatas of the dataset. The optionnal argument *filtered* ables you to apply (or not) the filters you previously defined. It is True by default. 


`get_dataset_properties` show you all your dataset entries

`save_parameters` save your class parameters in the yaml file. Note that this function is called by the `set` function to ensure we do not loose parameters.


`load_data` or `load_metadata` ables you to load (without filers) datas or metadatas.


## Description of sepcific methods 
### Add the sequences to the dataset
Once a sequence is added to the dataset, there is no need to re-add it. 

One can add a sequence to the dataset using the `add_sequence_to_dataset` method. This method use the usual `export_dataset_to_pickle` function, passing it the following arguments:
```
export_data_set_to_pickle(
            folder=seq_dir,
            ROI=self.raw_ROI,
            ROD=self.raw_ROD,
            find_arrival_times=self.fit_arrival_times,
            n_max_cycles=1e8,
            histogramm_width=self.fit_histogram_width,
            ROI_for_fit=self.fit_roi,
            width_saturation=self.fit_width_saturation,
            supplementary_rois=self.supplementary_rois,
            metadata=self.metadata_to_gather,
        )
```

**Reminder of the export_data_set_to_pickle arguments**

* folder: path
    *    Path to the folder containing all the .atoms files.
* ROI: dictionary
    * Region of Interest: we will only select atoms within this ROI. It can be empty, and in that case, we keep all atoms (and not none). Example: {"T": {"min": 300, "max": 350}}. The format of this dictionary must match the official format of an ROI (see the apply_ROI function for more details).
* ROD: dictionary
    * Region of Disinterest: we exclude all atoms within this region. WARNING: this function currently acts axis by axis.
* find_arrival_times: boolean
    * If True, fit the arrival time of the BEC. If False, do not perform the fit.
* n_max_cycles: int
    * Maximum number of cycles to select. Not often useful.
* histogram_width: float,
    * In ms, width of the bins in the histogram for fitting according to T and .times.
* width_saturation: float,
    * Saturation width during which there is no signal due to TDC saturation. Histogram points between tmax, the time at which the signal is maximal, and tmax + dt are removed and not considered in the fit.
* supplementary_rois: list of ROIs
    * Additional ROIs in which we want to count atoms (to apply filters during sequence analysis, for example). Data from these ROIs will be added to the datafbec arrival_times.
* metadata: list of string
    * List of metadata to load with the sequence parameters.