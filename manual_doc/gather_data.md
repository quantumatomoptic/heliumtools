# GatherData

This program is located in misc. It has quite a variety of functions that are very useful.

### Fit and Data Retrieval Functions

#### export_data_set_to_pickle
This function takes as input a sequence folder and a ROI as mandatory arguments. It loads all .atoms files of the corresponding folder and saves them into a single pickle file. If find_arrival_times is True, it will also fit arrival time of each cycle.
The function loops over all .atoms files in the folder.

#### fit_BEC_arrival_time
This function takes as input a dataframe with X, Y, and T + the path to the corresponding *.atoms* file.

#### export_metadatas_to_pickle
Export all metadatas of a folder to a pickle file. Returns a dataframe with metadatas in the metadata_list.\
Takes as arguments:\
``folder: str``\
``path of the folder from which you want to load metadatas.``\
``metadata_list: list, optional``\
``list of strings with metadatas you want to load. By default ["json"]``\
		Possibilities are 
* 'json' for saved parameters,
* 'pico' for picoscope treated datas,
* 'HAL fits' for HAL fits type datas,
* 'MCP stats' for statistics of the MCP HAL extenstion,
* 'all parameters' to load all the sequence parameters 

The function will return a dataframe with the following metadata if one element of the list is contained in the string given in front of each metadata type. This is defined in the *load_metadata* function.

* **Default HAL Parameter:** *".json parameters"* → for example 'param' will return parameters saved in the .json file
* **Picoscope**: *".picoscope_treated"* → for example, 'pico' will return picoscope3000 data
* **ALL parameters**: *"all parameters every parameters"* → for example, 'all param' will return parameters from the total sequence_parameter dictionary. Note that 'param' will not return ALL parameters but default HAL parameters since this option is checked first
* **HAL MCP stats:** *".mcp stats .mcp_stats"* → for example "MCP" will return parameters saved in .MCPstats folder
