#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the Tue Feb 20 2024 by Victor
#
# Developped by Victor, ...
#
# Last (big) change on the ... by ...
#
# Copyright (c) 2024 - Helium1@LCF
# ----------------------------------
#
"""
Content of dataset.py

Please document your code ;-).

"""

import os, logging, re
from heliumtools.misc.gather_data import export_data_set_to_pickle
import pandas as pd
from heliumtools.misc.logger import getLogger
from heliumtools.tools import data_filter

log = getLogger(__name__)
# set the desired level of warning : DEBUG / INFO / WARNING
log.setLevel(logging.DEBUG)
import yaml


class Dataset:
    def __init__(self, name):
        """Instanciation of a dataset defined by its folder name which is also
        its location. If it does not exists, the function creates the associated
        folder.

        Parameters
        ----------
        name : path
            path to the dataset folder. It might not be the sequence directory.
        """
        self.name = name
        self.sequences = []
        self.raw_ROI = {"T": [312, 328], "X": [-30, 0], "Y": [-15, 15]}

        self.raw_ROD = {}
        self.metadata_to_gather = ["picoscope", "bec", "param"]
        self.supplementary_rois = []
        self.filters = {}
        # for the fit of the
        self.fit_roi = {
            "T": {"min": 306.2, "max": 309.7},
            "X": [-30, 0],
            "Y": [-15, 15],
        }
        self.fit_arrival_times = False
        self.fit_histogram_width = 0.01
        self.fit_width_saturation = 0  # saturation for the fit of the BEC
        self._check_dataset_location()
        self.load_parameters()
        self.save_parameters()

    def set(self, **kwargs):
        self.__dict__.update(**kwargs)
        self.save_parameters()

    def remove(self, *args) -> None:
        """remove an element from the parameters"""
        self.save_parameters()
        for arg in args:
            if arg in self.__dict__ and (arg not in ["sequences", "name"]):
                if arg == "filters":
                    self.filters = {}

                self.__dict__.pop(arg)
        self.save_parameters(load_file_before_dump=False)

    def add_filter(self, key, val):
        self.filters[key] = val
        self.save_parameters()

    def load_parameters(self):
        """method that load the properties of the dataset store in the yaml file."""
        file_path = os.path.join(self.name, "properties.yml")
        try:
            with open(file_path, "r") as file:
                data = yaml.safe_load(file)
                old_name = data.pop("name", None)
                if old_name != self.name:
                    self._raise_path_warning(old_name)
                self.__dict__.update(data)
        except FileNotFoundError:
            log.info(
                "The file parameter does not exist. Creating it with default parameters."
            )
            self.save_parameters(load_file_before_dump=False)
        except Exception as e:
            msg = "An error occured while loading the properties.yml file. Please take a look at it. \n "
            msg += f"[ERROR] --> {e}"

            log.warning(msg)

    def _raise_path_warning(self, old_name) -> None:
        """warning function for the user so that he knows that he might have path issues if he tries to use path registered

        Parameters
        ----------
        old_name : old datset name
            string that is the path that WAS the dataset
        """
        msg = f"The path registered {old_name} does not match your path {self.name}. "
        msg += "This might be due to windows/unix path compatibility issues. The path to the dataset data is well defined however sequence location might be wrong."
        log.warning(msg)
        msg = "Note that this is not a problem if you do not need to update the atoms nor the metadata stored in the dataset folder. However, it can be an issue if you want to re-export datas to the dataset. If so, I advise you to change the sequences attributs to switch it to the real sequences you are aiming to gather. Good luck !"
        log.info(msg)
        log.info("The saved sequences are the following : {}".format(self.sequences))

    def save_parameters(self, load_file_before_dump=True) -> None:
        """Save the parameters of the class into the properties.yml file.
        This method makes sure that a parameter was not manually added or that an other
        file added a parameter. This can be ingored using the load_file_before_dump
        argument.

        Parameters
        ----------
        load_file_before_dump : bool, optional
            if we check the file before saving, by default True
        """
        file_path = os.path.join(self.name, "properties.yml")
        if load_file_before_dump:
            try:
                with open(file_path, "r") as file:
                    data = yaml.safe_load(file)
                for key, val in data.items():
                    if key not in self.__dict__.keys():
                        self.__dict__[key] = val
                        log.info(
                            f"{key} was in the configuration file and not in your dataset properties. Adding it."
                        )
            except Exception as e:
                msg = "Loading the previous properties file "
                msg += "to check variables before write failed."
                msg += f"I will overwrite the file anyway.Error is {e}"
                log.warn(msg)
        try:
            with open(file_path, "w") as file:
                yaml.dump(self.__dict__, file, default_flow_style=False)
        except Exception as e:
            msg = f"Failed to save parameters in {file_path}. Error is {e}"
            log.error(msg)

    def _check_dataset_location(self):
        if not os.path.exists(self.name):
            try:
                os.mkdir(self.name)
                msg = "Dataset was initialized in folder {}".format(self.name)
                log.info(msg)
            except Exception as e:

                msg = "Dataset initialization failed in folder {}. ".format(self.name)
                msg += "Please check the error log. \n " + str(e)
                log.critical(msg)
                raise ValueError(msg)

    def get_dataset_properties(self):
        log.info("Requiring dataset properties...")
        for key, val in self.__dict__.items():
            msg = key.replace("_", " ")
            msg += " " * max(20 - len(msg), 0)
            try:
                str_val = str(val)
                msg += ": {}".format(str_val[0:50])
                print(msg)
            except Exception:
                print(msg + ": output failed.")

    def get_sequence_id_from_seq_dir(self, seq_dir: str) -> str:
        """returns the sequence id string from the sequence

        Parameters
        ----------
        seq_dir : str
            sequence directory string

        Returns
        -------
        str
            _description_
        """
        seq_dir = str(seq_dir)
        seq_id = seq_dir.replace("/", "-").replace("\\", "-")
        match = re.search(r"\b\d{4}-\d{2}-\d{2}-\d{3}\b", seq_id)

        # Vérification si une correspondance a été trouvée
        if match:
            seq_id = match.group(0)
        else:
            msg = f"The sequence ID was not recognized from {seq_id}"
            log.error(msg)
            return None
        return seq_id

    def add_sequence_to_dataset(self, seq_dir, force_update=False) -> None:
        """function that add a sequence to the dataset

        Parameters
        ----------
        seq_dir : string
            path to the sequence to add to dataset
        find_arrival_time : bool, optional
            if one fits arrival time or not, by default False
        force_update : bool, optionnal
            if the sequence is already in the database, force to update
        """

        seq_dir = str(seq_dir)
        seq_id = self.get_sequence_id_from_seq_dir(seq_dir)
        if not seq_id:
            msg = f"Sequence {seq_dir} was not added to the database"
            log.warning(msg)
            return
        if not os.path.exists(seq_dir):
            msg = f"The sequence directory {seq_dir} does not exists"
            log.error(msg)
            return
        if seq_id in self.sequences:
            msg = f"Sequence {seq_dir} is already in the database (registered under {seq_id})."
            if not force_update:
                log.warning(msg)
                return
            else:
                msg += " Removing sequence is not yet implemented.... Sorry."
                log.warning(msg)
                return
        # we add the sequence to the dataset

        [atoms, params, bec_arr] = export_data_set_to_pickle(
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
        new_atoms = pd.read_pickle(os.path.join(seq_dir, "dataset.pkl"))
        old_atoms = self.load_data()
        self._save_data(pd.concat([old_atoms, new_atoms]))
        if bec_arr:
            new_metadata = pd.read_pickle(os.path.join(seq_dir, "arrival_times.pkl"))
        else:
            new_metadata = pd.read_pickle(os.path.join(seq_dir, "parameters.pkl"))

        old_meta = self.load_metadata()
        self._save_metadata(pd.concat([old_meta, new_metadata]))
        self.sequences.append(seq_id)
        self.save_parameters()
        log.info(f"Sequence {seq_dir} was succesfully added to the dataset.")
        #     return
        # except Exception as e:
        #     msg = f"The export of the sequence {seq_dir} failed. Error is {e}."
        #     log.error(msg)
        #     return

    def load_data(self) -> pd.DataFrame:
        """funcion that load data from the dataset.

        Returns
        -------
        pd.DataFrame
            data loaded from the dataset
        """
        data_dir = os.path.join(self.name, "data.pkl")
        if os.path.exists(data_dir):
            try:
                data = pd.read_pickle(data_dir)
                return data
            except Exception as e:
                msg = f"Failed to load data. Error is {e}."
                log.warning(msg)
        return pd.DataFrame({"X": [], "Y": [], "T": [], "Cycle": []})

    def _save_data(self, data):
        """function that save datas to the data archive that contains atoms.

        Parameters
        ----------
        data : pd.DataFrame
            pandas dataframe that contains atoms to be saved
        """
        data.to_pickle(os.path.join(self.name, "data.pkl"))

    def load_metadata(self) -> pd.DataFrame:
        """funcion that load metadata from the dataset.

        Returns
        -------
        pd.DataFrame
            data loaded from the dataset
        """
        data_dir = os.path.join(self.name, "metadata.pkl")
        if os.path.exists(data_dir):
            try:
                metadata = pd.read_pickle(data_dir)
                return metadata
            except Exception as e:
                msg = f"Failed to load metadata. Error is {e}."
                log.warning(msg)
        return pd.DataFrame({"Cycle": []})

    def _save_metadata(self, metadata):
        """function that save metadatas to the data archive that contains metadatas.

        Parameters
        ----------
        metadata : pd.DataFrame
            pandas dataframe that contains metadata to be saved
        """
        metadata.to_pickle(os.path.join(self.name, "metadata.pkl"))

    def clean_up_dataset(self):
        """function that cleans up the dataset."""
        try:
            os.remove(os.path.join(self.name, "metadata.pkl"))
            os.remove(os.path.join(self.name, "data.pkl"))
            self.sequences = []
            self.save_parameters()
            log.warning("I suppressed all data and metadatas of the dataset.")
        except OSError as e:
            print(f"Error while cleaning up the dataset: {e}")

    def load(self, filtered=True):
        """method that loads data and metadata from the dataset

        Parameters
        ----------
        filtered : bool, optional
            if we filter data we the saved filters, by default True

        Returns
        -------
        _type_
            _description_
        """
        data = self.load_data()
        metadata = self.load_metadata()
        if not filtered:
            return data, metadata
        return data_filter(data, metadata, self.filters)


## %% FOR TEST %%
if __name__ == "__main__":
    d = Dataset("/home/victor/Bureau/tmp/salut")
    d.get_dataset_properties()
