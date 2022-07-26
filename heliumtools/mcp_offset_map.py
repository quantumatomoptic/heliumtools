#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 22-7-2022 by Victor, using the code of Quentin
#
# Developped by Victor, ...
#
# Last (big) change on the ... by ...
#
# Copyright (c) 2022 - Helium1@LCF
# ----------------------------------
#
"""
Content of mcp_offset_map.py
-----------------------------

Cette classe sert d'interface avec la base de données permettant de sauvegarder les onnées. Nous avons ainsi affaire à deux tables : 
La table sequences qui conbtient toutes les séquences ajoutées à la base de données (chaque atome n'a pas un unique id donc on ne veut pas l'ajouter deux fois)
La table atoms contenant la position et l'offset.

----------- Table atoms -----------


--------- Table sequences ---------



----------- Table pixels ----------


#### DEVELOPPERS REMARKS ######

Follows some usefull commands to deal with sql database

*** Basis ***

self.database_name = /.../offset.db   # le nom de notre base de donnée
self.connexion = sqlite3.connect(self.database_name)   # cela permet de nous connecter 
self.cursor = self.connexion.cursor() # le curseur permettant de nous y connecter.

*** SQL and Pandas ***
pixels = pd.read_sql_query(("SELECT * FROM pixels", self.connexion)) # charge les données de la table 'pixels' de la base de donnée dans un pandas dataframe
dataframe.to_sql(name_of_table, self.connexion, if_exists='fail',index = True) 
# Always use index = False
# For atoms : use if_exists = 'append'   ;because we always want to append the atom (if we already add an atom with this offset, we have to add it and not to drop it).
# For pixels, it should be 'replace' 



Some note (v.g. on the 26 of July) :
We have a lot of speed issue when constructing the database.I tried to delete each variable 

"""

import os, re, json, copy
from pathlib import Path
from select import select
import pandas as pd
import numpy as np
import sqlite3 as db
from tqdm import tqdm


class MCP_offset_map:
    def __init__(self, path_to_folder="mcp_offset_map"):
        self.path_to_folder = path_to_folder
        database_name = "offsets.db"
        result_name = "result.pkl"
        self.database_name = os.path.join(self.path_to_folder, database_name)
        self.result_name = os.path.join(self.path_to_folder, result_name)
        self.time_resolution = 120.0e-12  # seconds
        self.dll_temp_lenght = 80.0e-9  # seconds
        self.mcp_diameter = 80  # mm

        # les séquences qui ont contribué à faire la carte d'offset sous la forme
        # {"2022/08/09/023": number of cycles}
        self.sequences_that_contributed = {}
        self.load_existing_datas()
        self.connect_to_database()

    def connect_to_database(self):
        """This methods connect to the offset database."""
        self.connexion = db.connect(self.database_name)
        self.cursor = self.connexion.cursor()

    def unconnect(self):
        """This method unconnect from the database"""
        self.connexion.close()

    def load_existing_datas(self):
        """
        Cette méthode vérifie que le dossier self.path_to_folder existe bien. Elle vérifie également que celui-ci contient bien les fichiers avec les données pour réaliser la carte d'offset. Si elle continet des fichiers, elle charge
        la carte des offsets moyens
        la carte des offset std
        la carte du nombre de coups ayant réalisé l'offset
        les directories de sséquences ayant contribué à la carte d'offset
        --> on ne charge pas la carte avec la valeur de chaque offset car je pense que ce sera trop gros.
        """
        # 1 : on verifie que le dossier renseigné existe bien
        if not os.path.isdir(self.path_to_folder):
            print("Your MCP folder does not exist : I create it and every needed files")
            os.mkdir(self.path_to_folder)
        ## Then we deal with the database :
        self.connect_to_database()
        # récupère les tables déjà existantes : si elles n'existent pas, on les crée.
        listofTables = self.cursor.execute(
            "SELECT name FROM sqlite_master WHERE type='table';"
        ).fetchall()
        # On récupère une liste de la forme suivante : [("sequences", ), ("atoms", ), ("pixels", )] si les tables existent et une liste vide sinon.
        if ("sequences",) not in listofTables:
            self.cursor.execute("CREATE TABLE sequences (name NVARCHAR(50))")
        if ("atoms",) not in listofTables:
            self.cursor.execute(
                "CREATE TABLE atoms (deltaX INT, deltaY INT, offset INT)"
            )
        if ("pixels",) not in listofTables:
            self.cursor.execute("CREATE TABLE pixels (deltaX INT, deltaY INT)")
        self.connexion.commit()  # il faut commit pour prendre en compte l'ajout
        self.unconnect()

        # On charge maintenant les fichiers pixels et results
        self.load_pixels()
        self.load_results()

    def load_pixels(self):
        self.connect_to_database()
        self.pixels = pd.read_sql_query("SELECT * FROM pixels", self.connexion)
        self.unconnect()

    def load_results(self):
        if os.path.isfile(self.result_name):
            self.result = pd.read_pickle(self.result_name)
        else:
            self.result = copy.deepcopy(self.pixels)
            self.result["mean"] = 0
            self.result["std"] = 0
            self.result["count"] = 0
            self.result["max"] = 0
            self.result["min"] = 0

    def save_result(self):
        """Save the result dataframe as pkl"""
        self.result.to_pickle(self.result_name)

    ##################################################################
    ########## METHDS TO ADD DATAS TO DATABASE
    ##################################################################
    def add_sequence_to_database(self, sequence_directory):
        """Ajoute la séquence sequence_directory à la base de donnée des cartes d'offset. Cette fonction cycle sur l'ensemble des .atoms du dossier de séquence, charge le fichier atoms, construit deltaX, deltaY et offset puis envoie ces numpy array (taille = au nombre de coups dans la séquence) aux fonction update_atom_table dans laquelle on ajoute la position d'arrivée de chaque atome plus la valeur d'offset trouvée. Une fios les atomes chargés dans la base de donnée, on appelle la fonction update_pixel_table dont on charge la base de donnée, on la compare avec l'ensemble de positions d'arrivée et on la rajoute. 

        Parameters
        ----------
        sequence_directory : string
            le répertoire de la séquence que l'on veut ajouter
        """
        print(sequence_directory)
        # from the sequence directory, we get back the normalized sequence name e.g. 2022/07/15/003.
        boolean, normalized_sequence_name = self.get_normalized_sequence_name(
            sequence_directory
        )
        # If the directory does not look like a sequence directory, we do nothing. (boolean is set to False.)
        if not boolean:
            return
        # check
        if self.check_sequences_table(normalized_sequence_name):
            print("This sequence was already in the database !!!")
            return
        # On récupère tous les .atoms de la séquence

        selected_files = self.select_atoms_in_directory(sequence_directory)
        for file in tqdm(selected_files):
            atoms = np.fromfile(file, dtype="uint64")
            events_list = np.array(
                atoms.reshape(int(len(atoms) / 4), 4).T, dtype="int64"
            )
            deltaX = np.array(events_list[1] - events_list[0], dtype="int64")
            deltaY = np.array(events_list[3] - events_list[2], dtype="int64")
            offset = np.array(
                events_list[0] + events_list[1] - events_list[2] - events_list[3],
                dtype="int64",
            )
            self.update_atom_table(deltaX, deltaY, offset)
            self.update_pixel_table(deltaX, deltaY)
            del atoms
            del events_list
            del offset
            del deltaY
            del deltaX
            

        self.update_sequences_table(normalized_sequence_name)
        print(f"Sequence {normalized_sequence_name} added to the database :-).")

    def update_atom_table(self, deltaX, deltaY, offset):
        """Update the atom tables adding new offset values. We never load the atome table entirely since it is too large. When we add atoms to the sql database, ALWAYS USE if_exists = 'append' so that we never ever forget atoms.

        Parameters
        ----------
        deltaX : numpy array of integers
            positions sur X
        deltaY : numpy array of integers
            positions sur Y
        offset : numpy array of in
            offset à la position correspondante
        """
        self.connect_to_database()  # connexion to database
        new_atoms = pixels_to_compare = pd.DataFrame(
            np.transpose(np.array([deltaX, deltaY, offset], dtype="int64")),
            columns=["deltaX", "deltaY", "offset"],
        )
        new_atoms.to_sql("atoms", self.connexion, index=False, if_exists="append")
        self.connexion.commit()  # il faut commit pour prendre en compte l'ajout
        self.unconnect()

    def update_pixel_table(self, deltaX, deltaY):
        """Loads the pixel table from the database, compare it to the two lists deltaX, deltaY and adds pixel that was not already stored in the table.

        Parameters
        ----------
        deltaX : numpy array of integers
            positions sur X
        deltaY : numpy array of integers
            positions sur Y
        """
        # Load all pixel from the pixel database
        self.connect_to_database()  # connexion to database
        pixels = pd.read_sql_query("SELECT * FROM pixels", self.connexion)
        pixels_to_compare = pd.DataFrame(
            np.transpose(np.array([deltaX, deltaY], dtype="int64")),
            columns=["deltaX", "deltaY"],
        )
        # on vient de charger 2 dataframe pandas. on veut rajouter dans pixels les couples (X,Y) de pixels_to_compare qui ne sont pas dans pixels.
        pixels = pixels.merge(
            pixels_to_compare, on=["deltaX", "deltaY"], how="outer"
        ).drop_duplicates()
        pixels.to_sql("pixels", self.connexion, index=False, if_exists="replace")
        self.connexion.commit()  # il faut commit pour prendre en compte l'ajout
        self.unconnect()

    def update_sequences_table(self, normalized_sequence_name):
        """ajoute normalized_sequence_name to the sequences table

        Parameters
        ----------
        normalized_sequence_name : string like object
            nom normalisé de la séquence (eg '2022/07/15/020')
        """
        self.connect_to_database()  # connexion to database
        self.cursor.execute(
            f"INSERT INTO sequences VALUES ('{normalized_sequence_name}')"
        )
        self.connexion.commit()  # il faut commit pour prendre en compte l'ajout
        self.unconnect()

    def check_sequences_table(self, normalized_sequence_name):
        """Vérifie si normalized_sequence_name est dans la table sequences de la base de données.

        Parameters
        ----------
        normalized_sequence_name : string
            le nom normalisé de la séquence

        Returns
        -------
        True if the sequence was already add to database, False otherwise
        """
        self.cursor.execute(
            f"SELECT EXISTS (SELECT * FROM sequences WHERE name = '{normalized_sequence_name}')"
        )

        result = self.cursor.fetchall()[0][
            0
        ]  # if there is something already, le curseur va renvoyer [(1,)] et s'il n'y a rien : [(0,)].
        # On a donc result = 1 --> il y a déjà une séquence : on renvoie True
        # result = 0 --> ill n'ya apas encore de séquence : on renvoie false.
        self.unconnect()
        if result == 0:
            return False
        else:
            return True

    def select_atoms_in_directory(self, sequence_directory):
        """Select all the .atoms files in the sequence_directory. If the directory is empty, returns empty list.

        Parameters
        ----------
        sequence_directory : string or path like object
            the directory in which we want to obtains all .atoms files

        Returns
        -------
        list of string/path like objects
            liste des fichiers .atoms du sequence_directory.
        """
        ext = [".atoms"]
        sequence_directory = Path(sequence_directory)
        selected_files = sorted(
            [
                path.as_posix()
                for path in filter(
                    lambda path: path.suffix in ext, sequence_directory.glob("*")
                )
            ]
        )
        return selected_files

    def get_normalized_sequence_name(self, sequence_directory):
        """Retourne le nom normalisé de la séquence ajoutée à la base de données. Par exemple, si on donne à cette méthode le string "/home/victor/gus_data/2022/08/01/004" cette fonction retourne "2022/08/01/004" qui est le nom normalisé de la séquence (ne dépend pas du disque sur lequel on la regarde)

        Parameters
        ----------
        sequence_directory : strin or pathlike object
            directory of the sequence we look for

        Returns
        -------
        _type_
            _description_
        """
        pattern = "20[0-9][0-9]/[0-9][0-9]/[0-9][0-9]/[0-9][0-9][0-9]"
        directory = re.findall(pattern, str(sequence_directory))
        if not directory:
            print(
                "Sorry but I cannot extract the date and the reference of the sequence hence I will not add it to my databse."
            )
            return False, ""
        normalized_sequence_name = directory[0]
        return True, normalized_sequence_name

    ##################################################################
    ########## METHDS TO UPDATE THE MCP MAPS
    ##################################################################

    def update_result(self):
        """This methods update the result dataframe that contains statisctics properties of our mcp map.
        The procedure to update it is the following :

        1. we load all the possible pixels of the MCP
        2. We loop over the possible pixels and load from the database the offsets associated to the pixel location.
        3. For each pixel we compute statistical properties aka mean, counts, std, max and min.
        4. we save it into a pkl file.
        """
        ###
        self.load_pixels()
        # INITIALISATION
        self.result = copy.deepcopy(self.pixels)
        self.result["mean"] = 0
        self.result["std"] = 0
        self.result["count"] = 0
        self.result["max"] = 0
        self.result["min"] = 0
        self.connect_to_database()
        serideltaX_serie = self.result.deltaX.unique()
        self.result = self.result.set_index(["deltaX", "deltaY"])

        for value in tqdm(serideltaX_serie):
            a = pd.read_sql_query(
                f"SELECT deltaY, offset FROM atoms WHERE deltaX = {value};",
                self.connexion,
            )
            a = a.groupby("deltaY")
            df = a.mean()
            df.columns = df.columns.str.replace("offset", "mean")

            df["std"] = a.std()
            df["count"] = a.count()
            df["max"] = a.max()
            df["min"] = a.min()
            df.reset_index(inplace=True)
            df = df.fillna(0)
            df["deltaX"] = value
            df = df.set_index(["deltaX", "deltaY"])
            self.result.update(df)
        self.result.reset_index(inplace=True)
        """
        # On parcourt chaque pixel du MCP
        for i, row in tqdm(self.result.iterrows()):
            self.cursor.execute(
                f"SELECT offset FROM atoms WHERE deltaX = {row.deltaX} AND deltaY = {row.deltaY};"
            )
            offset_datas = self.cursor.fetchall()
            offset_datas = np.array(offset_datas)
            self.result.at[i, "mean"] = offset_datas.mean()
            self.result.at[i, "std"] = offset_datas.std()
            self.result.at[i, "count"] = offset_datas.size
            self.result.at[i, "max"] = offset_datas.max()
            self.result.at[i, "min"] = offset_datas.min()
        """
        self.unconnect()

        self.save_result()

    ##################################################################
    ########## METHDS TO SHOW THE MCP MAPS
    ##################################################################
    def show_detectivity(self):
        mean = np.array(
            self.result.pivot(index="deltaX", columns="deltaY", values="mean").fillna(0)
        )
        import matplotlib.pyplot as plt

        plt.imshow(mean)
        plt.show()


if __name__ == "__main__":
    mcp_map = MCP_offset_map("/home/victor/mcpmaps/")
    mcp_map.add_sequence_to_database("/mnt/manip_E/2022/07/15/007")
    # mcp_map.update_result()
    
