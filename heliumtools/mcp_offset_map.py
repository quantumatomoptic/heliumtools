#!/usr/bin/env python
# -*- mode:Python; coding: utf-8 -*-

# ----------------------------------
# Created on the 22-7-2022 by Victor
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
La table atoms contient la position et l'offset de chaque atome détecté. Elle est donc assez grosse. 

#### USER REMARKS ######

This class is usefull when one requires to build MCP offset maps. It works as follows :
--> the only argument you give to this class is the path (unix string like object) to the folder in which you want to store the data or in which the data already exist. 
--> If it is empty, it will creat a empty SQLite database. First thing you have to do is to add sequences to this database. To do so, use the instruction 'add_sequence_to_database('Path/To/your/sequences'). In this case, the programm tries to get the reference of this sequence (of the form yyy/mm/dd/seq) and search into its first database (whos name is sequences) if this sequences was already add to the database or not. If not, it adds it (the 'atoms' table), looping over all cycles of the sequence. Such an operation takes a few minutes for really large sequences.
--> What we add to the atom database are the following elements : 
    + deltaX : position on X of the atom,
    + deltaY : position on Y of the atom,
    + offset : offset of the atom, X2+X1-Y2-Y2
--> Once you've built the database, you can compute the statistical properties of your MCP. To do so, you just have to call update_result() : this function will loop over all X pixels of the MCP and load atoms that where found at this X position : we calculate the number of count, the mean offset, std, min and max at each Y position and store it into a dataframe named 'result'. This dataframe (size : 1500x1500x8) is then store into a pkl file whose name is result.pkl. Note that this step is quite long (an hour typically). 
Note also that if your compute have huge memmory, you might be able to load all the atom database in one shot using :
    mcp_map = MCP_offset_map("your/folder/database")
    mcp_map.connect_to_database()
    all_atoms = pd.read_sql_query("SELECT * FROM atoms", mcp_map.connexion)


#### DEVELOPPERS REMARKS ######

Follows some usefull commands to deal with sql database

*** Basis to use SQL ***

self.database_name = /.../offset.db   # le nom de notre base de donnée
self.connexion = sqlite3.connect(self.database_name)   # cela permet de nous connecter 
self.cursor = self.connexion.cursor() # le curseur permettant de nous y connecter.
NB : les 2 étapes ci-dessus sont contenues dans 'linstruction self.connect_to_database() qui est appelée chaque fois qu'on veut s'y connecter (d'où le nom...).

Pour envoyer une requète à la base de données, on utilise l'instruction execute du curseur : 
self.cursor.execute("CREATE TABLE sequences (name NVARCHAR(50))") # on crée une table avec une seule colonne du type string de max 50 caractères.

Pour récupérer le résultat d'une requète SQL, il faut envoyer la requète avec notre commande puis "fetch" pour récupérer le résultat. 
self.cursor.execute(f"SELECT EXISTS (SELECT * FROM sequences WHERE name = '2022/07/157007')")
resultat = self.cursor.fetchall()
# Dans ce cas, la commande résultat sera une liste de tuples. Puisque la commande ci-dessous appelle une réponse booléenne, la forme de résultat sera [(0,),] ou [(1,),] 


*** SQL and Pandas ***
pixels = pd.read_sql_query(("SELECT * FROM pixels", self.connexion)) # charge les données de la table 'pixels' de la base de donnée dans un pandas dataframe
dataframe.to_sql(name_of_table, self.connexion, if_exists='fail',index = True) 
# Always use index = False
# For atoms : use if_exists = 'append'   ;because we always want to append the atom (if we already add an atom with this offset, we have to add it and not to drop it).
"""

import os, re, json, copy
from pathlib import Path
from select import select
import pandas as pd
import numpy as np
import sqlite3 as db
from tqdm import tqdm
from pylab import cm
import matplotlib.pyplot as plt


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
        Dans l'idée de ce code, on stock tous les données dans un répertoire dont l'adresse est self.path_to_folder. Celui-ci doit contenir :
        une base de donnée "offset.db",
        un fichier result au format pkl.

        Cette méthode s'assure que le dossier existe et si non, elle le crée. Elle se connecte ensuite à la base de donnée et récupère l'ensemble des nom des tables qui y sont stockée (listofTables) puis


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
        # On récupère une liste de la forme suivante : [("sequences", ), ("atoms", ),] si les tables existent et une liste vide sinon.
        if ("sequences",) not in listofTables:
            self.cursor.execute("CREATE TABLE sequences (name NVARCHAR(50))")
        if ("atoms",) not in listofTables:
            self.cursor.execute(
                "CREATE TABLE atoms (deltaX INT, deltaY INT, offset INT)"
            )
        self.connexion.commit()  # il faut commit pour prendre en compte l'ajout
        self.unconnect()

        # On charge maintenant le fichier results
        self.load_results()

    def load_results(self):
        if os.path.isfile(self.result_name):
            self.result = pd.read_pickle(self.result_name)
        else:
            self.result = pd.DataFrame(
                columns=["deltaX", "deltaY", "mean", "std", "counts", "max", "min"]
            )

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
        self.connect_to_database()
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
        new_atoms = pd.DataFrame(
            np.transpose(np.array([deltaX, deltaY, offset], dtype="int64")),
            columns=["deltaX", "deltaY", "offset"],
        )
        new_atoms.to_sql("atoms", self.connexion, index=False, if_exists="append")
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
        self.connect_to_database()
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
        print(" --------------------------------------------------- ")
        print(" /!\ WARNING /!\ ")
        print(" THIS PROGRAMM IS REALLY SLOW")
        print("Are you sure the result is not already available ?!?")
        ###
        self.load_pixels()
        # INITIALISATION

        self.result = pd.DataFrame(
            columns=["deltaX", "deltaY", "mean", "std", "counts", "max", "min"]
        )

        self.load_pixels()
        self.connect_to_database()
        print(self.pixelsX)
        for value in tqdm(self.pixelsX.deltaX):
            print(f"SELECT deltaY, offset FROM atoms WHERE deltaX = {value};")
            a = pd.read_sql_query(
                f"SELECT deltaY, offset FROM atoms WHERE deltaX = {value};",
                self.connexion,
            )
            a = a.groupby("deltaY")
            df = a.mean()
            df.columns = df.columns.str.replace("offset", "mean")

            df["std"] = a.std()
            df["counts"] = a.count()
            df["max"] = a.max()
            df["min"] = a.min()
            df.reset_index(inplace=True)  # --> deltaY
            df = df.fillna(0)
            df["deltaX"] = value

            self.result = pd.concat([self.result, df])
            del df

        self.unconnect()

        self.save_result()

    ##################################################################
    ########## METHDS TO SHOW THE MCP MAPS
    ##################################################################
    def show_detectivity(self):
        """montre un image 2D du nombre de coups par pixel"""
        counts = np.array(
            self.result.pivot(index="deltaX", columns="deltaY", values="counts").fillna(
                0
            )
        )
        fig, ax = plt.subplots()
        shw = ax.imshow(counts)
        bar = plt.colorbar(shw)
        plt.title("Gain map")
        plt.show()

    def show_dead_pixels(self):
        """montre les pixels qui n'ont reçu aucun atome"""
        maxi_counts = self.result.counts.max()
        counts = np.array(
            self.result.pivot(index="deltaX", columns="deltaY", values="counts").fillna(
                3 * maxi_counts
            )
        )
        counts = np.floor(counts / (2 * maxi_counts))
        fig, ax = plt.subplots()
        shw = ax.imshow(counts)
        plt.title("Pixels with zero-detectivity")
        plt.show()

    def show_offset(self):
        """représente la carte des offsets"""
        offset_mean = self.result["mean"].mean()
        print(offset_mean)
        offsets = np.array(
            self.result.pivot(index="deltaX", columns="deltaY", values="mean").fillna(0)
        )
        fig, ax = plt.subplots()
        shw = ax.imshow(offsets)
        bar = plt.colorbar(shw)
        plt.title("Offset map (mean is {:.0f})".format(offset_mean))
        plt.show()

    def show_stds(self):
        """représente la carte des déviations standards des offsets"""
        stds = np.array(
            self.result.pivot(index="deltaX", columns="deltaY", values="std").fillna(0)
        )
        fig, ax = plt.subplots()
        shw = ax.imshow(stds)
        bar = plt.colorbar(shw)
        plt.title("Offset stds map ")
        plt.show()

    def show_offset_cuts_alongY(self, coupe1=-400, coupe2=0, coupe3=200):
        """Jolie figure avec 3 coupes selon Y"""
        counts = np.array(
            self.result.pivot(index="deltaY", columns="deltaX", values="counts").fillna(
                0
            )
        )
        offset = np.array(
            self.result.pivot(index="deltaY", columns="deltaX", values="mean").fillna(0)
        )
        std = np.array(
            self.result.pivot(index="deltaY", columns="deltaX", values="std").fillna(0)
        )
        offsetX = self.result.deltaX.min()
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
        # Première figure : détectivité
        shw = ax1.imshow(counts)
        # ax.set_title("Gain map")
        from matplotlib.patches import Rectangle

        ax1.add_patch(Rectangle((coupe2 - offsetX - 3, 0), 6, 1400))
        ax1.add_patch(Rectangle((coupe3 - offsetX - 3, 0), 6, 1400))

        # Deuxièmes figures
        coupes = [coupe1, coupe2, coupe3]
        coupes_axes = [ax2, ax3, ax4]
        my_palette = [
            "indianred",
            "steelblue",
            "olivedrab",
            "goldenrod",
            "darkslategrey",
        ]
        print(self.result.deltaX.drop_duplicates())
        for i in range(len(coupes)):
            print(i)
            color = my_palette[i]
            cut = coupes[i]
            ax1.add_patch(Rectangle((cut - offsetX - 3, 0), 6, 1400, color=color))

            datas = self.result[self.result["deltaX"] == cut]
            # print(datas)
            ax = coupes_axes[i]
            x = datas["deltaY"]
            y = datas["mean"]
            u_y = datas["std"].fillna(0)

            # ax.errorbar(x, y, yerr=u_y, fmt=".")
            ax.plot(x, y, color, label=f"x={cut}")
            ax.plot(x, y - u_y, alpha=0.5, color=color)
            ax.plot(x, y + u_y, alpha=0.5, color=color)
            ax.set_ylim(0, 80)
            ax.set_xlabel("Position on x")
            ax.set_ylabel("Offset (with resolution)")
            ax.legend()
        plt.show()

    def show_offset_distribution(
        self, points=[(300, 300), (-150, 150), (-200, -200), (100, -100)]
    ):
        """affiche la distribution des offset pour l'ensemble des points (deltaX, deltaY) de la liste points

        Parameters
        ----------
        points : list of tuples of integers (deltaX, deltaY)
            points dont on veut afficher la distribution, by default [(300, 300), (600,500), (-200, -200), (100, 100)]
        """
        from math import sqrt, ceil

        n_cols = int(sqrt(len(points)))
        n_lines = int(ceil(len(points) / n_cols))
        fig, axs = plt.subplots(
            n_lines, n_cols, figsize=(10, 6), constrained_layout=True
        )
        self.connect_to_database()
        for ax, p in zip(axs.flat, points):
            deltaX = p[0]
            deltaY = p[1]
            offset = np.random.normal(size=100)
            offset = pd.read_sql_query(
                f"SELECT offset FROM atoms WHERE deltaX = {deltaX} AND deltaY = {deltaY};",
                self.connexion,
            )
            ax.hist(offset, bins=20, label=str(p))
            ax.legend()
        self.unconnect()
        plt.show()


if __name__ == "__main__":
    mcp_map = MCP_offset_map("/home/victor/mcpmaps/")
    # mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/00/010") # --> FAKE SEQUENCE !!!!!
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/15/007")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/15/029")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/18/004")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/18/005")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/19/004")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/19/005")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/19/007")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/19/008")
    mcp_map.add_sequence_to_database("/home/victor/gus_data/2022/07/19/014")
    # mcp_map.update_result()
    # mcp_map.show_detectivity()
    # mcp_map.show_dead_pixels()
    # mcp_map.show_offset()
    # mcp_map.show_stds()
    # mcp_map.show_offset_cuts_alongY()
    # mcp_map.show_offset_distribution()
