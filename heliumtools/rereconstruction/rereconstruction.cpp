// ----------------------------------
//  Created on the Wed Nov 02 2022 by Victor
//
//  Developped by Victor, ...
//
//  Last (big) change on the ... by ...
//
//  Copyright (c) 2022 - Helium1@LCF
// ----------------------------------
//

/*
Ce code permet de re-reconstruire les atomes d'un répertoire. L'utilisation est la suivante :
on écrit dans un fichier de configuration
*/
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// Autre
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/types.h>
#include <dirent.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>
using namespace std;
typedef uint64_t timedata;
struct atomdata
{
    timedata TX1, TX2, TY1, TY2;
};

struct paramstruct // reconstruction parameters
{
    double res = 120E-12;                  // TDC resolution (seconds) -> 0.12e-9
    double evgate = 85.00E-9;              // time gate to identify matching detection events
    double atgate = 5.00E-9;               // time gate to identify matching events that should be counted as an atom
    double deltaTgate = 10.00E-9;          // time gate maximale une fois qu'on a obtenu un quadruplet
    std::string seq_folder = "salut";      // le fichier depuis lequel on veut charger la séquence
    std::string new_seq_folder = "salut2"; // le nouveau dossier dans lequel on veut enregistrer les données
    int reconstruction_number = 1;
    bool use_offset_map = true;
    bool keep_all_potential_atoms = true;
};

typedef std::vector<std::string> stringvec;

//
// Definition des fonctions utilisées dans le code
bool read_configuration_file(string filepath, paramstruct &parameters_p);
void get_cycles(const std::string directory, stringvec &cycles_p);
void load_offset_map(string offset_map_path, std::vector<std::vector<int>> &offset_p);
list<timedata> load_timedata(std::string path_to_time);
void reconstruction1(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params);
/*
==============================================================================================================
FONCTION PRINCIPALE DU CODE
==============================================================================================================
*/
int main()
{
    string configuration_file = "conf.txt";
    string offset_map_path = "Offset.txt";
    paramstruct params;
    bool failure = false;

    // ####
    // Lecture du fichier de configuration conf.txt
    // ####

    failure = read_configuration_file(configuration_file, params); // la fonction agit sur param
    if (failure)
        exit(0);
    // ####
    // Chargement de la carte d'offset
    // ####
    std::vector<std::vector<int>> offset;
    load_offset_map(offset_map_path, offset);

    cout << params.seq_folder << "\n";

    // ####
    // Creation du fichier de séquence params.new_seq_folder
    // ####
    // const char *dirnameChars = params.new_seq_folder.c_str(); // conversion en char pour créer le dossier de séquence
    // if (mkdir(dirnameChars, 0777) == -1)
    // {
    //     cerr
    //         << "Error while creating the new folder :  " << strerror(errno) << endl;
    //     exit(0);
    // }
    // else
    //     cout << "\n Directory created";

    // ####
    // On récupère la liste de tous les cycles qu'il faut rereconstruire.
    // ####
    stringvec cycles;
    get_cycles(params.seq_folder, cycles);

    // ####
    // On parcourt tous les cycles de la séquence et on rereconstruit les atomes.
    // ####
    // Définition des variables dont nous aurons l'usage plus tard.
    list<timedata> X1, X2, Y1, Y2;
    vector<atomdata> atoms;
    map<string, float> reconstruction_statistics; // reconstruction_statistics is a dictionnary like

    cout << "\n ... Starting to load .times files... \n";
    for (auto i = cycles.begin(); i != cycles.end(); ++i)
    { // ------ LOAD TIMES ------
        X1 = load_timedata(params.seq_folder + "/" + *i + ".timesx1");
        X2 = load_timedata(params.seq_folder + "/" + *i + ".timesx2");
        Y1 = load_timedata(params.seq_folder + "/" + *i + ".timesy1");
        Y2 = load_timedata(params.seq_folder + "/" + *i + ".timesy2");
        reconstruction_statistics["Number of X1"] = X1.size();
        reconstruction_statistics["Number of X2"] = X2.size();
        reconstruction_statistics["Number of Y1"] = Y1.size();
        reconstruction_statistics["Number of Y2"] = Y2.size();
        cout << "Starting reconstruction for cycle " << *i << "...";
        if (params.reconstruction_number == 1)
        {
            reconstruction1(X1, X2, Y1, Y2, atoms, params);
            cout << atoms.front().TX1 << "\n";
        }
    }
}

/*
==============================================================================================================
FONCTIONS DE RECONSTRUCTION APPELÉES DANS LE CODE
==============================================================================================================
*/
void reconstruction1(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params)
{
    // gate X, gatY->bulb width
    timedata gateY = timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);
    atoms_p.push_back(atomdata{0, 0, 0, 0});
}
/*
==============================================================================================================
AUTRES FONCTIONS APPELÉES DANS LE CODE
==============================================================================================================
*/

list<timedata> load_timedata(std::string path_to_time)
{ // Cette fonction charge des données binaires du fichier dont
    // l'adresse est path_to_time et retourne le contenu dans une liste.
    //  https://stackoverflow.com/questions/59718584/reading-a-binary-file-and-interpret-as-integers
    uint64_t event;
    list<timedata> buffer;
    ifstream file(path_to_time, ios::in | ios::binary);
    if (file.is_open())
    {
        while (file.read((char *)&event, sizeof(event)))
        {
            buffer.push_back(event);
        }
    }
    return buffer;
}

bool String2Convert(std::string var)
{ // convertie un string vers un booléen.
    if ((var == "true") || (var == "TRUE") || (var == "True") || (var == "1"))
    {
        return true;
    }
    else
    {
        return false;
    }
}
bool read_configuration_file(string filepath, paramstruct &parameters_p)
{ // Cette fonction lit le fichier de configuration filepath (à priori conf.txt). Elle lit les lignes une à une et attribut les valeurs dans l'ordre aux éléments de la structure de type parameters_structure. L'ordre des paramètres est donc très important.
    // args :
    //      filepath : chemin vers le fichier de configuration. À priori 'conf.txt'.
    //      parameters_p : pointeur vers une structure du type paramètre. À chaque ligne lue, on remplit la valeur du paramètre correspondant.

    std::ifstream myfile;
    myfile.open(filepath);
    // on lit le fichier ligne à ligne s'il existe
    if (!(myfile))
    {
        std::cout << "I  CANNOT OPEN THE CONFIG FILE\n";
        return true;
    }
    std::cout << "Startign to read the config file " + filepath + "\n";
    std::cout
        << "\n \n === Parameters for the sequence === \n \n";
    std::string myline;
    // LIGNE 1 : le dossier initial de la séquence.
    std::getline(myfile, myline);
    std::cout << "Treating atoms from folder ";
    std::cout << myline << '\n';
    parameters_p.seq_folder = myline;

    // LIGNE 2 : le dossier où on va sauver les .atoms
    std::getline(myfile, myline);
    if (!(myfile))
    {
        std::cout << "Error during reading the conf file : no line 2 -->  where do we save atoms ???";
        return true;
    }
    std::cout << "Saving atoms in " << myline << '\n';
    parameters_p.new_seq_folder = myline;

    // LIGNE 3 : le numéro du programme de reconstruction qu'on veut utiliser
    std::getline(myfile, myline);
    if (!(myfile))
    {
        std::cout << "Error during reading the conf file : no line 3 -->  which program do we want to use ? 1 is default";
        return true;
    }
    std::cout << "The programm number chosen is " << myline << "\n";
    parameters_p.reconstruction_number = stoi(myline);

    // // LIGNE 3 : si on veut utiliser la carte d'offset
    // std::getline(myfile, myline);
    // if (!(myfile))
    // {
    //     std::cout << "Error during reading the conf file : no line 3";
    //     return true;
    // }
    // parameters_p.use_offset_map = String2Convert(myline);
    // std::cout << "Use offset map boolen --> " << parameters_p.use_offset_map << "\n";

    // // LIGNE 4 : si on veut garder tous les atomes potentiels
    // std::getline(myfile, myline);
    // if (!(myfile))
    // {
    //     std::cout << "Error during reading the conf file : no line 4 \n";
    //     return true;
    // }
    // parameters_p.keep_all_potential_atoms = String2Convert(myline);
    // std::cout << "Keep ALL potential atoms boolean --> " << parameters_p.keep_all_potential_atoms << "\n";

    std::cout << "\n  \n";
    return false;
}

void get_cycles(const std::string directory, stringvec &cycles_p)
{ // retourne l'ensemble des fichiers .atoms du dossier directory et peuple le vecteur de string cycles (pointeur en argument).
    std::string extension = ".atoms";
    DIR *dirp = opendir(directory.c_str());
    struct dirent *
        dp;
    while ((dp = readdir(dirp)) != NULL)
    {
        std::string filename = dp->d_name;
        std::size_t ind = filename.find(extension); // Find the starting position of extension in the string
        if (ind != std::string::npos)
        {
            filename.erase(ind, extension.length()); // erase function takes two parameter, the starting index in the string from where you want to erase characters and total no of characters you want to erase.
            cycles_p.push_back(filename);
        }
    }
    closedir(dirp);
}

void load_offset_map(string offset_map_path, std::vector<std::vector<int>> &offset_p)
{ // Load the Offset Map : location is in the git repository of the TDC driver.
    std::ifstream file(offset_map_path);
    const std::string &delim = " \t";
    string line;
    string strnum;
    // clear first
    offset_p.clear();
    cout << "Reading offset map " + offset_map_path + "..." << endl;
    // parse line by line
    while (getline(file, line))
    {
        offset_p.push_back(vector<int>());

        for (string::const_iterator i = line.begin(); i != line.end(); ++i)
        {
            // If i is not a delim, then append it to strnum
            if (delim.find(*i) == string::npos)
            {
                strnum += *i;
                if (i + 1 != line.end()) // If it's the last char, do not continue
                    continue;
            }

            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if (strnum.empty())
                continue;

            // If we reach here, we got a number. Convert it to double.
            double number;

            istringstream(strnum) >> number;
            offset_p.back().push_back(number);

            strnum.clear();
        }
    }
    std::cout << "Done. \n \n";
}