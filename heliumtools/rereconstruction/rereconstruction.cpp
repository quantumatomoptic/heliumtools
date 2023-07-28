// ----------------------------------
//  Created on the Wed Nov 02 2022 by Victor
//
//  Developped by Victor, ...
//
//  Last (big) change on the ... by ...
//  12/07/23 (Victor) It happened on the 03/07, TDC time lists were not sorted and this caused the algorithm (3 & 4) not to work properly. To avoid this, we now sort the list when loading them.
//
//  Copyright (c) 2022 - Helium1@LCF
// ----------------------------------
//

/*
This code ables to reconstruct atoms a posteriori. It should not be mandatory to recompile it to make it work on a linux computer. One should just change the configuration file conf.txt and enter
Firts line : sequence directory in which there are .times
Second line : directory (not yet created) in which new .atoms will be stored
Third line : the programm number the user wants to use
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
using namespace std::chrono; // to measure the execution time of the programm

typedef uint64_t timedata;
struct atomdata
{
    timedata TX1, TX2, TY1, TY2;
};

struct paramstruct // reconstruction parameters
{
    double res = 120E-12;                  // TDC resolution (seconds) -> 0.12e-9 ; 275e-12 old MCP
    double evgate = 85.00E-9;              // time gate to identify matching detection events --> retard line.
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
bool reconstruction1(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params);
bool reconstruction2(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params,
                     std::vector<std::vector<int>> &offset_p);
bool reconstruction3(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params,
                     std::vector<std::vector<int>> &offset_p);
bool reconstruction4(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params);

bool reconstruction5(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params);

void WriteStatistics(string filename, map<string, float> &reconstruction_statistics);
void WriteAtoms(string filename, vector<atomdata> &atoms);
void copy_paste_file(string input_file, string output_file);
list<int> get_offset_of_atoms(vector<atomdata> &atoms_p, std::vector<std::vector<int>> &offset_p);
void WriteOffsets(string filename, list<int> &offset_values_p);
timedata AbsDiff(timedata T1, timedata T2);
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
    const char *dirnameChars = params.new_seq_folder.c_str(); // conversion en char pour créer le dossier de séquence
    if (mkdir(dirnameChars, 0777) == -1)
    {
        cerr
            << "Error while creating the new folder :  " << strerror(errno) << endl;
        exit(0);
    }
    else
        cout << "\n Directory created";
    string inputfile = "rereconstruction.cpp";
    string outputfile = params.seq_folder + "/rereconstruction.cpp";
    copy_paste_file("rereconstruction.cpp", params.new_seq_folder + "/rereconstruction.cpp");
    copy_paste_file("conf.txt", params.new_seq_folder + "/conf.txt");
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
        string inputfile = params.seq_folder + "/" + *i;
        string outputfile = params.new_seq_folder + "/" + *i;
        X1 = load_timedata(inputfile + ".timesx1");
        X1.sort();
        X2 = load_timedata(inputfile + ".timesx2");
        X2.sort();
        Y1 = load_timedata(inputfile + ".timesy1");
        Y1.sort();
        Y2 = load_timedata(inputfile + ".timesy2");
        Y2.sort();
        reconstruction_statistics["Number of X1"] = X1.size();
        reconstruction_statistics["Number of X2"] = X2.size();
        reconstruction_statistics["Number of Y1"] = Y1.size();
        reconstruction_statistics["Number of Y2"] = Y2.size();
        // ------- START RECONSTRUCTION --------
        // Choose reconstruction number between :
        // 1 : reconstruction until september 2022 (no offset map, first atom that match is chosen)
        // 2 : reconstruction used from september 2022 with offset map but same reconstruction code.
        cout << "\n         --------------------------------------------" << endl;
        cout << "         Starting reconstruction for cycle " << *i << "... " << endl;
        cout << "         --------------------------------------------" << endl;
        cout << "Input file : " + inputfile << endl;
        cout << "Output file : " + outputfile << endl;
        copy_paste_file(inputfile + ".json", outputfile + ".json");
        copy_paste_file(inputfile + ".sequence_parameters", outputfile + ".sequence_parameters");
        auto t_start = high_resolution_clock::now();
        if (params.reconstruction_number == 1)
        {
            reconstruction1(X1, X2, Y1, Y2, atoms, params);
        }
        if (params.reconstruction_number == 2)
        {
            reconstruction2(X1, X2, Y1, Y2, atoms, params, offset);
        }
        if (params.reconstruction_number == 3)
        {
            reconstruction3(X1, X2, Y1, Y2, atoms, params, offset);
        }
        if (params.reconstruction_number == 4)
        {
            reconstruction4(X1, X2, Y1, Y2, atoms, params);
        }
        if (params.reconstruction_number == 5)
        {
            reconstruction5(X1, X2, Y1, Y2, atoms, params);
        }
        auto t_end_reconstruction = high_resolution_clock::now();
        auto duration_reconstruction = duration_cast<milliseconds>(t_end_reconstruction - t_start);
        reconstruction_statistics["Reconstruction program number"] = float(params.reconstruction_number);
        reconstruction_statistics["Reconstruction duration"] = float(duration_reconstruction.count()) / 1000; // en secondes
        float atoms_size = atoms.size();
        reconstruction_statistics["Number of atoms"] = atoms.size();
        reconstruction_statistics["Reconstruction rate of X1"] = round(100 * atoms_size / reconstruction_statistics["Number of X1"]);
        reconstruction_statistics["Reconstruction rate of X2"] = round(100 * atoms_size / reconstruction_statistics["Number of X2"]);
        reconstruction_statistics["Reconstruction rate of Y1"] = round(100 * atoms_size / reconstruction_statistics["Number of Y1"]);
        reconstruction_statistics["Reconstruction rate of Y2"] = round(100 * atoms_size / reconstruction_statistics["Number of Y2"]);
        WriteAtoms(outputfile + ".atoms", atoms);
        WriteStatistics(outputfile + ".stats", reconstruction_statistics);
        // list<int> offset_of_atoms;
        // offset_of_atoms = get_offset_of_atoms(atoms, offset);
        // WriteOffsets(outputfile + ".offsets", offset_of_atoms);

        // Clean up
        X1.clear();
        X2.clear();
        Y1.clear();
        Y2.clear();
        atoms.clear();
    }
}

/*
==================================================
FONCTIONS DE RECONSTRUCTION APPELÉES DANS LE CODE
==================================================
*/
timedata AbsDiff(timedata T1, timedata T2)
{
    if (T1 > T2)
        return T1 - T2;
    else
        return T2 - T1;
}

bool reconstruction1(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params)
{
    cout << "Reconstruction 1 : Ziyad style" << endl;
    // gate X, gatY->bulb width
    timedata gateY = timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);
    while (X1_p.begin() != X1_p.end())
    {
        // We first definitively get rid of all events occurring before the bulb start on X2, Y1, Y2
        // By construction, they will never match with a later event on X1
        if (*X1_p.begin() > gateX)
            while (X2_p.begin() != X2_p.end() && *X2_p.begin() < (*X1_p.begin() - gateX))
                X2_p.erase(X2_p.begin());
        if (*X1_p.begin() > gateY)
        {
            while (Y1_p.begin() != Y1_p.end() && *Y1_p.begin() < (*X1_p.begin() - gateY))
                Y1_p.erase(Y1_p.begin());
            while (Y2_p.begin() != Y2_p.end() && *Y2_p.begin() < min((*X1_p.begin() - gateY), *X1_p.begin()))
                Y2_p.erase(Y2_p.begin());
        }
        // We then look for events that can correspond to an atom
        // As soon as we have found such events, we erase them from the list and keep going along X1
        bool atomfound = false;
        auto searchX2 = X2_p.begin(), searchY1 = Y1_p.begin(), searchY2 = Y2_p.begin();
        while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX) && !atomfound)
        {
            while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY) && !atomfound)
            {
                while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY) && !atomfound)
                {
                    timedata TX1 = *X1_p.begin();
                    timedata TX2 = *searchX2;
                    timedata TY1 = *searchY1;
                    timedata TY2 = *searchY2;

                    timedata dTX = AbsDiff(TX1, TX2);
                    timedata dTY = AbsDiff(TY1, TY2);
                    timedata TX = (TX1 + TX2) / 2;
                    timedata TY = (TY1 + TY2) / 2;
                    // distance to the MCP centre
                    timedata dist = timedata(sqrt(pow(dTX, 2) + pow(dTY, 2)));

                    // time difference between events on X and Y
                    // its a curious way of calculating an absolute value
                    // but we want to make sure that we don't loose precision (time coded on 64bits)
                    timedata dT = AbsDiff(TX, TY);

                    // an atom would be inside the MCP radius and fall atoms the same time on X and Y
                    if (dist < MCPdiameter && dT < deltaT)
                    {
                        atoms_p.push_back(atomdata{TX1, TX2, TY1, TY2});
                        atomfound = true;
                    }
                    ++searchY2;
                }
                ++searchY1;
            }
            ++searchX2;
        }
        if (atomfound)
        {
            X2_p.erase(prev(searchX2));
            Y1_p.erase(prev(searchY1));
            Y2_p.erase(prev(searchY2));
        }

        X1_p.erase(X1_p.begin());
    }
    if (atoms_p.empty())
    {
        return false;
    }
    else
        return true;
}

bool reconstruction2(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params,
                     std::vector<std::vector<int>> &offset_p)
{
    cout << "Reconstruction 2 : Jean-Louis Soudure style" << endl;
    // gate X, gatY->bulb width
    timedata gateY = timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);
    while (X1_p.begin() != X1_p.end())
    {
        // We first definitively get rid of all events occurring before the bulb start on X2, Y1, Y2
        // By construction, they will never match with a later event on X1
        if (*X1_p.begin() > gateX)
            while (X2_p.begin() != X2_p.end() && *X2_p.begin() < (*X1_p.begin() - gateX))
                X2_p.erase(X2_p.begin());
        if (*X1_p.begin() > gateY)
        {
            while (Y1_p.begin() != Y1_p.end() && *Y1_p.begin() < (*X1_p.begin() - gateY))
                Y1_p.erase(Y1_p.begin());
            while (Y2_p.begin() != Y2_p.end() && *Y2_p.begin() < min((*X1_p.begin() - gateY), *X1_p.begin()))
                Y2_p.erase(Y2_p.begin());
        }
        // We then look for events that can correspond to an atom
        // As soon as we have found such events, we erase them from the list and keep going along X1
        bool atomfound = false;
        auto searchX2 = X2_p.begin(), searchY1 = Y1_p.begin(), searchY2 = Y2_p.begin();
        while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX) && !atomfound)
        {
            while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY) && !atomfound)
            {
                while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY) && !atomfound)
                {
                    timedata TX1 = *X1_p.begin();
                    timedata TX2 = *searchX2;
                    timedata TY1 = *searchY1;
                    timedata TY2 = *searchY2;

                    timedata dTX = AbsDiff(TX1, TX2);
                    timedata dTY = AbsDiff(TY1, TY2);
                    timedata TX = (TX1 + TX2) / 2;
                    timedata TY = (TY1 + TY2) / 2;
                    // distance to the MCP centre
                    timedata dist = timedata(sqrt(pow(dTX, 2) + pow(dTY, 2)));

                    // time difference between events on X and Y
                    // its a curious way of calculating an absolute value
                    // but we want to make sure that we don't loose precision (time coded on 64bits)
                    timedata dT = AbsDiff(TX, TY);

                    // an atom would be inside the MCP radius and fall atoms the same time on X and Y
                    if (dist < MCPdiameter && dT < deltaT)
                    { // Once we KNOW that the possible atom is on the MCP, we compute its offset value and
                        // compare it to the MCP offset reference map.

                        // time difference between events on X and Y
                        // its a curious way of calculating an absolute value
                        // but we want to make sure that we don't loose precision (time coded on 64bits)
                        timedata S = TX1 + TX2 - TY1 - TY2;

                        int X = 708 - TX1 + TX2;
                        int Y = 708 - TY1 + TY2;
                        if (AbsDiff(offset_p[X][Y], S) < 5) // 5 ??? --> this need to be set using resolution map.
                        /*if (AbsDiff(S, 0) < 80)*/
                        {
                            atoms_p.push_back(atomdata{TX1, TX2, TY1, TY2});
                            atomfound = true;
                        }
                    }
                    ++searchY2;
                }
                ++searchY1;
            }
            ++searchX2;
        }
        if (atomfound)
        {
            X2_p.erase(prev(searchX2));
            Y1_p.erase(prev(searchY1));
            Y2_p.erase(prev(searchY2));
        }

        X1_p.erase(X1_p.begin());
    }
    if (atoms_p.empty())
    {
        return false;
    }
    else
        return true;
}

bool reconstruction3(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params,
                     std::vector<std::vector<int>> &offset_p)
{
    cout << "Reconstruction 3 : Recovering All Potential Atoms" << endl;
    // gate X, gatY->bulb width
    timedata gateY = 2 * timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = 2 * timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);

    while (X1_p.begin() != X1_p.end())
    {
        // We first definitively get rid of all events occurring before the bulb start on X2, Y1, Y2
        // By construction, they will never match with a later event on X1
        if (*X1_p.begin() > gateX)
            while (X2_p.begin() != X2_p.end() && *X2_p.begin() < (*X1_p.begin() - gateX))
                X2_p.erase(X2_p.begin());
        if (*X1_p.begin() > gateY)
        {
            while (Y1_p.begin() != Y1_p.end() && *Y1_p.begin() < (*X1_p.begin() - gateY))
                Y1_p.erase(Y1_p.begin());
            while (Y2_p.begin() != Y2_p.end() && *Y2_p.begin() < min((*X1_p.begin() - gateY), *X1_p.begin()))
                Y2_p.erase(Y2_p.begin());
        }
        // We then look for events that can correspond to an atom
        // As soon as we have found such events, we erase them from the list and keep going along X1
        auto searchX2 = X2_p.begin();
        // cout << "#######################" << endl;
        // cout << "X1" << *X1_p.begin() << endl;
        while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX))
        {
            if (*X1_p.begin() == 83123996)
            {
                cout << "Salut " << *searchX2 << endl;
            };
            auto searchY1 = Y1_p.begin();
            while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY))
            {
                auto searchY2 = Y2_p.begin();
                while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY))
                {
                    timedata TX1 = *X1_p.begin();
                    timedata TX2 = *searchX2;
                    timedata TY1 = *searchY1;
                    timedata TY2 = *searchY2;
                    timedata dTX = AbsDiff(TX1, TX2);
                    timedata dTY = AbsDiff(TY1, TY2);
                    timedata TX = (TX1 + TX2) / 2;
                    timedata TY = (TY1 + TY2) / 2;
                    // distance to the MCP centre
                    timedata dist = timedata(sqrt(pow(dTX, 2) + pow(dTY, 2)));

                    // time difference between events on X and Y
                    // its a curious way of calculating an absolute value
                    // but we want to make sure that we don't loose precision (time coded on 64bits)
                    timedata dT = AbsDiff(TX, TY);

                    // an atom would be inside the MCP radius and fall atoms the same time on X and Y
                    if (dist < MCPdiameter && dT < deltaT)
                    { // Once we KNOW that the possible atom is on the MCP, we compute its offset value and
                        // compare it to the MCP offset reference map.

                        // time difference between events on X and Y
                        // its a curious way of calculating an absolute value
                        // but we want to make sure that we don't loose precision (time coded on 64bits)
                        timedata S = TX1 + TX2 - TY1 - TY2;

                        int X = 708 - TX1 + TX2;
                        int Y = 708 - TY1 + TY2;
                        if (AbsDiff(offset_p[X][Y], S) < 5) // 5 ??? --> this need to be set using resolution map.
                        /*if (AbsDiff(S, 0) < 80)*/
                        {
                            atoms_p.push_back(atomdata{TX1, TX2, TY1, TY2});
                        }
                    }
                    ++searchY2;
                }
                ++searchY1;
            }
            ++searchX2;
        }
        X1_p.erase(X1_p.begin());
    }
    if (atoms_p.empty())
    {
        return false;
    }
    else
        return true;
}

bool reconstruction4(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params)
{
    cout << "Reconstruction 4 : Recovering All Potential Atoms Without Any Offset Filter" << endl;
    // gate X, gatY->bulb width
    timedata gateY = timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);

    while (X1_p.begin() != X1_p.end())
    {
        // We first definitively get rid of all events occurring before the bulb start on X2, Y1, Y2
        // By construction, they will never match with a later event on X1
        if (*X1_p.begin() > gateX)
            while (X2_p.begin() != X2_p.end() && *X2_p.begin() < (*X1_p.begin() - gateX))
                X2_p.erase(X2_p.begin());
        if (*X1_p.begin() > gateY)
        {
            while (Y1_p.begin() != Y1_p.end() && *Y1_p.begin() < (*X1_p.begin() - gateY))
                Y1_p.erase(Y1_p.begin());
            while (Y2_p.begin() != Y2_p.end() && *Y2_p.begin() < min((*X1_p.begin() - gateY), *X1_p.begin()))
                Y2_p.erase(Y2_p.begin());
        }
        // We then look for events that can correspond to an atom
        // As soon as we have found such events, we erase them from the list and keep going along X1

        auto searchX2 = X2_p.begin();
        while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX))
        {
            auto searchY1 = Y1_p.begin();
            while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY))
            {
                auto searchY2 = Y2_p.begin();
                while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY))
                {
                    timedata TX1 = *X1_p.begin();
                    timedata TX2 = *searchX2;
                    timedata TY1 = *searchY1;
                    timedata TY2 = *searchY2;
                    timedata dTX = AbsDiff(TX1, TX2);
                    timedata dTY = AbsDiff(TY1, TY2);
                    timedata TX = (TX1 + TX2) / 2;
                    timedata TY = (TY1 + TY2) / 2;
                    // distance to the MCP centre
                    timedata dist = timedata(sqrt(pow(dTX, 2) + pow(dTY, 2)));

                    // time difference between events on X and Y
                    // its a curious way of calculating an absolute value
                    // but we want to make sure that we don't loose precision (time coded on 64bits)
                    timedata dT = AbsDiff(TX, TY);

                    // an atom would be inside the MCP radius and fall atoms the same time on X and Y
                    if (dist < MCPdiameter && dT < deltaT)
                    { // Once we KNOW that the possible atom is on the MCP, we compute its offset value and
                        // compare it to the MCP offset reference map.
                        if (dT > deltaT)
                        {
                            cout << "dT < deltaT" << endl;
                        }
                        // time difference between events on X and Y
                        // its a curious way of calculating an absolute value
                        // but we want to make sure that we don't loose precision (time coded on 64bits)
                        timedata S = TX1 + TX2 - TY1 - TY2;

                        int X = 708 - TX1 + TX2;
                        int Y = 708 - TY1 + TY2;
                        atoms_p.push_back(atomdata{TX1, TX2, TY1, TY2});
                    }
                    ++searchY2;
                }
                ++searchY1;
            }
            ++searchX2;
        }
        X1_p.erase(X1_p.begin());
    }
    if (atoms_p.empty())
    {
        return false;
    }
    else
        return true;
}

bool reconstruction5(list<timedata> &X1_p,
                     list<timedata> &X2_p,
                     list<timedata> &Y1_p,
                     list<timedata> &Y2_p,
                     vector<atomdata> &atoms_p,
                     paramstruct params)
{
    cout << "Reconstruction 5 : Recovering Only Isolated Atoms for very dilute clouds." << endl;
    // gate X, gatY->bulb width
    timedata gateY = timedata((params.evgate + params.atgate) / params.res);
    timedata gateX = timedata((params.evgate + params.atgate) / params.res);
    // MCPdiameter -> events can only correspond to an atom if they can be traced back to a position inside the MCP radius
    timedata MCPdiameter = timedata(params.evgate / params.res);
    // deltaT -> events can only correspond to an atom if the times on X and Y are close to each other
    timedata deltaT = timedata(params.deltaTgate / params.res);
    auto atoms_alone = 0;
    auto atoms_duplicata = 0;

    while (X1_p.begin() != X1_p.end())
    {
        // We first definitively get rid of all events occurring before the bulb start on X2, Y1, Y2
        // By construction, they will never match with a later event on X1
        if (*X1_p.begin() > gateX)
            while (X2_p.begin() != X2_p.end() && *X2_p.begin() < (*X1_p.begin() - gateX))
                X2_p.erase(X2_p.begin());
        if (*X1_p.begin() > gateY)
        {
            while (Y1_p.begin() != Y1_p.end() && *Y1_p.begin() < (*X1_p.begin() - gateY))
                Y1_p.erase(Y1_p.begin());
            while (Y2_p.begin() != Y2_p.end() && *Y2_p.begin() < min((*X1_p.begin() - gateY), *X1_p.begin()))
                Y2_p.erase(Y2_p.begin());
        }

        // This part is specific to this reconstruction : we want to make sure that there are no other candidate that could be reconstructed
        auto count = 0;
        // We then look for events that can correspond to an atom
        // As soon as we have found such events, we erase them from the list and keep going along X1
        atomdata atom;
        timedata TX1 = 0;
        timedata TX2 = 0;
        timedata TY1 = 0;
        timedata TY2 = 0;
        auto searchX2 = X2_p.begin();
        bool atomfound = false;
        while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX))
        {
            auto searchY1 = Y1_p.begin();
            while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY))
            {
                auto searchY2 = Y2_p.begin();
                while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY))
                {
                    count++;
                    timedata TX1 = *X1_p.begin();
                    timedata TX2 = *searchX2;
                    timedata TY1 = *searchY1;
                    timedata TY2 = *searchY2;
                    timedata dTX = AbsDiff(TX1, TX2);
                    timedata dTY = AbsDiff(TY1, TY2);
                    timedata TX = (TX1 + TX2) / 2;
                    timedata TY = (TY1 + TY2) / 2;
                    // distance to the MCP centre
                    timedata dist = timedata(sqrt(pow(dTX, 2) + pow(dTY, 2)));

                    // time difference between events on X and Y
                    // its a curious way of calculating an absolute value
                    // but we want to make sure that we don't loose precision (time coded on 64bits)
                    timedata dT = AbsDiff(TX, TY);

                    // an atom would be inside the MCP radius and fall atoms the same time on X and Y
                    if (dist < MCPdiameter && dT < deltaT)
                    { // Once we KNOW that the possible atom is on the MCP, we compute its offset value and
                        // compare it to the MCP offset reference map.
                        if (dT > deltaT)
                        {
                            cout << "dT < deltaT" << endl;
                        }
                        // time difference between events on X and Y
                        // its a curious way of calculating an absolute value
                        // but we want to make sure that we don't loose precision (time coded on 64bits)
                        timedata S = TX1 + TX2 - TY1 - TY2;
                        atomfound = true;
                        int X = 708 - TX1 + TX2;
                        int Y = 708 - TY1 + TY2;
                        atom = atomdata{TX1, TX2, TY1, TY2};
                    }
                    ++searchY2;
                }
                ++searchY1;
            }
            ++searchX2;
        }
        if (atomfound)
        {
            if (count == 1)
            {
                atoms_p.push_back(atom);
                atoms_alone++;
            }
            else if (count > 1)
            {
                atoms_duplicata = atoms_duplicata + count;
            }
        }

        X1_p.erase(X1_p.begin());
    }
    cout << "Nombre d'atomes reconstruits : " << atoms_alone << endl;
    cout << "Nombre d'atomes non pris en compte : " << atoms_duplicata << endl;
    if (atoms_p.empty())
    {
        return false;
    }
    else
        return true;
}

/*
======================================
AUTRES FONCTIONS APPELÉES DANS LE CODE
======================================
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

void WriteStatistics(string filename, map<string, float> &reconstruction_statistics)
// This functions save on the disk the statistics of the reconstruction programm.
//
// Arguments
// ---------
// filename : string,
//		name of the file one wants to save. It should be "path/to/directory/name". The extension of the file will be .stats
// reconstruction_statistics : map (kind of dictionnary) of string and float
//		the statistics of the code, for exemple ["reconstruction time":1.23 , "X2 lenght": 1023, etc... ]
// 		first must be string, second float.

{
    ofstream statfile;
    statfile.open(filename);
    map<string, float>::iterator iter;
    for (iter = reconstruction_statistics.begin(); iter != reconstruction_statistics.end(); iter++)
    {
        statfile << (*iter).first << " : " << (*iter).second << "\n";
        cout << (*iter).first << " ---> " << (*iter).second << "\n";
    }

    statfile.close();
}

void WriteAtoms(string filename, vector<atomdata> &atoms_p)
/*
filename : of the form /path/where/the/code/save/the/file WITHOUT the extension
atoms : vecteur de atomdata aka 4 temps x1, x2, y1, y2
*/
{

    ofstream file;
    file.open(filename, ofstream::out | ofstream::binary);

    ostringstream message;

    // Serialize the atom data
    std::vector<timedata> buffer;
    buffer.reserve(4 * atoms_p.size());
    for (auto i = atoms_p.begin(); i != atoms_p.end(); ++i)
    {
        buffer.push_back((*i).TX1);
        buffer.push_back((*i).TX2);
        buffer.push_back((*i).TY1);
        buffer.push_back((*i).TY2);
    }

    file.write(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(timedata));

    file.close();
    // WriteLog(message.str());
}

void WriteOffsets(string filename, list<int> &offset_values_p)
// This function write the pointer list offset_value into a file name file.
{

    ofstream file;
    file.open(filename, ofstream::out | ofstream::binary);
    ostringstream message;
    std::vector<timedata> buffer;
    buffer.reserve(offset_values_p.size());
    std::copy(begin(offset_values_p), end(offset_values_p), back_inserter(buffer));

    file.write(reinterpret_cast<const char *>(buffer.data()), buffer.size() * sizeof(timedata));
    file.close();
}

void copy_paste_file(string input_file, string output_file)
/*Copy the file input_file and paste into output_file*/
{
    std::ifstream path(input_file);
    std::ofstream writer(output_file);
    if (path)
    {
        if (writer)
        {
            path.seekg(0, path.end);
            long sizee = path.tellg();
            path.seekg(0);

            char *buffer = new char[sizee];

            path.read(buffer, sizee);
            writer.write(buffer, sizee);

            delete[] buffer;
            std::cout << "[Copy] " << input_file << " into " << output_file << std::endl;
            path.close();
            writer.close();
        }
        else
        {
            std::cout << "[WARNING] Impossible to write in " << output_file << std::endl;
        }
    }
    else
    {
        std::cout << "[WARNING] Impossible to read " << input_file << std::endl;
    }
}

list<int> get_offset_of_atoms(vector<atomdata> &atoms_p, std::vector<std::vector<int>> &offset_p)
{
    list<int> offset_of_atoms;
    for (auto vectorit = atoms_p.begin(); vectorit != atoms_p.end(); ++vectorit)
    {
        // cout << (*vectorit).TX1 << endl;
        int S = (*vectorit).TX1 + (*vectorit).TX2 - (*vectorit).TY1 - (*vectorit).TY2;
        int X = 708 - (*vectorit).TX1 + (*vectorit).TX2;
        int Y = 708 - (*vectorit).TY1 + (*vectorit).TY2;
        int offset_diff = S - offset_p[X][Y];
        offset_of_atoms.push_back(offset_diff);
    }
    return offset_of_atoms;
}
