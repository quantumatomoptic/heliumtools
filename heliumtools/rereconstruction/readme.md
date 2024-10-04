# From time signal to atoms : the rereconstruction module

The  [rereconstruction program](../heliumtools/rereconstruction/rereconstruction.cpp) ables one to reconstruct *a posteriori* a sequence if we saved all the raw datas from the experiment. It can be run on a linux computer. In order to launch the program, one must just change directory towards the location of the code, configure the reconstruction parameters (see [below](#choosing-the-rereconstruction-parameters)) and enter :
```
./rereconstruction
```
You can also modify the rereconstruction program. To do so, edit the `.cpp` file and then compile it. For example:
```
g++ rereconstruction.cpp -o rereconstruction
```
and then just run your brand-new program. Of course, do not forget to commit your modification, comment your code and update this readme file!

### Choosing the rereconstruction parameters
When it is run, the rereconstruction program reads the following files 
* [Offset.txt](../heliumtools/rereconstruction/Offset.txt) : the offset map of the experiment,
* [conf.txt](../heliumtools/rereconstruction/conf.txt) : the parameter file for the rereconstruction process.

The configuration file must have the following lines :
* Line 1: The sequence folder we want to re-reconstruct : for example `/mnt/manip_E/2022/11/23/003`,
* Line 2: The folder in which we want to create the .atoms files. This folder should not exists.
* Line 3: The reconstruction program number we want to use. See the available reconstruction program below.
* Line 4: The maximum deviation allowed for the offset from the offset map.
* Line 5: The offset we want to add to ALL offset difference values. 
            formula for the offset difference is
                    off_diff = S - offset_p[X][Y] + params.offset_offset;


The available reconstructions are the following
* Reconstruction 1 : Old reconstruction program [Ziyad style](#reconstruction1-program)
* Reconstruction 2 : New reconstruction program. Recover all potential atoms (calling reconstruction3) and then sort atoms by their offset_diff value and keep only the lowest offset if two potential atoms have a "same" value in their column.
* Reconstruction 3 : Recovering All Potential Atoms taking into account offset maps.
* Reconstruction 4 : Recovering Only Isolated Atoms for very dilute clouds --> use this to perform offset maps.


## Reconstruction programs 
### Reconstruction1 program
This is the reconstruction used since the new MCP electronics, which is the one we used until July 2022. The offset test was just: is the offset lower than 10. The 4 lists X1, X2, Y1, Y2 are ordered times and α, β, γ, and δ are respectively elements of each list. \
T is the temporal diameter of the MCP: T = 80 ns (see [Q. Marolleau PhD thesis](https://theses.fr/2022UPASP166) for more details).
```cython
type X1, X2, Y1, Y2  = lists of ordered times
atoms = [] // empty list
while X1.is_not_empty():
    α = X1.begin()
    {X2, Y1, Y2}.suppress_all_events_before(α – T); //they will never be used since X1 is time ordered
    selected_events = {[α], X2, Y1, Y2}.select_all_events_before(α+T);
    for each quadruplet Q = {α, β, γ, δ} in selected_events:
	    test1 = is_the_atom_on_MCP(Q);
	    test2 = is_offset_zero(Q); // check the atom offset
        if test1 & test2 :
            atoms.append({α, β, γ, δ});
            X2.remove(β); Y1.remove(γ); Y2.remove(δ);
            break
    X1.remove(α); 
```

NB : sur le code ci-dessus, il est écrit qu'on test tous les quadruplets Q =  {α, β, γ, δ} de la liste des évennements sélectionnés. Ce n'est en fait pas vrai (regarder les codes en détail pour une plus grande explication). En effet, on initialise le début de la liste Y1 et Y2 en dehors de la boucle X2 donc on va bien visiter plusieurs couples de Y1 et Y2 pour la première valeur de X2 mais seulement un seul pour les autres valeurs de X2.
Dans le code, concrètement, la boucle est ainsi : 
```
auto searchX2 = X2_p.begin(), searchY1 = Y1_p.begin(), searchY2 = Y2_p.begin();
while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX) && !atomfound)
{
    while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY) && !atomfound)
    {
        while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY) && !atomfound)
        {...}}}
```
alors que si on avait voulu explorer TOUS les quadruplets Q, alors il aurait fallu écrire :
```
auto searchX2 = X2_p.begin();
while (searchX2 != X2_p.end() && *searchX2 < (*X1_p.begin() + gateX) && !atomfound)
{	
	auto  searchY1 = Y1_p.begin();
    while (searchY1 != Y1_p.end() && *searchY1 < (*X1_p.begin() + gateY) && !atomfound)
    {	
    	auto searchY2 = Y2_p.begin();
        while (searchY2 != Y2_p.end() && *searchY2 < (*X1_p.begin() + gateY) && !atomfound)
        {...}}}
```
### Reconstruction2 program
This is the reconstruction used between July 2022 and February 2023. It follows the same principle as the program above, adding a criterion on the offset value using the offset maps made in July 2022. I think this program doesn't really have any interest and it's only here to keep track of the program. 
### Reconstruction3 program 
Unlike the previous ones, this reconstruction *tests all possible quadruplets*, and as soon as a quadruplet meets the offset criterion, it adds it to the list of atoms.
```cython
type X1, X2, Y1, Y2  = lists // of ordered times 
atoms = [] // empty list
T = 1.3 * MCP_diameter // take a security margin
while X1.is_not_empty():
    α = X1.begin()
    {X2, Y1, Y2}.suppress_all_events_before(α – T); //they will never be used since X1 is time ordered
    selected_events = {[α], X2, Y1, Y2}.select_all_events_before(α+T);
    for each quadruplet Q = {α, β, γ, δ} in selected_events:
	    test1 = is_the_atom_on_MCP(Q); // the position of the detected atom must be on the MCP
	    test2 = is_offset_zero_using_offset_map(Q); // we use offset map now
        if test1 & test2 :
            atoms.append({α, β, γ, δ});
            // X2.remove(β); Y1.remove(γ); Y2.remove(δ); --> we do not remove any X2, Y1 and Y2
            // break we do not break but keep for looking atoms
    X1.remove(α); 
```
### Reconstruction4 program
This reconstruction is a copy of the number 3 BUT it does not test the offset so it reconstructs A LOT of atoms. 

### Reconstruction5 program
This code was added on July, 20th of 2023.
In this code, we make sure to reconstruct an atom if and only if he is the only candidates. This means that no signal on any channel should have been recorded between the signal time and ±90 ns.

$$
\forall Z_J \in \{X_1, X_2, Y_1, Y_2\}, \quad \forall \alpha \in Z_J, \left| \alpha - Z_J \right| < 90 \, \text{ns}
$$


This code is usefull to acquire an offset map as we do not mix atoms. 


