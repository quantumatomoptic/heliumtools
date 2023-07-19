# Code de rereconstruction
Ce code s'utilise avec le fichier de configuration `conf.txt` dans lequel sont remplis les paramètres suivants :
* ligne 1 : le dossier de séquence que nous voulons re-reconstruire,
* ligne 2 : le dossier dans lequel on veut créer les .atoms
* ligne 3 : le numéro du programme que l'on veut utiliser


In the following we give some explanations about how those codes work. 


## Reconstruction 1
Il s'agit de la reconstruction utilisée depuis la nouvelle électronique du MCP, c'est celui-ci que nous utilisions jusque Juillet 2022.
The 4 lists X1, X2, Y1, Y2 are ordored times and α, β, γ and δ are respectively elements of each list.
T is the temporal diameter of the MCP : T = 80 ns (see Q. Marolleau PhD thesis for more details).
```
type X1, X2, Y1, Y2  = lists of ordered times
atoms = [] // empty list
while X1.is_not_empty():
    α = X1.begin()
    {X2, Y1, Y2}.suppress_all_events_before(α – T); //they will never be used since X1 is time orderd
    selected_events = {[α], X2, Y1, Y2}.select_all_events_before(α+T);
    for each quadruplet Q = {α, β, γ, δ} in selected_events:
	    test1 = is_the_atom_on_MCP(Q);
	    test2 = is_offset_zero(Q);
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

## Reconstruction 2
Il s'agit de la reconstruction utilisée entre juillet 2022 et février 2023. Celle-ci reprend le même principe que le programme ci-dessus en rajoutant un critère sur la valeur d'offset en utilisant les cartes d'offset réalisées en juillet 2022. Je pense que ce programme n'a pas vraiment d'intérêt et il n'est là que pour ne pas perdre la trace du programme. 

## Reconstruction 3
Contrairement aux précédentes, cette reconstruction test tous les quadruplets possibles et dès qu'un quadruplet valide le critère d'offset, elle l'ajoute à la liste des atomes.
```
type X1, X2, Y1, Y2  = lists of ordered times 
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


## Reconstruction 4
This reconstruction is a copy of the number 3 BUT it does not test the offset so it reconstructs A LOT of atoms. 


