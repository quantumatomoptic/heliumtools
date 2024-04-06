
The default fit class of HAL was copied and pasted in the heliumtools repository so that any user is not compelled to recode a fit function. In the present page, we show how to use the fit class inspired by HAL fit class. 

#### Exemple
The following code fit a damped oscillation
```python3
import numpy as np
import matplotlib.pyplot as plt
from heliumtools.fit.dampedoscillation import DampedOscillation1DFit


fit = DampedOscillation1DFit(x = x_data, z = z_data)
fit.x_unit = "a.u." #unité de x
fit.z_unit = "a.u." #unité de z
fit.do_guess()  # guess data
#fit.guess = [ ... ]
fit.do_fit()
fit.compute_values()
text = fit.return_result_as_string()
print(text)

yth = fit.eval(x_data)
plt.plot(x_data, z_data, ".", label="exp")
plt.plot(x_data, yth, label="fit")
plt.xlabel(f"x axis ({fit.x_unit})")
plt.ylabel(f"y axis ({fit.z_unit})")
plt.legend()
plt.show()
```



#### Arguments de la classe fit
*x* = [datas along x], default is [ ]\
*z* = [datas along z], default is [ ]\
*guess* : list de la taille du nombre d'argument que prend _fitfunc, None par défaut.\
*x_unit* and *z_unit* : string, unité des de x et z.\
*formula_help* : un string qui donne la forme des paramètres.\
*parameters_help* : un string qui donne le nom des paramètres. \
*name* : le nom du fit\
*popt* : paramètres optimaux\
*pcov* : matrice de covariance\
*perr* : erreur sur les paramètres\

#### Methods of the fit class
*_fitfunc* : la fonction qui définie le fit : prend en argument  \
*do_fit* : pour faire le fit . Si le guess p0 n'est pas donné,c'est par la routine de fit.Any keyword option is passed to the fitting routine scipy.optimize.curve_fit\
*do_guess* : faire un guess automatique. \
*eval(x_data)* : évaluer sur x_data en utilisant les résultats du fit. \
*compute_values* : calculer les values (quantités physiques).\
*return_result_as_string* : retourne un grand string pour printer les résultat du fit comme c'est fait dans HAL. \

