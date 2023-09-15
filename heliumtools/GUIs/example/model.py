import numpy as np

from heliumtools.GUIs.base.parameter_model import GlobalModelForParameters


def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * standard_deviation**2))


class Model(GlobalModelForParameters):
    """
    The model inherits the GlobalModelForParameters class from heliumtools.misc.
    """

    def __init__(self, parameters={}, data=[], metadata={}):
        self.mean = 2
        self.std = 4
        self.amplitude = 3
        self.x_range = [-20, 20]
        # The parameter dictionary given in the initialisation will overwrite
        # other parameters
        super().__init__(parameters=parameters, file=__file__)
