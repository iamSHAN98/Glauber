import numpy as np

class Centrality(object) :

    def __init__(self, cent_min, cent_max) :
        self.name = r"%.0f-%.0f %%" % (cent_min, cent_max)
        self.quantity = { 'nch' : [] }

    def check_quantity(self, key) :
        return key in self.quantity

    def add_quantity(self, key) :
        if self.check_quantity(key) :
            raise NameError(key, "already defined")
        self.quantity[key] = [] 

    def add(self, key, val) :
        self.quantity[key].append(val)    

    def clear(self, key) :
        self.quantity[key] = []

    def clear_all(self) :
        for key in self.quantity :
            self.quantity[key] = []

    def erase(self) :
        self.quantity = {}

    def get_name(self) :
        return self.name

    def get(self, key) :
        return self.quantity[key]

    def get_mean(self, key) :
        return np.mean(self.quantity[key])

    def get_maximum(self, key) :
        return np.amax(self.quantity[key])

    def get_minimum(self, key) :
        return np.amin(self.quantity[key])

    def get_error(self, key) :
        return np.std(self.quantity[key])