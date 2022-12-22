"""
Class object to hold information germane to bitci 
wavefunctions
"""
import copy as copy

class Bitciwfn:
    """Class Bitciwfn holds configuration and file name
       information required to keep track of a bitci
       wavefunction"""
    def __init__(self):
        # class variables
        # List of numbers of configurations per irrep
        self.nconf         = {'adiabatic' : None, 'diabatic' : None}
        # bitci configuration scratch file names
        self.conf_name     = {'adiabatic' : None, 'diabatic' : None}
        # bitci configuration scratch file numbers
        self.conf_units    = {'adiabatic' : None, 'diabatic' : None}
        # Lists of bitci eigenvector scratch file names (one per irrep)
        # for various representations
        self.ci_name       = {'adiabatic' : None, 'diabatic' : None}
        # List of bitci eigenvector scratch file numbers (one per irrep)
        # for various representations
        self.ci_units      = {'adiabatic' : None, 'diabatic' : None}
        # List of bitci Q-space energy correction scratch file numbers
        # (one per irrep)
        self.eq_units      = None
        # class label
        self.label         = 'default'

    def copy(self):
        """create of deepcopy of self"""
        new = Bitciwfn()

        var_dict = {key:value for key,value in self.__dict__.items()
                   if not key.startswith('__') and not callable(key)}

        for key, value in var_dict.items():
            setattr(new, key, copy.deepcopy(value))

        return new

    #
    def set_nconf(self, nconf, rep='adiabatic'):
        """Sets the numbers of configurations"""
        self.nconf[rep] = nconf
        return

    #
    def set_confunits(self, confunits, rep='adiabatic'):
        """Sets the bitci configuration scratch file numbers"""
        self.conf_units[rep] = confunits
        return

    #
    def set_ciunits(self, ciunits, rep='adiabatic'):
        """Adds the list of bitci eigenvector scratch file numbers"""
        self.ci_units[rep] = ciunits
        return

    #
    def set_confname(self, confname, rep='adiabatic'):
        """Sets the bitci configuration scratch file names"""
        self.conf_name[rep] = confname
        return

    #
    def set_ciname(self, ciname, rep='adiabatic'):
        """Adds the list of bitci eigenvector scratch file names"""
        self.ci_name[rep] = ciname
        return

    #
    def set_equnits(self, equnits):
        """Adds the list of bitci Q-space energy correction scratch
        file numbers"""
        self.eq_units = equnits
        return
