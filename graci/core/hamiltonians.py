"""
DFT/MRCI Hamiltonian information
"""

# aliases used for backwards compatibility with old input files
aliases = {
    'canonical'       : 'abinitio',
    'grimme_standard' : 'grimme',
    'grimme_short'    : 'grimme_short',
    'lyskov_standard' : 'r2016',
    'lyskov_short'    : 'r2016_short',
    'heil17_standard' : 'r2017',
    'heil17_short'    : 'r2017_short',
    'heil18_standard' : 'r2018',
    'heil18_short'    : 'r2018_short',
    'test_exp'        : 'qe8'
}

# references
references = {
    'grimme'        : 'S. Grimme & M. Waletzke, J. Chem. Phys., 111, 5645 (1999)',
    'grimme_short'  : 'S. Grimme & M. Waletzke, J. Chem. Phys., 111, 5645 (1999)',
    'r2016'         : 'I. Lyskov et al., J. Chem. Phys, 144, 034104 (2016)',
    'r2016_short'   : 'I. Lyskov et al., J. Chem. Phys, 144, 034104 (2016)',
    'r2017'         : 'A. Heil & C. Marian, J. Chem. Phys, 147, 194104 (2017)',
    'r2017_short'   : 'A. Heil & C. Marian, J. Chem. Phys, 147, 194104 (2017)',
    'r2018'         : 'A. Heil et al., J. Chem. Phys., 149, 164106 (2018)',
    'r2018_short'   : 'A. Heil et al., J. Chem. Phys., 149, 164106 (2018)',
    'r2022'         : 'D. Dombrowski et al., J. Phys. Chem. A, 127, 2011 (2023)',
    'r2022_short'   : 'D. Dombrowski et al., J. Phys. Chem. A, 127, 2011 (2023)',
    'qe8'           : 'Unpublished',
    'qe8_short'     : 'Unpublished',
    'cvs-qe8'       : 'Unpublished',
    'cvs-qe8_short' : 'Unpublished'
    }

# pretty names for printing
desel  = ' \u03b4E\u209B\u2091\u2097'
pretty = {
    'grimme'        : 'Grimme,' + desel + ' = 1.0',
    'grimme_short'  : 'Grimme,' + desel + ' = 0.8',
    'r2016'         : 'R2016,' + desel + ' = 1.0',
    'r2016_short'   : 'R2016,' + desel + ' = 0.8',
    'r2017'         : 'R2017,' + desel + ' = 1.0',
    'r2017_short'   : 'R2016,' + desel + ' = 0.8',
    'r2018'         : 'R2018,' + desel + ' = 1.0',
    'r2018_short'   : 'R2016,' + desel + ' = 0.8',
    'r2022'         : 'R2022,' + desel + ' = 1.0',
    'r2022_short'   : 'R2022,' + desel + ' = 0.8',
    'qe8'           : 'QE8,' + desel + ' = 1.0',
    'qe8_short'     : 'QE8,' + desel + ' = 0.8',
    'cvs-qe8'       : 'CVS-QE8,' + desel + ' = 1.0',
    'cvs-qe8_short' : 'CVS-QE8,' + desel + ' = 0.8'
}

# Intended use XC functionals
xc_intended = {
    'grimme'        : 'BHLYP',
    'grimme_short'  : 'BHLYP',
    'r2016'         : 'BHLYP',
    'r2016_short'   : 'BHLYP',
    'r2017'         : 'BHLYP',
    'r2017_short'   : 'BHLYP',
    'r2018'         : 'BHLYP',
    'r2018_short'   : 'BHLYP',
    'r2022'         : 'BHLYP',
    'r2022_short'   : 'BHLYP',
    'qe8'           : 'QTP17',
    'qe8_short'     : 'QTP17',
    'cvs-qe8'       : 'QTP17',
    'cvs-qe8_short' : 'QTP17'
}
