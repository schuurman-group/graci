"""
Solvent information
"""

# Dielectric constants for common solvents
dielectric = {
    'water'                           : 78.3553,
    'acetonitrile'                    : 35.6880,
    'methanol'                        : 32.6130,
    'ethanol'                         : 24.8520,
    'isoquinoline'                    : 11.0000,
    'quinoline'                       :  9.1600,
    'chloroform'                      :  4.7113,
    'diethylether'                    :  4.2400,
    'dichloromethane'                 :  8.9300,
    'dichloroethane'                  : 10.1250,
    'carbontetrachloride'             :  2.2280,
    'benzene'                         :  2.2706,
    'toluene'                         :  2.3741,
    'chlorobenzene'                   :  5.6968,
    'nitromethane'                    : 36.5620,
    'heptane'                         :  1.9113,
    'cyclohexane'                     :  2.0165,
    'aniline'                         :  6.8882,
    'acetone'                         : 20.4930,
    'tetrahydrofuran'                 :  7.4257,
    'dimethylsulfoxide'               : 46.8260,
    'argon'                           :  1.4300,
    'krypton'                         :  1.5190,
    'xenon'                           :  1.7060,
    'n-octanol'                       :  9.8629,
    '1,1,1-trichloroethane'           :  7.0826,
    '1,1,2-trichloroethane'           :  7.1937,
    '1,2,4-trimethylbenzene'          :  2.3653,
    '1,2-dibromoethane'               :  4.9313,
    '1,2-ethanediol'                  : 40.2450,
    '1,4-dioxane'                     :  2.2099,
    '1-bromo-2-methylpropane'         :  7.7792,
    '1-bromooctane'                   :  5.0244,
    '1-bromopentane'                  :  6.2690,
    '1-bromopropane'                  :  8.0496,
    '1-butanol'                       : 17.3320,
    '1-chlorohexane'                  :  5.9491,
    '1-chloropentane'                 :  6.5022,
    '1-chloropropane'                 :  8.3548,
    '1-decanol'                       :  7.5305,
    '1-fluorooctane'                  :  3.8900,
    '1-heptanol'                      : 11.3210,
    '1-hexanol'                       : 12.5100,
    '1-hexene'                        :  2.0717,
    '1-hexyne'                        :  2.6150,
    '1-iodobutane'                    :  6.1730,
    '1-iodohexadecane'                :  3.5338,
    '1-iodopentane'                   :  5.6973,
    '1-iodopropane'                   :  6.9626,
    '1-nitropropane'                  : 23.7300,
    '1-nonanol'                       :  8.5991,
    '1-pentanol'                      : 15.1300,
    '1-pentene'                       :  1.9905,
    '1-propanol'                      : 20.5240,
    '2,2,2-trifluoroethanol'          : 26.7260,
    '2,2,4-trimethylpentane'          :  1.9358,
    '2,4-dimethylpentane'             :  1.8939,
    '2,4-dimethylpyridine'            :  9.4176,
    '2,6-dimethylpyridine'            :  7.1735,
    '2-bromopropane'                  :  9.3610,
    '2-butanol'                       : 15.9440,
    '2-chlorobutane'                  :  8.3930,
    '2-heptanone'                     : 11.6580,
    '2-hexanone'                      : 14.1360,
    '2-methoxyethanol'                : 17.2000,
    '2-methyl-1-propanol'             : 16.7770,
    '2-methyl-2-propanol'             : 12.4700,
    '2-methylpentane'                 :  1.8900,
    '2-methylpyridine'                :  9.9533,
    '2-nitropropane'                  : 25.6540,
    '2-octanone'                      :  9.4678,
    '2-pentanone'                     : 15.2000,
    '2-propanol'                      : 19.2640,
    '2-propen-1-ol'                   : 19.0110,
    '3-methylpyridine'                : 11.6450,
    '3-pentanone'                     : 16.7800,
    '4-heptanone'                     : 12.2570,
    '4-methyl-2-pentanone'            : 12.8870,
    '4-methylpyridine'                : 11.9570,
    '5-nonanone'                      : 10.6000,
    'aceticacid'                      :  6.2528,
    'acetophenone'                    : 17.4400,
    'a-chlorotoluene'                 :  6.7175,
    'anisole'                         :  4.2247,
    'benzaldehyde'                    : 18.2200,
    'benzonitrile'                    : 25.5920,
    'benzylalcohol'                   : 12.4570,
    'bromobenzene'                    :  5.3954,
    'bromoethane'                     :  9.0100,
    'bromoform'                       :  4.2488,
    'butanal'                         : 13.4500,
    'butanoicacid'                    :  2.9931,
    'butanone'                        : 18.2460,
    'butanonitrile'                   : 24.2910,
    'butylamine'                      :  4.6178,
    'butylethanoate'                  :  4.9941,
    'carbondisulfide'                 :  2.6105,
    'cis-1,2-dimethylcyclohexane'     :  2.0600,
    'cis-decalin'                     :  2.2139,
    'cyclohexanone'                   : 15.6190,
    'cyclopentane'                    :  1.9608,
    'cyclopentanol'                   : 16.9890,
    'cyclopentanone'                  : 13.5800,
    'decalin-mixture'                 :  2.1960,
    'dibromomethane'                  :  7.2273,
    'dibutylether'                    :  3.0473,
    'diethylamine'                    :  3.5766,
    'diethylsulfide'                  :  5.7230,
    'diiodomethane'                   :  5.3200,
    'diisopropylether'                :  3.3800,
    'dimethyldisulfide'               :  9.6000,
    'diphenylether'                   :  3.7300,
    'dipropylamine'                   :  2.9112,
    'e-1,2-dichloroethene'            :  2.1400,
    'e-2-pentene'                     :  2.0510,
    'ethanethiol'                     :  6.6670,
    'ethylbenzene'                    :  2.4339,
    'ethylethanoate'                  :  5.9867,
    'ethylmethanoate'                 :  8.3310,
    'ethylphenylether'                :  4.1797,
    'fluorobenzene'                   :  5.4200,
    'formamide'                       : 108.940,
    'formicacid'                      : 51.1000,
    'hexanoicacid'                    :  2.6000,
    'iodobenzene'                     :  4.5470,
    'iodoethane'                      :  7.6177,
    'iodomethane'                     :  6.8650,
    'isopropylbenzene'                :  2.3712,
    'm-cresol'                        : 12.4400,
    'mesitylene'                      :  2.2650,
    'methylbenzoate'                  :  6.7367,
    'methylbutanoate'                 :  5.5607,
    'methylcyclohexane'               :  2.0240,
    'methylethanoate'                 :  6.8615,
    'methylmethanoate'                :  8.8377,
    'methylpropanoate'                :  6.0777,
    'm-xylene'                        :  2.3478,
    'n-butylbenzene'                  :  2.3600,
    'n-decane'                        :  1.9846,
    'n-dodecane'                      :  2.0060,
    'n-hexadecane'                    :  2.0402,
    'n-hexane'                        :  1.8819,
    'nitrobenzene'                    : 34.8090,
    'nitroethane'                     : 28.2900,
    'n-methylaniline'                 :  5.9600,
    'n-methylformamide-mixture'       : 181.560,
    'n,n-dimethylacetamide'           : 37.7810,
    'n,n-dimethylformamide'           : 37.2190,
    'n-nonane'                        :  1.9605,
    'n-octane'                        :  1.9406,
    'n-pentadecane'                   :  2.0333,
    'n-pentane'                       :  1.8371,
    'n-undecane'                      :  1.9910,
    'o-chlorotoluene'                 :  4.6331,
    'o-cresol'                        :  6.7600,
    'o-dichlorobenzene'               :  9.9949,
    'o-nitrotoluene'                  : 25.6690,
    'o-xylene'                        :  2.5454,
    'pentanal'                        : 10.0000,
    'pentanoicacid'                   :  2.6924,
    'pentylamine'                     :  4.2010,
    'pentylethanoate'                 :  4.7297,
    'perfluorobenzene'                :  2.0290,
    'p-isopropyltoluene'              :  2.2322,
    'propanal'                        : 18.5000,
    'propanoicacid'                   :  3.4400,
    'propanonitrile'                  : 29.3240,
    'propylamine'                     :  4.9912,
    'propylethanoate'                 :  5.5205,
    'p-xylene'                        :  2.2705,
    'pyridine'                        : 12.9780,
    'sec-butylbenzene'                :  2.3446,
    'tert-butylbenzene'               :  2.3447,
    'tetrachloroethene'               :  2.2680,
    'tetrahydrothiophene-s,s-dioxide' : 43.9620,
    'tetralin'                        :  2.7710,
    'thiophene'                       :  2.7270,
    'thiophenol'                      :  4.2728,
    'trans-decalin'                   :  2.1781,
    'tributylphosphate'               :  8.1781,
    'trichloroethene'                 :  3.4220,
    'triethylamine'                   :  2.3832,
    'xylene-mixture'                  :  2.3879,
    'z-1,2-dichloroethene'            :  9.2000
}
