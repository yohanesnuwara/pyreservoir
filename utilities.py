"""
Basic Utilities for Calculations
@author: Yohanes Nuwara
@email: ign.nuwara97@gmail.com
"""

def convert(x, from_='c', to_='f'):
    """
    Convert metric / SI units to oilfield units
    """

    # temperature
    if from_ == 'c' and to_ == 'f':
        return((x * 9 / 5) + 32)
    if from_ == 'c' and to_ == 'k':
        return(x + 273.15)
    if from_ == 'c' and to_ == 'r':
        return((x * 9 / 5) + 491.67)
    if from_ == 'f' and to_ == 'c':
        return((x - 32) * 5 / 9)
    if from_ == 'f' and to_ == 'k':
        return(((x - 32) * 5 / 9) + 273.15)
    if from_ == 'f' and to_ == 'r':
        return(x + 459.67)
    if from_ == 'k' and to_ == 'c':
        return(x - 273.15)
    if from_ == 'k' and to_ == 'f':
        return((x - 273.15) * (9 / 5) + 32)
    if from_ == 'k' and to_ == 'r':
        return(x * 9 / 5)
    if from_ == 'r' and to_ == 'c':
        return((x - 491.67) * 5 / 9)
    if from_ == 'r' and to_ == 'f':
        return(x - 459.67)
    if from_ == 'r' and to_ == 'k':
        return(x * 5 / 9)

    # pressure
    if from_ == 'atm' and to_ == 'psi':
        return(x * 14.6959)
    if from_ == 'pa' and to_ == 'psi':
        return(x * 0.000145038)
    if from_ == 'bar' and to_ == 'psi':
        return(x * 14.5038)
    if from_ == 'lbf/ft2' and to_ == 'psi':
        return(x / 144)
    if from_ == 'dyne/cm2' and to_ == 'psi':
        return(x * 68947.6)

    # mass
    if from_ == 'kg' and to_ == 'lbm':
        return(x * 2.20462)

    # length
    if from_ == 'm' and to_ == 'ft':
        return(x * 3.28084)
    if from_ == 'mile' and to_ == 'ft':
        return(x * 5280)

    # area
    if from_ == 'm2' and to_ == 'ft2':
        return(x * 10.7639)
    if from_ == 'acre' and to_ == 'ft2':
        return(x * 43560)
    if from_ == 'ha' and to_ == 'ft2':
        return(x * 107639)

    # volume
    if from_ == 'm3' and to_ == 'ft3':
        return(x * 35.3147)
    if from_ == 'acre-ft' and to_ == 'ft3':
        return(x * 43559.9)
    if from_ == 'ft3' and to_ == 'bbl':
        return(x * 0.178108)
    if from_ == 'bbl' and to_ == 'ft3':
        return(x * 5.61458)
    if from_ == 'gal' and to_ == 'bbl':
        return(x * 0.02)
    if from_ == 'gal' and to_ == 'ft3':
        return(x * 0.133681)

    # permeability
    if from_ == 'm2' and to_ == 'md':
        return(x * 9.869233E+13)
    if from_ == 'ft2' and to_ == 'md':
        return((x / 10.764) * 9.869233E+13)

def nomenclatures(parameter):
    """
    Dictionary for inputs and outputs, giving the descriptions, and the oilfield units

    Input:
    parameter = the input-output parameter, string

    Output:
    description = the description of the parameter
    unit = the oilfield unit of the parameter
    """
    description = {"Bg": "gas formation volume factor",
                   "Bo": "oil formation volume factor",
                   "Bw": "water formation volume factor"}
    unit = {"Bg": "RB/scf",
            "Bo": "RB/STB",
            "Bw": "RB/STB"}

    description = description[parameter]
    unit = unit[parameter]

    return (description, unit)

# testing
print(nomenclatures("Bw"))
