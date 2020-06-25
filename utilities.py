def nomenclatures(parameter):
    """
    Dictionary for inputs and outputs, giving the descriptions, and the oilfield units
    
    Input:
    parameter = the input-output parameter, string
    
    Output:
    description = the description of the parameter
    unit = the oilfield unit of the parameter
    """
    description = {"Bg": "gas formation volume factor", "Bo": "oil FVF", "Bw": "water FVF"}
    unit = {"Bg": "RB/scf", "Bo": "RB/STB", "Bw": "RB/STB"}

    description = description[parameter]
    unit = unit[parameter]

    return(description, unit)

# testing
print(nomenclatures("Bg"))
