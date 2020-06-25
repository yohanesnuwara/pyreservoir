def nomenclatures(parameter):
    nomen = {"Bg": "gas formation volume factor", "Bo": "oil FVF", "Bw": "water FVF"}
    units = {"Bg": "RB/scf", "Bo": "RB/STB", "Bw": "RB/STB"}

    whatis = nomen[parameter]
    unit = units[parameter]

    return(whatis, unit)

# testing
print(nomenclatures("Bg"))
