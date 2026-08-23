def Chem1():
    syt = float(
        input("What system is it:\n1=Adiabatic\n2=Isothermal\n3=Isochoric\n4=Isobaric\n5=Free Space expansion\n"))
    while syt > 5 or syt < 1:
        syt = int(
            input("What system is it:\n1=Adiabatic\n2=Isothermal\n3=Isochoric\n4=Isobaric\n5=Free Space expansion\n"))
    if syt == 1.0:
        que = str(input("What do you want to Calculate i.e ft for finaltemperature: "))
        dar = que.replace(" ", "")
        Cal = dar.lower()
        par = []
        if Cal == "ft":
            try:
                num = int(input("How many parameters are given: "))
                while num >= 5:
                    num = int(input("How many parameters are given: "))
            except:
                print("Invalid input! Enter a valid number less than 5")
            else:
                for i in range(1, num + 1):
                    parameter = str(input(f"Enter no {i} parameter :"))
                    no_space = parameter.replace(" ", "")
                    para = no_space.lower()
                    par.append(para)
                if "it" in par and "iv" in par and "fv" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    initialtemperature = float(input("Enter the Initial Temperature(in Kelvin) i.e:it: "))
                    initialvolume = float(input("Enter the Initial Volume(in centimeter cube):i.e iv: "))
                    finalvolume = float(input("Enter the final volume(in centimeter cube) :i.e fv: "))
                    z = (initialvolume / finalvolume)
                    a = float(Gamma - 1)
                    b = z ** a
                    finalltemperature = initialtemperature * b
                    roundabout = round(finalltemperature, 2)
                    print(f"Final Tmperature {roundabout} K")
                elif "ip" in par and "fp" in par and "it" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    initialtemperature = float(input("Enter the Initial Temperature(in Kelvin):i.e it: "))
                    initialpressure = float(input("Enter the Initial Pressure(in bar):i.e ip: "))
                    finalpressure = float(input("Enter the final Pressure(in bar):i.e fp: "))
                    z = (finalpressure / initialpressure)
                    a = float(Gamma - 1)
                    p = z ** a
                    x = (initialtemperature) ** Gamma
                    y = (p) * (x)
                    finaltemperature = y ** (1 / Gamma)
                    roundabout = round(finaltemperature, 2)
                    print(f"Final Tmperature {roundabout} K")
        elif Cal == "it":
            try:
                num = int(input("How many parameters are given: "))
                while num >= 5:
                    num = int(input("How many parameters are given: "))
            except:
                print("Invalid input! Enter a valid number less than 5")
            else:
                for i in range(1, num + 1):
                    parameter = str(input(f"Enter no {i} parameter :"))
                    no_space = parameter.replace(" ", "")
                    para = no_space.lower()
                    par.append(para)
                if "ft" in par and "iv" in par and "fv" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    finaltemperature = float(input("Enter the final Temperature(in Kelvin): "))
                    initialvolume = float(input("Enter the Initial Volume(in centimeter cube): "))
                    finalvolume = float(input("Enter the final volume(in centimeter cube): "))
                    z = (initialvolume / finalvolume)
                    a = float(Gamma - 1)
                    b = z ** a
                    initialtemperature = finaltemperature / b
                    roundabout = round(initialtemperature, 2)
                    print(f"Final Temperature {roundabout} K")
                elif "ip" in par and "fp" in par and "ft" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    finaltemperature = float(input("Enter the Initial Temperature(in Kelvin): "))
                    initialpressure = float(input("Enter the Initial Pressure(in bar): "))
                    finalpressure = float(input("Enter the final Pressure(in bar): "))
                    z = (initialpressure / finalpressure)
                    a = float(Gamma - 1)
                    p = z ** a
                    x = (finaltemperature) ** Gamma
                    y = (p) * (x)
                    initialtemperature = y ** (1 / Gamma)
                    roundabout = round(initialtemperature, 2)
                    print(f"Initial Tmperature {roundabout} K")
        elif Cal == "iv":
            try:
                num = int(input("How many parameters are given: "))
                while num >= 5:
                    num = int(input("How many parameters are given: "))
            except:
                print("Invalid input! Enter a valid number less than 5")
            else:
                for i in range(1, num + 1):
                    parameter = str(input(f"Enter no {i} parameter :"))
                    no_space = parameter.replace(" ", "")
                    para = no_space.lower()
                    par.append(para)
                if "ft" in par and "it" in par and "fv" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    finaltemperature = float(input("Enter the final Temperature(in Kelvin): "))
                    initialtemperature = float(input("Enter the Initial temperature(in centimeter cube): "))
                    finalvolume = float(input("Enter the final volume(in centimeter cube): "))
                    z = (finaltemperature / initialtemperature)
                    a = float(Gamma - 1)
                    b = (finalvolume) ** a
                    c = z * b
                    initialvolume = c ** (1 / a)
                    roundabout = round(initialvolume, 2)
                    print(f"Initial Volume {roundabout} cm^3")
                elif "fv" in par and "ip" in par and "fp" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    initialpressure = float(input("Enter the Initial Pressure(in N/m^2): "))
                    finalvolume = float(input("Enter the Final Volume(in meter cube): "))
                    finalpressure = float(input("Enter the final pressure(in N/m^2): "))
                    z = (finalpressure / initialpressure)
                    b = (finalvolume) ** (Gamma)
                    c = z * b
                    initialvolume = c ** (1 / Gamma)
                    roundabout = round(initialvolume, 2)
                    print(f"Initial Volume {roundabout} m^3")
        elif Cal == "fv":
            try:
                num = int(input("How many parameters are given: "))
                while num >= 5:
                    num = int(input("Should be less than 5 parameters: "))
            except:
                print("Invalid input! Enter a valid number less than 5")
            else:
                for i in range(1, num + 1):
                    parameter = str(input(f"Enter no {i} parameter :"))
                    no_space = parameter.replace(" ", "")
                    para = no_space.lower()
                    par.append(para)
                if "ft" in par and "iv" in par and "it" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    finaltemperature = float(input("Enter the final Temperature(in Kelvin): "))
                    initialtemperature = float(input("Enter the Initial Temperature(in kelvin): "))
                    initialvolume = float(input("Enter the initial volume(in centimeter cube): "))
                    z = (initialtemperature / finaltemperature)
                    a = float(Gamma - 1)
                    b = (initialvolume) ** a
                    c = z * b
                    finalvolume = c ** (1 / a)
                    roundabout = round(finalvolume, 2)
                    print(f"Initial Volume {roundabout} cm^3")
                elif "iv" in par and "ip" in par and "fp" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    initialpressure = float(input("Enter the Initial Pressure(in N/m^2): "))
                    initialvolume = float(input("Enter the Initial Volume(in meter cube): "))
                    finalpressure = float(input("Enter the final pressure(in N/m^2): "))
                    z = (initialpressure / finalpressure)
                    b = (initialvolume) ** (Gamma)
                    c = z * b
                    finalvolume = c ** (1 / Gamma)
                    roundabout = round(finalvolume, 2)
                    print(f"Final Volume {roundabout} m^3")
        elif Cal == "fp":
            try:
                num = int(input("How many parameters are given: "))
                while num >= 5:
                    num = int(input("Should be less than 5 parameters: "))
            except:
                print("Invalid input! Enter a valid number less than 5")
            else:
                for i in range(1, num + 1):
                    parameter = str(input(f"Enter no {i} parameter :"))
                    no_space = parameter.replace(" ", "")
                    para = no_space.lower()
                    par.append(para)
                if "ft" in par and "ip" in par and "it" in par:
                    print(
                        "Always remember for \nmonoatomic gas Cv=(3/2)R, Cp=Cv+R\nFor diatomic gas Cv=(5/2)R, Cp = Cv+R\nfor polyatomic gas Cv=3R,Cp=Cv+R")
                    Gamma = float(input("Enter the value of γ =Cp/Cv: "))
                    finaltemperature = float(input("Enter the final Temperature(in Kelvin): "))
                    initialtemperature = float(input("Enter the Initial Temperature(in kelvin): "))
                    initialpressure = float(input("Enter the Initial Pressure(in bar): "))
                    z = (finaltemperature / initialtemperature)
                    a = float(Gamma - 1)
                    b = (initialpressure) ** a
                    c = z ** Gamma
                    x = c * b
                    finalpressure = x ** (1 / a)
                    roundabout = round(finalpressure, 2)
                    print(f"Final Pressure:  {roundabout} bar")


if __name__ == "__main__":
    Chem1()