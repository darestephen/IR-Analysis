import csv
from re import search
from tkinter import filedialog, Tk
import os
import platform
import io
from pandas.errors import EmptyDataError
import jcamp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
from bs4 import BeautifulSoup
import sqlite3
import steve
from pathlib import Path
import Fourier
import Get_jdx
import Chem_Cal
import Multiple_Compounds

def MainMenu():
    while True:
        try:
            select = int(input(
                "\nMain Menu:\n"
                "1: Download JCAMP IR spectrum (NIST)\n"
                "2: Perform Fourier Transform on IR data\n"
                "3: Chemistry calculations\n"
                #"4: Identify an Unknown Compound\n"
                "4: Identify Multiple Compounds\n"
                "0: Exit\n"
                "Enter your choice: "
            ))

            if select == 1:
                Get_jdx.process_multiple_compounds()

            elif select == 2:
                Fourier.process_jdx_file()

            elif select == 3:
                Chem_Cal.Chem1()
            #elif select == 4:
             #   steve.Unknown_compound_identifier()

            elif select == 4:
                Multiple_Compounds.Multiple() 
            elif select == 0:
                print("Logging out... 👋")
                break

            else:
                print("Invalid selection. Try again.")
        except ValueError:
            print("❌ Please enter a valid number.")
if __name__ == "__main__":
    MainMenu()


