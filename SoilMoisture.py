from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Felix
path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20230517_20240517_11180_17tA_20240517.zip"
#Bettina

#Theresa


#read in the data:
ismn_data = ISMN_Interface(path, parallel=False)
#print(ismn_data)

#Select a Station:
station = ismn_data['SCAN']['MccrackenMesa']



i=0