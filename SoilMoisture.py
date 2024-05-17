%pip install ismn
import ismn
from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Felix
#path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20230517_20240517_11180_17tA_20240517.zip"
#Bettina
path = r"C:\Users\betti\OneDrive\STUDIUM\SS24\Python für Geowissenschaften\SoftwareProject\Data_separate_files_header_20140517_20240517_11181_y6B1_20240517.zip"
#Theresa


#read in the data:
ismn_data = ISMN_Interface(path, parallel=False)
#print(ismn_data)

#Select a Station:
station_nam = "MccrackenMesa"  
    
station = ismn_data['SCAN'][station_nam]
station.metadata


sensor = ismn_data['SCAN'][station_nam]['Hydraprobe-Digital-Sdi-12-(2.5-Volt)_soil_moisture_0.101600_0.101600']
print(sensor.metadata.to_pd())
ax = sensor.data.soil_moisture[sensor.data.soil_moisture >= 0].plot(figsize=(12,4))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Soil Moisture [$m^3 m^{-3}$]")
ax.set_title("Soil Moisture at station " + station_nam)

sensor = ismn_data['SCAN'][station_nam]['Pulse-Count_precipitation_0.000000_0.000000']
print(sensor.metadata.to_pd())
ax = sensor.data.precipitation[sensor.data.precipitation >= 0].plot(figsize=(12,4))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Precipitation")
ax.set_title("Precipitation at station " + station_nam)


i=0