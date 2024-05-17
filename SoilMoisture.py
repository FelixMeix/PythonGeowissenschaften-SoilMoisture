from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#Felix
path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20140517_20240517_11180_l6Xf_20240517.zip"
#Bettina

#Theresa


#read in the data:
ismn_data = ISMN_Interface(path, parallel=False)
#print(ismn_data)

#Select a Station:
station = ismn_data['SCAN']['MccrackenMesa']

#Select Sensor:
sensor = ismn_data['SCAN']['MccrackenMesa']['Hydraprobe-Digital-Sdi-12-(2.5-Volt)_soil_moisture_0.050800_0.050800']
#print(sensor.metadata.to_pd())

#plot the original data
ax = sensor.data.soil_moisture[sensor.data.soil_moisture >= 0].plot(figsize=(12,4))
ax.set_xlabel("Time [year]")
ax.set_ylabel("Soil Moisture [$m^3 m^{-3}$]")
#plt.savefig('firstplot.png')

i=0