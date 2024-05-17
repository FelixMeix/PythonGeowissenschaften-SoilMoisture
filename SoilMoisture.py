from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr

# Felix
path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20140517_20240517_11180_l6Xf_20240517.zip"
# Bettina
#path = r"C:\Users\betti\OneDrive\STUDIUM\SS24\Python für Geowissenschaften\SoftwareProject\Data_separate_files_header_20140517_20240517_11181_y6B1_20240517.zip"
# Theresa


# read in the data:
ismn_data = ISMN_Interface(path, parallel=False)

# Select a Station and sensors:
station_nam = "MccrackenMesa"
station = ismn_data['SCAN'][station_nam]

sensor_sm = ismn_data['SCAN'][station_nam]['Hydraprobe-Digital-Sdi-12-(2.5-Volt)_soil_moisture_0.101600_0.101600']
sensor_pc = ismn_data['SCAN'][station_nam]['Pulse-Count_precipitation_0.000000_0.000000']

sm_filter = sensor_sm.data.soil_moisture[sensor_sm.data.soil_moisture >= 0]
pc_filter = sensor_pc.data.precipitation[sensor_pc.data.precipitation >= 0]

#plot the original timeseries:
# fig, ax1 = plt.subplots(figsize=(12,4))
#
# pc_filter.plot(ax=ax1, label='Precipitation', linewidth=0.8)
#
# ax2 = ax1.twinx()
# sm_filter.plot(ax=ax2, alpha=0.65, label='Soil Moisture', color='darkgreen', linewidth=0.5)
#
#
# ax1.set_xlabel("Time [year]")
# ax1.set_ylabel("Precipitation [mm]")
# ax2.set_ylabel("Soil Moisture [$m^3 m^{-3}$]")
# ax1.set_title("Soil Moisture and Precipitation at station " + station_nam)
#
# lines_1, labels_1 = ax1.get_legend_handles_labels()
# lines_2, labels_2 = ax2.get_legend_handles_labels()
# ax1.legend(lines_1 + lines_2, labels_1 + labels_2, bbox_to_anchor=(1.08, 0.5), loc="center left")
#
# plt.grid(alpha=0.4)
# plt.tight_layout()
# plt.show()

#Function for SM-Model
sm = sm_filter.values
pc = pc_filter.values #ToDo: pc tag um eins verschieben

start_date = max(sm_filter.index.min(), pc_filter.index.min())
end_date = min(sm_filter.index.max(), pc_filter.index.max())

# Serien auf den gemeinsamen Zeitbereich beschränken
common_index = pd.date_range(start_date, end_date)
sm_filter_synced = sm_filter.reindex(common_index)
pc_filter_synced = pc_filter.reindex(common_index)

# Kombinieren der synchronisierten Serien in einen DataFrame
combined_df = pd.DataFrame({'sm_filter': sm_filter_synced, 'pc_filter': pc_filter_synced})

# Löschen der Zeilen mit fehlenden Werten
combined_df_cleaned = combined_df.dropna()

# Getrennte Serien extrahieren
sm_filter_cleaned = combined_df_cleaned['sm_filter']
pc_filter_cleaned = combined_df_cleaned['pc_filter']
def sm_prediction():
    sm = sm_filter_cleaned.values
    pc = pc_filter_cleaned.values

    def error_function(lam, sm, pc):
        sm_pred = sm * lam + pc
        return np.sqrt(np.mean((sm_pred - sm) ** 2))

    initial_guess = 0.5
    result = minimize(error_function, initial_guess, args=(sm, pc))#, bounds=(0, 1))

    lam = result.x[0]

    sm_pred = sm * lam + pc

    corr = spearmanr(sm, sm_pred) #ToDo: corr ev anpassen --> anschauen

    return sm_pred, corr, lam

#sm = sm_filter_cleaned.values
#pc = pc_filter_cleaned.values

predicted_sm, corr, lam = sm_prediction()


i = 0