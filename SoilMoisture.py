from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr
from datetime import datetime, timedelta

#Fragen:
#ToDo: Haben wir die Formel richtig verstanden, dass für SM immer die Daten verwendet werden und nicht die prediction von vorigen Tag?
#ToDo: Welchen SM sensor sollen wir verwenden?
#ToDo: Wenn Werte fehlen, z.B. eine Woche, müssten wir doch wieder bei einem gemessenen Startwert beginnen?


# Felix
path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20140517_20240517_11180_l6Xf_20240517.zip"
# Bettina
#path = r"C:\Users\betti\OneDrive\STUDIUM\SS24\Python für Geowissenschaften\SoftwareProject\Data_separate_files_header_20140517_20240517_11181_y6B1_20240517.zip"
# Theresa


# read in the data:
ismn_data = ISMN_Interface(path, parallel=False)

station_nam = "MccrackenMesa"

# function to select sensors and filter data
def station_filtered(station_nam):

    station = ismn_data['SCAN'][station_nam]

    sensor_sm = ismn_data['SCAN'][station_nam][1] #sm sensors are usually second to last
    sensor_pc = ismn_data['SCAN'][station_nam][0] #precipitation sensor is usually first

    #get lon lat from station:
    lon, lat = sensor_sm.metadata.to_dict()['longitude'][0][0], sensor_sm.metadata.to_dict()['latitude'][0][0]
    #if number not as the first entry: use filter np.array(df_2["lon"].values[0][0])[np.array(df_2["lon"].values[0][0]) != None][0] .values[0][0][0]

    sm_filter = sensor_sm.data.soil_moisture[sensor_sm.data.soil_moisture >= 0]
    pc_filter = sensor_pc.data.precipitation[sensor_pc.data.precipitation >= 0]

    sm_filter = sm_filter.resample("D").sum()
    pc_filter = pc_filter.resample("D").sum()

    return sm_filter, pc_filter, lon, lat


def align_timestamps(sm, pc):

    start_date = sm.index.min()  # start at first sm entry because we need pc from the next day
    end_date = pc.index.max()

    common_index_sm = pd.date_range(start_date, end_date - pd.Timedelta(days=1))
    common_index_pc = pd.date_range(start_date + pd.Timedelta(days=1), end_date)

    sm_synced = sm.reindex(common_index_sm)
    pc_synced = pc.reindex(common_index_pc)
    sm = sm.reindex(common_index_pc) # measured sm for pc timeframe to compare to predictions later

    sm_aligned = sm_synced.values
    pc_aligned = pc_synced.values# / 1000 # mm to m

    df_aligned = pd.DataFrame({"sm": sm.values, "sm_t_minus_1": sm_aligned, "pc_t": pc_aligned}, index=common_index_pc)

    return df_aligned
    #i9=0
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


# Function to optimise loss factor and calculate predictions in one step
def sm_prediction(station_nam):
    sm, pc, lon, lat = station_filtered(station_nam)
    df_aligned = align_timestamps(sm, pc)

    def error_function(lam, sm, pc):
        sm_pred = sm * lam + pc
        return np.sqrt(np.mean((sm_pred - sm) ** 2)) #ToDo: In Paper schauen, welche error function wir nehmen sollen

    initial_guess = 0.5 # does not make a difference if 0, 0.5 or 1
    result = minimize(error_function, initial_guess, args=(df_aligned["sm_t_minus_1"], df_aligned["pc_t"]), bounds=[(0, 1)])

    lam = result.x[0]

    sm_pred = []

    for i, sm in enumerate(df_aligned["sm_t_minus_1"]):
         if i == 0:
             pred = df_aligned["sm_t_minus_1"][i] * lam + df_aligned["pc_t"][i]
             sm_pred.append(pred)
         else:
             pred = sm_pred[-1] * lam + df_aligned["pc_t"][i]
             sm_pred.append(pred)

    sm_pred = np.array(sm_pred)

    #sm_pred = df_aligned["sm_t_minus_1"] * lam + df_aligned["pc_t"]
    df_aligned['sm_pred'] = sm_pred

    corr = spearmanr(df_aligned["sm"], sm_pred)
    rmse = np.sqrt(np.mean((sm_pred - df_aligned["sm"])**2))

    df_stations = pd.DataFrame({"station": [station_nam], "lon": [lon], "lat": [lat], "lamda": [lam], "spearman": [corr], "rmse": [rmse],
                               "sm_pred": [np.array(sm_pred)], "sm": [df_aligned["sm"].values]})


    return df_aligned, df_stations


df_1, df_2 = sm_prediction(station_nam)


#plot sm and sm_predict:
fig2, ax3 = plt.subplots(figsize=(12,4))

ax3.plot(df_1.index, df_1.sm_pred, label='Soil Moisture Prediction')
ax3.plot(df_1.index, df_1.sm, label='Soil Moisture Measured')
ax3.plot(df_1.index, df_1.pc_t, label='Precipitation')

plt.legend(loc='best')
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()
i=0
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