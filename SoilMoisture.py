from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr, pearsonr, gamma #, expon
from datetime import timedelta
#%matplotlib inline



#Fragen:

# Felix
#path = r"C:\Users\felix\OneDrive\Dokumente\TU Wien\Geowissenschaften-Python\Data_separate_files_header_20140517_20240517_11180_l6Xf_20240517.zip"
# Bettina
#path = r"C:\Users\betti\OneDrive\STUDIUM\SS24\Python für Geowissenschaften\SoftwareProject\Data_separate_files_header_20140517_20240517_11181_y6B1_20240517.zip"
# Theresa
path = r"/Users/theresa/Documents/UIW/Master/Python-Programmierung für Geowissenschaften/Data_separate_files_header_20140517_20240517_11182_NAWF_20240517.zip"

# read in the data:
ismn_data = ISMN_Interface(path, parallel=False)



station_nam = "MccrackenMesa"
station_nam = "Mason#1"
#station_nam = "MtVernon"
#ismn_data['SCAN']['Mason#1']
#%%

# function to imput missing data based on Gamma distribution
def imput_missing(data):
    n_missing = np.sum(np.isnan(data))
    mask = np.isnan(data)
    
    n = np.sum([~mask])

    available = data[~mask]
    print("n_missing: ",n_missing)
    a, loc, scale = gamma.fit(available)
    gamma_sample = gamma.rvs(a, loc=loc, scale=scale, size=n_missing)
    #loc, scale = expon.fit(available)
    #sample = expon.rvs(loc=loc, scale=scale, size=n_missing)
    
    data[mask] = gamma_sample
    
    #plt.hist(data, bins=100, density=True, alpha=0.6, color='g', label='Data')

    # # Plot the PDF of the fitted gamma distribution
    # x = np.linspace(0, data.max(), 100)
    # pdf = gamma.pdf(x, a, loc=loc, scale=scale)
    # plt.plot(x, pdf, 'r-', lw=2, label='Fitted Gamma PDF')
    
    return data, n

# function to select sensors and filter data
def station_filtered(station_nam):

    station = ismn_data['SCAN'][station_nam]
    
    sens = station.sensors
    
    target_string = "precipitation"
    prec_sensor = [sensor for sensor in sens if target_string in sensor][0]
    
    target_string = "soil_moisture"
    sm_sensor = [sensor for sensor in sens if target_string in sensor][0]

    sensor_sm = ismn_data['SCAN'][station_nam][sm_sensor] #sm sensors are usually second to last
    sensor_pc = ismn_data['SCAN'][station_nam][prec_sensor] #precipitation sensor is usually first

    #get lon lat from station:
    lon, lat = sensor_sm.metadata.to_dict()['longitude'][0][0], sensor_sm.metadata.to_dict()['latitude'][0][0]
    #if number not as the first entry: use filter np.array(df_2["lon"].values[0][0])[np.array(df_2["lon"].values[0][0]) != None][0] .values[0][0][0]

    pc_insight = sensor_pc.data.precipitation
    
    sm_filter = sensor_sm.data.soil_moisture[sensor_sm.data.soil_moisture >= 0]
    pc_filter = sensor_pc.data.precipitation[sensor_pc.data.precipitation >= 0]
    
    # imput missing hourly values
    start_date = sm_filter.index.min()
    end_date = pc_filter.index.max()
    common_index = pd.date_range(start_date, end_date, freq="h")
    
    sm_filter, n_sm = imput_missing(sm_filter.reindex(common_index))
    pc_filter, n_pc = imput_missing(pc_filter.reindex(common_index))
    
    #sm_filter = sm_filter.resample("D").mean()
    #pc_filter = pc_filter.resample("D").sum()

    return sm_filter, pc_filter, lon, lat, n_sm, n_pc

   
#sm, pc, lon, lat = station_filtered(station_nam)

def align_timestamps(sm, pc):

    start_date = sm.index.min()  # start at first sm entry because we need pc from the next day
    end_date = pc.index.max()

    #common_index_sm = pd.date_range(start_date, end_date - pd.Timedelta(days=1))
    #common_index_pc = pd.date_range(start_date + pd.Timedelta(days=1), end_date)

    common_index_sm = pd.date_range(start_date, end_date - pd.Timedelta(hours=1), freq='H')
    common_index_pc = pd.date_range(start_date + pd.Timedelta(hours=1), end_date, freq='H')

    sm_synced = sm.reindex(common_index_sm)
    pc_synced = pc.reindex(common_index_pc)
    sm = sm.reindex(common_index_pc) # measured sm for pc timeframe to compare to predictions later

    sm_aligned = sm_synced.values
    pc_aligned = pc_synced.values# / 1000 # mm to m
    
    #sm_aligned = imput_missing(sm_aligned)
    #pc_aligned = imput_missing(pc_aligned)

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


# function to rescale predicted soil moisture to m^3/m^3 * 100
def rescale_sm(sm, sm_pred):
    #: x_scaled = (x – mean(x))/std(x) * std(y) + mean(y)
    x = sm_pred
    y = sm
    x_scaled = (x - np.mean(x))/np.std(x) * np.std(y) + np.mean(y)
    return x_scaled

#function to optimize lambda by maximizing the pearson correlation
#also results in lam = 1 for all stations
def optimize_lam(sm, pc):
    max_corr = 0
    for lam in np.arange(0, 1.01, 0.01):
        sm_pred = sm * lam + pc
        corr = pearsonr(sm_pred, sm)[0]
        if corr > max_corr:
            max_corr = corr
            optimized_lam = lam
    return optimized_lam

# Function to optimise loss factor and calculate predictions in one step
def sm_prediction(station_nam):
    sm, pc, lon, lat, n_sm, n_pc = station_filtered(station_nam)
    df_aligned = align_timestamps(sm, pc)

    def error_function(lam, sm, pc):
        pc_rescaled = rescale_sm(sm, pc) #rescale precipitation to  m^3/(m^3*100) to avoid unit problem
        sm_pred = sm * lam + pc_rescaled
        return np.sqrt(np.mean((sm_pred - sm) ** 2)) # does not work bc unit-dependent
        #return (1 - (spearmanr(sm_pred, sm)[0]))

    initial_guess = 0.5 # most lamda stay at 0.5 if initial guess
    result = minimize(error_function, initial_guess, args=(df_aligned["sm_t_minus_1"], df_aligned["pc_t"]), bounds=[(0, 1)])
    lam = result.x[0]
    
    # test with other methods, results always in lam=initial guess or lam=1 (or lam >> 1 when bounds don't work for method)
    #lam_l = []
    #methods = ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP', 'trust-constr']
    #for meth in methods:
     #   result = minimize(error_function, initial_guess, args=(df_aligned["sm_t_minus_1"], df_aligned["pc_t"]), bounds=[(0, 1)], method=meth)
      #  lam_l.append(result.x[0])
    #lam = lam_l[0]
    
    #lam = optimize_lam(df_aligned["sm_t_minus_1"], df_aligned["pc_t"])

    sm_pred = []

    for i, sm in enumerate(df_aligned["sm_t_minus_1"]):
         if i == 0:
             #pred = df_aligned["sm_t_minus_1"][i] * lam + df_aligned["pc_t"][i]
             pred = 0 # mit 0 anfangen als Startwert
             sm_pred.append(pred)
         else:
             pred = sm_pred[-1] * lam + df_aligned["pc_t"].iloc[i]
             sm_pred.append(pred)

    sm_pred = np.array(sm_pred)

    #sm_pred = df_aligned["sm_t_minus_1"] * lam + df_aligned["pc_t"]
    df_aligned['sm_pred'] = sm_pred

    corr = pearsonr(df_aligned["sm"], sm_pred)
    #rmse = np.sqrt(np.mean((sm_pred - df_aligned["sm"])**2)) #ToDo rescaling weil unterschiedliche Einheiten

    rescaled_sm = rescale_sm(df_aligned["sm"].values, df_aligned["sm_pred"].values)
    
    df_aligned["sm_pred_rescaled"] = rescaled_sm
    
    df_stations = pd.DataFrame({"station": [station_nam], "lon": [lon], "lat": [lat], "lamda": [lam], "pearson": [corr],"n_sm" : [n_sm], "n_pc" : [n_pc], #"rmse": [rmse],
                               "sm_pred_rescaled": [rescaled_sm], "sm_pred": [df_aligned["sm_pred"].values], "sm": [df_aligned["sm"].values]})


    return df_aligned, df_stations

#%%


df_1, df_2 = sm_prediction(station_nam)

#%%
#plot sm and sm_predict:
fig2, ax3 = plt.subplots(figsize=(12,4))

#ax3.plot(df_1.index, df_1.sm_pred, label='Soil Moisture Prediction [m³/m³ * 100]')
ax3.plot(df_1.index, df_1.pc_t, label='Precipitation [mm]', c="blue")
ax3.set_ylabel("mm")

ax3_2 = ax3.twinx()
ax3_2.plot(df_1.index, df_1.sm, label='Soil Moisture Measured [m³/m³ * 100]', c="green")
ax3_2.set_ylabel("m³/m³ * 100")
ax3_2.plot(df_1.index, df_1.sm_pred_rescaled, label='Soil Moisture Prediction [m³/m³ * 100]', c="red")
#ax3.plot(df_1.index, df_1.sm, label='Soil Moisture Measured [mm]', c="green")

#ax3_2.set_ylabel("mm")

ax3.set_title("Station: " + station_nam)

lines_1, labels_1 = ax3.get_legend_handles_labels()
lines_2, labels_2 = ax3_2.get_legend_handles_labels()
ax3.legend(lines_1 + lines_2, labels_1 + labels_2, bbox_to_anchor=(1.08, 0.5), loc="center left")

#plt.legend(loc='best')
plt.grid(alpha=0.4)
plt.tight_layout()
#plt.show()

#%%
#Plot for one year:

fig4, ax4 = plt.subplots(figsize=(12,4))

df_1_isel = df_1.loc['2017-01-01 00:00:00':'2017-12-31 23:00:00']

#ax3.plot(df_1.index, df_1.sm_pred, label='Soil Moisture Prediction [m³/m³ * 100]')
ax4.plot(df_1_isel.index, df_1_isel.pc_t, label='Precipitation [mm]', c="blue")
ax4.set_ylabel("mm")

ax4_2 = ax4.twinx()
ax4_2.plot(df_1_isel.index, df_1_isel.sm, label='Soil Moisture Measured [m³/m³ * 100]', c="green")
ax4_2.set_ylabel("m³/m³ * 100")
ax4_2.plot(df_1_isel.index, df_1_isel.sm_pred_rescaled, label='Soil Moisture Prediction [m³/m³ * 100]', c="red")
#ax3.plot(df_1.index, df_1.sm, label='Soil Moisture Measured [mm]', c="green")

#ax3_2.set_ylabel("mm")

ax4.set_title("Station: " + station_nam)

lines_1, labels_1 = ax4.get_legend_handles_labels()
lines_2, labels_2 = ax4_2.get_legend_handles_labels()
ax4.legend(lines_1 + lines_2, labels_1 + labels_2, bbox_to_anchor=(1.08, 0.5), loc="center left")

#plt.legend(loc='best')
plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()

#%%
# loop over all stations
stations = ['AAMU-jtg', 'Abrams', 'AdamsRanch#1', 'Alcalde', 'AlkaliMesa', 'AllenFarms', 
            'Ames', 'AshValley', 'BeasleyLake', 'Beaumont', 'BlueCreek', 'BodieHills', 
            'BraggFarm', 'BroadAcres', 'Buckhorn', 'BuffaloJump', 'BusbyFarm', 'Bushland#1', 
            'CMRBLTAR-MO', 'CacheJunction', 'CarverFarm', 'CaveValley', 'CentraliaLake', 
            'Charkiln', 'ChickenRidge', 'Circleville', 'CochoraRanch', 'ConradAgRc', 
            'CookFarmFieldD', 'Cper', 'CrescentLake#1', 'Crossroads', 'Cullman-NAHRC', 
            'DeathValleyJCT', 'DeeRiverRanch', 'DeepSprings', 'DesertCenter', 'Dexter', 
            'DoeRidge', 'Dugway', 'EagleLake', 'Eastland', 'EastviewFarm', 'ElsberryPMC', 
            'Enterprise', 'Ephraim', 'ErosDataCenter', 'Essex', 'EvergladesArs', 'FordDryLake', 
            'FortAssiniboine#1', 'FortReno#1', 'FrenchGulch', 'Geneva#1', 'GlacialRidge', 
            'GoodwinCreekPasture', 'GoodwinCreekTimber', 'Goshute', 'Grantsville', 'GreenRiver', 
            'GrouseCreek', 'HalsCanyon', 'HarmsWay', 'HartselleUSDA', 'Hodges', 'Holden', 
            'HubbardBrook', 'Hytop', 'IsbellFarms', 'JohnsonFarm', 'Jordan', 'JordanValleyCwma', 
            'JornadaExpRange', 'JournaganRanch', 'Kingsville', 'KnoxCity', 'KoptisFarms', 
            'Ku-nesa', 'KyleCanyon', 'LINDSAY', 'Levelland', 'Lind#1', 'LittleRedFox', 
            'LittleRiver', 'Livingston-UWA', 'LosLunasPmc', 'LovellSummit', 'LovelockNnr', 
            'LyeBrook', 'LynhartRanch', 'MahantangoCk', 'MammothCave', 'Mandan#1', 
            'Manderfield', 'MarbleCreek', 'MarkTwainHS', 'MascomaRiver', 'Mason#1', 'Mayday', 
            'McalisterFarm', 'MccrackenMesa', 'Milford', 'Moccasin', 'MollyCaren#1', 'MonoclineRidge', 
            'Morgan', 'MorrisFarms', 'MountMansfield', 'MountainHome', 'MtVernon', 'NPiedmontArec', 
            'Nephi', 'NorthIssaquena', 'Nunn#1', 'Onward', 'OrchardRangeSite', 'Panguitch', 'ParkValley', 
            'PeeDee', 'PerdidoRivFarms', 'Perthshire', 'Phillipsburg', 'PineNut', 'PorterCanyon', 
            'PowderMill', 'PowellGardens', 'PrairieView#1', 'Price', 'Princeton#1', 
            'ReynoldsHomestead', 'Riesel', 'RiverRoadFarms', 'RockSpringsPa', 'RogersFarm#1', 
            'SHELDON', 'SanAngelo', 'SandHollow', 'SandyRidge', 'Schell-Osage', 'Scott', 
            'SellersLake#1', 'Selma', 'Sevilleta', 'ShadowMtns', 'ShagbarkHills', 'ShawNatureReserve', 
            'SheepSpringsWX', 'Shenandoah', 'Sidney', 'SilverCity', 'Spickard', 'SplitMountain', 
            'Spooky', 'StanleyFarm', 'Starkville', 'Stephenville', 'Stubblefield', 'SudduthFarms', 
            'SunleafNursery', 'TABLEMOUNTAIN', 'TNCFortBayou', 'Tidewater#1', 'TidewaterArec', 
            'Torrington#1', 'TroughSprings', 'TuleValley', 'Tunica', 'Tuskegee', 'TwinPinesConservationArea', 
            'UAPBCampus-PB', 'UAPBDewitt', 'UAPBEarle', 'UAPBLonokeFarm', 'UAPBMarianna', 
            'UAPBPointRemove', 'Uvalde', 'UwPlatteville', 'Vance', 'Vermillion', 'Vernon', 'Violett', 
            'WTARS', 'Wabeno#1', 'Wakulla#1', 'WalnutGulch#1', 'Watkinsville#1', 'Wedowee', 'Weslaco', 
            'WestSummit', 'WillowWells', 'YoumansFarm']

result = pd.DataFrame(columns=["station", "lon", "lat", "lamda", "pearson","sm_pred", "sm"]) #"rmse"
for station in stations:
    if station == 'Sidney':             # to avoid FitErrors at these stations
        print('Station Sidney passed')
    elif station == 'Violett':
        print('Station Violett passed')
    else:
        name = station
        station = ismn_data['SCAN'][name]
        sens = station.sensors
        if len([sensor for sensor in sens if "precipitation" in sensor])==0:
            print(name)
            pass
        else:
            df_1, df_2 = sm_prediction(name)
            result = pd.concat([result, df_2], axis=0, ignore_index=True)
#%%

print(result['lamda'][30:60])            
print(min(result['lamda']))
print(max(result['lamda']))
np.median(result['lamda'])


#%%
###
result["corr_coef"] = [result["pearson"][row][0] for row in range(len(result["lon"]))] 

lat_north = 50
lat_south = 23
lon_west = -120
lon_east = -75

import cartopy.crs as ccrs
import cartopy.feature as cfeature

#ax = plt.axes(projection = ccrs.LambertConformal())
fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.LambertConformal()})


ax.set_extent([lon_west, lon_east, lat_south, lat_north])
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

sc = ax.scatter(result["lon"], result["lat"], transform=ccrs.PlateCarree(), c=result["corr_coef"], cmap="plasma") #, label="SCAN", edgecolors="k"
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.5)
cbar.set_label('Correlation Coefficient')

ax.set_title("Pearson: Measured Precipitation and Predicted Soil Moisture")
gl = ax.gridlines(draw_labels = True, x_inline = False, y_inline = False, alpha = 0.5, linestyle = "--")
gl.top_labels = False
gl.right_labels = False

#ax.legend()




#i=0
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
