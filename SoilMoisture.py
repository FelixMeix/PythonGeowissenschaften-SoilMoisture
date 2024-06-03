from ismn.interface import ISMN_Interface
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import minimize
from scipy.stats import spearmanr, pearsonr, gamma #, expon
from datetime import timedelta
#%matplotlib inline
from IPython import get_ipython
get_ipython().run_line_magic('matplotlib', 'inline')



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

# station with the least missing values:
<<<<<<< HEAD
#station_nam = "Hytop" # sm n_missing:  1; pc n_missing:  0 #only 10 months
=======
#station_nam = "Hytop" # sm n_missing:  1; pc n_missing:  0; only up to 2015
>>>>>>> refs/remotes/origin/main


station_nam = "TidewaterArec" # sm n_m missing: 24; pc n_missing: 25; up to 2019/06/20

#station_nam = "MtVernon" # sm n_m missing: 34; pc n_missing: 40; up to 2024 (second choice)
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
    
    data_imputed = data.copy()
    data_imputed[mask] = gamma_sample
    
    #plt.hist(data, bins=100, density=True, alpha=0.6, color='g', label='Data')

    # # Plot the PDF of the fitted gamma distribution
    #x = np.linspace(0, data.max(), 100)
    #pdf = gamma.pdf(x, a, loc=loc, scale=scale)
    #plt.plot(x, pdf, 'r-', lw=2, label='Fitted Gamma PDF')
    
    return data_imputed, n

# function to select sensors and filter data
def station_filtered(station_nam):

    station = ismn_data['SCAN'][station_nam]
    
    sens = station.sensors
    
    target_string = "precipitation"
    prec_sensor = [sensor for sensor in sens if target_string in sensor][-1] #changed from 0 to -1 to have pulse count pc sensor (newer data)
    
    target_string = "soil_moisture"
    sm_sensor = [sensor for sensor in sens if target_string in sensor][0]

    sensor_sm = ismn_data['SCAN'][station_nam][sm_sensor] #sm sensors are usually second to last
    sensor_pc = ismn_data['SCAN'][station_nam][prec_sensor] #precipitation sensor is usually first

    #get lon lat from station:
    lon, lat = sensor_sm.metadata.to_dict()['longitude'][0][0], sensor_sm.metadata.to_dict()['latitude'][0][0]
    #if number not as the first entry: use filter np.array(df_2["lon"].values[0][0])[np.array(df_2["lon"].values[0][0]) != None][0] .values[0][0][0]

    pc_insight = sensor_pc.data.precipitation
    
    sm_filter = sensor_sm.data.soil_moisture[sensor_sm.data.soil_moisture >= 0]
    pc_filter = sensor_pc.data.precipitation[sensor_pc.data.precipitation >= 0][sensor_pc.data.precipitation < 100]
    
    # imput missing hourly values
    start_date = sm_filter.index.min()
    end_date = pc_filter.index.max()
    common_index = pd.date_range(start_date, end_date, freq="h")
    
    
    pc_filter_unimputed = pc_filter.reindex(common_index)
    #print(station_nam)
    #print("sm")
    sm_filter, n_sm = imput_missing(sm_filter.reindex(common_index))
    #print("pc")
    pc_filter, n_pc = imput_missing(pc_filter.reindex(common_index))
    
    #sm_filter = sm_filter.resample("D").mean()
    #pc_filter = pc_filter.resample("D").sum()

    return sm_filter, pc_filter, pc_filter_unimputed, lon, lat, n_sm, n_pc

   
#sm, pc, lon, lat = station_filtered(station_nam)

def align_timestamps(sm, pc, pc_unimputed):

    start_date = sm.index.min()  # start at first sm entry because we need pc from the next day
    end_date = pc.index.max()

    #common_index_sm = pd.date_range(start_date, end_date - pd.Timedelta(days=1))
    #common_index_pc = pd.date_range(start_date + pd.Timedelta(days=1), end_date)

    common_index_sm = pd.date_range(start_date, end_date - pd.Timedelta(hours=1), freq='H')
    common_index_pc = pd.date_range(start_date + pd.Timedelta(hours=1), end_date, freq='H')

    sm_synced = sm.reindex(common_index_sm)
    pc_synced = pc.reindex(common_index_pc)
    sm = sm.reindex(common_index_pc) # measured sm for pc timeframe to compare to predictions later
    pc_unimputed_synced = pc_unimputed.reindex(common_index_pc)

    sm_aligned = sm_synced.values
    pc_aligned = pc_synced.values# / 1000 # mm to m
    pc_unimputed_aligned = pc_unimputed_synced.values
    
    #sm_aligned = imput_missing(sm_aligned)
    #pc_aligned = imput_missing(pc_aligned)

    df_aligned = pd.DataFrame({"sm": sm.values, "sm_t_minus_1": sm_aligned, "pc_t": pc_aligned, "pc_unimputed": pc_unimputed_aligned}, index=common_index_pc)

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
    print(station_nam)
    sm, pc, pc_unimputed, lon, lat, n_sm, n_pc = station_filtered(station_nam)
    df_aligned = align_timestamps(sm, pc, pc_unimputed)

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

    corr_pearson = pearsonr(df_aligned["sm"], sm_pred)
    corr_spearman = spearmanr(df_aligned["sm"], sm_pred)
    #rmse = np.sqrt(np.mean((sm_pred - df_aligned["sm"])**2)) #unten mit rescaled

    rescaled_sm = rescale_sm(df_aligned["sm"].values, df_aligned["sm_pred"].values)
    
    rmse = np.sqrt(np.mean((rescaled_sm - df_aligned["sm"].values)**2))
    
    df_aligned["sm_pred_rescaled"] = rescaled_sm
    
    df_stations = pd.DataFrame({"station": [station_nam], "lon": [lon], "lat": [lat], "lamda": [lam], "pearson": [corr_pearson], "spearman": [corr_spearman],"n_sm" : [n_sm], "n_pc" : [n_pc], "rmse": [rmse],
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
ax3_2.plot(df_1.index, df_1.sm*100*100, label='Soil Moisture Measured [m³/m³ * 100]', c="green")
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
ax4_2.plot(df_1_isel.index, df_1_isel.sm*100, label='Soil Moisture Measured [m³/m³ * 100]', c="green")
ax4_2.set_ylabel("m³/m³ * 100")
ax4_2.plot(df_1_isel.index, df_1_isel.sm_pred_rescaled*100, label='Soil Moisture Prediction [m³/m³ * 100]', c="red")
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

result = pd.DataFrame(columns=["station", "lon", "lat", "lamda", "pearson", "spearman","sm_pred", "sm"]) #"rmse"
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

#print(result['lamda'][:30])            
print(min(result['lamda']))
print(max(result['lamda']))
print(np.median(result['lamda']))

# find and compare two stations with different lamda

low_lamda = result.loc[result["n_sm"] > 10000]
low_lamda = low_lamda["station"].loc[low_lamda["lamda"].idxmin()]
high_lamda = result["station"].iloc[np.argmax(result["lamda"])]

subresult = result.loc[result["station"].isin([low_lamda, high_lamda])]

import cartopy.crs as ccrs
import cartopy.feature as cfeature

lat_north = 50
lat_south = 23
lon_west = -120
lon_east = -75

fig = plt.figure(figsize=(10, 6))#, subplot_kw={'projection': ccrs.LambertConformal()})

ax = plt.axes(projection=ccrs.LambertConformal())
ax.set_extent([lon_west, lon_east, lat_south, lat_north])
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

#sc = ax.scatter(subresult["lon"], subresult["lat"], transform=ccrs.PlateCarree(), color="red", ec="k", s=10) #, label="SCAN", edgecolors="k"

ax.set_title("Two stations with highest and lowest lamda")

station = low_lamda
for station in [low_lamda, high_lamda]:
    res = result.loc[result["station"].isin([station])]
    ax.plot(res['lon'], res['lat'], 'ro', transform=ccrs.PlateCarree())
    ax.annotate(
        f"Station: {station}\n loss factor: {round(float(res['lamda']),2)}\n spearman: {round(float(res['corr_coef']),2)}\n rmse: {round(float(res['rmse']),2)}\n n sm: {int(res['n_sm'])}",
        xy=(res['lon'], res['lat']),
        xycoords= ccrs.LambertConformal()._as_mpl_transform(ax),
        xytext=(10,10),
        textcoords='offset points',
        arrowprops=dict(facecolor='black', shrink=0.05),
        bbox=dict(boxstyle='round,pad=0.5', edgecolor='black', facecolor='wheat', alpha=0.5)
        )
    #plt.text(res["lon"], res["lat"], f"Station: {station}\n loss factor: {round(float(res['lamda']),2)}\n spearman: {round(float(res['corr_coef']),2)}\n rmse: {round(float(res['rmse']),2)}\n n sm: {int(res['n_sm'])}",
     #        horizontalalignment="left", transform=ccrs.PlateCarree())

gl = ax.gridlines(draw_labels = True, x_inline = False, y_inline = False, alpha = 0.5, linestyle = "--")
gl.top_labels = False
gl.right_labels = False

plt.savefig("high_low_lamda_cartopy.jpg", dpi=300)

# and plot for 1 yr
fig, ax = plt.subplots(1,2,figsize=(12,4))

for i,station_nam in enumerate([low_lamda, high_lamda]):
    df_1, df_2 = sm_prediction(station_nam)
    
    df_1 = df_1.loc['2017-01-01 00:00:00':'2017-12-31 23:00:00']
    
    ax[i].plot(df_1.index, df_1.sm_pred, label='Soil Moisture Prediction [mm]', c="green")
    ax[i].set_ylabel("mm")
    ax[i].plot(df_1.index, df_1.pc_t, label='Precipitation [mm]', c="red", linewidth=0.5, linestyle="--")
    
    ax[i].set_title("Station: " + station_nam + f"\n lamda = {round(float(df_2['lamda']),4)}")

    ax[i].legend(loc='best')
    ax[i].grid(alpha=0.4)
plt.tight_layout()
plt.show()

plt.savefig("high_low_lamda_prec_sm_mm.jpg", dpi=300)

#%%
###
result["corr_coef"] = [result["spearman"][row][0] for row in range(len(result["lon"]))] 

# cartopy plot for correlation  #added number of measurements as size

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.LambertConformal()})


ax.set_extent([lon_west, lon_east, lat_south, lat_north])
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

sc = ax.scatter(result["lon"], result["lat"], transform=ccrs.PlateCarree(), 
                c=result["corr_coef"], cmap="plasma", s=result["n_sm"]/1000,
                ec="k") #, label="SCAN", edgecolors="k"
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.5)
cbar.set_label('Correlation Coefficient')

plt.title("Spearman: Measured Soil Moisture and Predicted Soil Moisture \n scaled after number of available measurements")
#plt.suptitle("scaled after number of available measurements", y=1)
gl = ax.gridlines(draw_labels = True, x_inline = False, y_inline = False, alpha = 0.5, linestyle = "--")
gl.top_labels = False
gl.right_labels = False

plt.savefig("Spearman_cartopy.jpg", dpi=300)

# scatterplot for correlation coefficient and lamda

fig, ax = plt.subplots(figsize=(12,4))

#ax3.plot(df_1.index, df_1.sm_pred, label='Soil Moisture Prediction [m³/m³ * 100]')
ax.scatter(result["lamda"], result["corr_coef"], c="blue")
ax.set_ylabel("spearman correlation coefficient")
ax.set_xlabel("loss coefficient lamda") #r$\lamda$

ax.set_title("all stations - spearman correlation over lamda")

plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()

plt.savefig("corr_lamda_scatter.jpg", dpi=300)

# cartopy for rmse

fig, ax = plt.subplots(figsize=(10, 6), subplot_kw={'projection': ccrs.LambertConformal()})

ax.set_extent([lon_west, lon_east, lat_south, lat_north])
ax.add_feature(cfeature.COASTLINE)
ax.add_feature(cfeature.STATES)
ax.add_feature(cfeature.LAND)
ax.add_feature(cfeature.OCEAN)

sc = ax.scatter(result["lon"], result["lat"], transform=ccrs.PlateCarree(), c=result["rmse"]*100, 
                cmap="gist_rainbow_r", ec="k", vmin=0, vmax=95) #, label="SCAN", edgecolors="k"
cbar = plt.colorbar(sc, ax=ax, orientation='vertical', shrink=0.5)
cbar.set_label('Root Mean Squared Error [m³/m³ * 100]')

#ax.set_title("Pearson: Measured Precipitation and Predicted Soil Moisture")
ax.set_title("RSME: Measured Precipitation and Predicted Soil Moisture")
gl = ax.gridlines(draw_labels = True, x_inline = False, y_inline = False, alpha = 0.5, linestyle = "--")
gl.top_labels = False
gl.right_labels = False

plt.savefig("RMSE_cartopy.jpg", dpi=300)


# histogram for lamda 

fig, ax = plt.subplots(figsize=(12,4))

ax.hist(result["lamda"], ec="k", fc="lightblue", bins=20, density=False, cumulative=False)
ax.set_ylabel("absolute frequency")
ax.set_xlabel("loss coefficient lamda")

ax.set_title("all stations - lamda ")

plt.grid(alpha=0.4)
plt.tight_layout()
plt.show()

plt.savefig("Lamda_histo.jpg", dpi=300)

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



#%%
# for chosen station set percentage of pc values to nan artifically

station_nam = "TidewaterArec" # sm n_m missing: 24; pc n_missing: 25; up to 2019/06/20
df_1, df_2 = sm_prediction(station_nam)
ds_pc = df_1['pc_unimputed']

percentage = [10, 20, 30, 50, 70]
ds_del_l = []

for perc in percentage:
    ds_del = ds_pc.copy()
    factor = int((len(ds_pc)*perc)/100)
    ds_sample = ds_del.sample(factor)
    ds_del.loc[ds_sample.index] = np.nan
    print(ds_del.isnull().sum()/len(ds_del)) # to verify
    ds_del_l.append(ds_del)

# create df with all different NaN percentages
df_del = pd.concat([ds for ds in ds_del_l], axis=1)
df_del.columns =[f'{str(perc)}% NaN' for perc in percentage]


#%%

# refill nan values in 3 different ways:

# 1. gamma distribution

# function to fill nan values based on Gamma distribution
def refill_gamma(ds):
    n_missing = ds.isnull().sum()
    available = ds[ds.notnull()]
    
    #available = pd.to_numeric(available, errors='coerce') # to avoid input error
    
    a, loc, scale = gamma.fit(available)
    gamma_sample = gamma.rvs(a, loc=loc, scale=scale, size=n_missing)
    
    ds_filled = ds.copy()
    ds_filled[ds_filled.isnull()] = gamma_sample
    
    #x = np.linspace(0, np.max(ds_filled), 100)
    #pdf = gamma.pdf(x, a, loc=loc, scale=scale)
    #plt.plot(x, pdf, 'r-', lw=2, label='Fitted Gamma PDF')

    return ds_filled

ds_gamma_refilled_l = []
for column in df_del.columns:
    ds_gamma_refilled_l.append(refill_gamma(df_del[column]))

# df with all gamma refilled series
df_gamma_refilled = pd.concat([ds for ds in ds_gamma_refilled_l], axis=1)
df_gamma_refilled.plot()

#%%
# 2. machine learning

from sklearn import metrics
from sklearn.linear_model import LinearRegression
from sklearn.svm import SVR
from sklearn.tree import DecisionTreeRegressor

ds_sm = df_1['sm']

ds_y = df_del['10% NaN'].copy()
NAN_indices = ds_y[ds_y.isnull()].index
ds_x = ds_sm.copy()


# Training Data
y_trn = ds_y.drop(NAN_indices , axis=0).copy()
X_trn = ds_x.drop(NAN_indices , axis=0).copy()

# Testing Data
y_tst = ds_pc.loc[NAN_indices]      # y_tst are the original data
X_tst = ds_x.loc[NAN_indices]

print('Training set shape : ', X_trn.shape)
print('Testing set shape : ', X_tst.shape)

# Convert Training and Testing dataframes to arrays to make them ready for modeling
X_trn = X_trn.to_numpy().reshape(-1,1)
y_trn = y_trn.to_numpy().reshape(-1,1)
X_tst = X_tst.to_numpy().reshape(-1,1)
y_tst = y_tst.to_numpy().reshape(-1,1)

# function to check accuracy
def acc_chk (y_tst , y_hat):
    # Mean absolute error
    MAE = metrics.mean_absolute_error(y_tst ,y_hat)
    
    # Mean Suared Error (MSE)
    MSE = metrics.mean_squared_error(y_tst ,y_hat)
    
    # R2-score
    R2_sc = metrics.r2_score(y_tst ,y_hat)
    
    print("Mean absolute error: %.2f" % MAE)
    print("Mean Squared Error : %.2f" % MSE)
    print("R2-score           : %.2f" % R2_sc )
    
    return MAE , MSE , R2_sc

# Linear Regression

LR = LinearRegression(fit_intercept=True,copy_X=True,n_jobs=-1)
LR.fit(X_trn ,y_trn)
yhat_NANLR = LR.predict(X_tst)

print ('Accuracy scores for filling NAN with Linear Regression model are:')
NANLR_MAE , NANLR_MSE , NANLR_R2_sc = acc_chk(y_tst , yhat_NANLR)

print('LR Model Train Score is : ' , LR.score(X_trn, y_trn))
print('LR Model Test Score is : ' , LR.score(X_tst, y_tst))


# Support Vector Regressor

SVR = SVR(C = 1.0 ,epsilon=0.1,kernel = 'rbf') # it also can be : linear, poly, rbf, sigmoid, precomputed
SVR.fit(X_trn, y_trn)

yhat_NANSVR = LR.predict(X_tst)

print ('Accuracy scores for filling NAN with Support Vector Regressor model are:')
NANSVR_MAE , NANSVR_MSE , NANSVR_R2_sc = acc_chk(y_tst , yhat_NANSVR)

print('SVR Model Train Score is : ' , SVR.score(X_trn, y_trn))
print('SVR Model Test Score is  : ' , SVR.score(X_tst, y_tst))


# Decision Tree Regressor

DTR = DecisionTreeRegressor( max_depth=3)
DTR.fit(X_trn, y_trn)

yhat_NANDTR = DTR.predict(X_tst)

print ('Accuracy scores for filling NAN with Decision Tree Regressor model are:')
NANDTR_MAE , NANDTR_MSE , NANDTR_R2_sc = acc_chk(y_tst , yhat_NANDTR)

print('DTR Model Train Score is : ' , DTR.score(X_trn, y_trn))
print('DTR Model Test Score is  : ' , DTR.score(X_tst, y_tst))

#%%
# different ML strategy: 
# KNN Imputer
from sklearn.impute import KNNImputer

def refill_kNN(ds):   

    df = ds.to_frame()
    
    # Create an instance of KNNImputer with desired number of neighbors
    imputer = KNNImputer(n_neighbors=2)
    
    # Fit the imputer and transform the DataFrame
    imputed_ds = pd.DataFrame(imputer.fit_transform(df), columns=df.columns, index=ds.index)
    
    return imputed_ds

ds_kNN_refilled_l = []
for column in df_del.columns:
    ds_kNN_refilled_l.append(refill_kNN(df_del[column]))

# df with all gamma refilled series
df_kNN_refilled = pd.concat([ds for ds in ds_kNN_refilled_l], axis=1)
df_kNN_refilled.plot()


#%%
# 3. set to 0
def refill_0(ds):
    ds_filled = ds.copy()
    ds_filled[ds_filled.isnull()] = 0.0
    return ds_filled

ds_0_refilled_l = []
for column in df_del.columns:
    ds_0_refilled_l.append(refill_0(df_del[column]))

# df with all 0 refilled series
df_0_refilled = pd.concat([ds for ds in ds_0_refilled_l], axis=1)
df_0_refilled.plot()


#%%
# change sm_prediction function to use on refilled pc dataframes
# Function to optimise loss factor and calculate predictions in one step
def sm_prediction_2(df_aligned, ds_refilled):
    
    def error_function(lam, sm, pc):
        pc_rescaled = rescale_sm(sm, pc) #rescale precipitation to  m^3/(m^3*100) to avoid unit problem
        sm_pred = sm * lam + pc_rescaled
        return np.sqrt(np.mean((sm_pred - sm) ** 2)) # does not work bc unit-dependent

    initial_guess = 0.5 # most lamda stay at 0.5 if initial guess
    result = minimize(error_function, initial_guess, args=(df_aligned["sm_t_minus_1"], ds_refilled), bounds=[(0, 1)])
    lam = result.x[0]
    

    sm_pred = []

    for i, sm in enumerate(df_aligned["sm_t_minus_1"]):
         if i == 0:
             pred = 0 # mit 0 anfangen als Startwert
             sm_pred.append(pred)
         else:
             pred = sm_pred[-1] * lam + df_aligned["pc_t"].iloc[i]
             sm_pred.append(pred)

    sm_pred = np.array(sm_pred)
    
    df_aligned_2 = df_aligned.copy()
    df_aligned_2['sm_pred'] = sm_pred

    corr_pearson = pearsonr(df_aligned_2["sm"], sm_pred)
    corr_spearman = spearmanr(df_aligned_2["sm"], sm_pred)
    #rmse = np.sqrt(np.mean((sm_pred - df_aligned["sm"])**2)) #unten mit rescaled

    rescaled_sm = rescale_sm(df_aligned_2["sm"].values, df_aligned_2["sm_pred"].values)
    
    rmse = np.sqrt(np.mean((rescaled_sm - df_aligned_2["sm"].values)**2))
    
    df_aligned_2["sm_pred_rescaled"] = rescaled_sm
    
    df_stations = pd.DataFrame({"station": [station_nam], "lamda": [lam], "pearson": [corr_pearson], "spearman": [corr_spearman], "rmse": [rmse],
                               "sm_pred_rescaled": [rescaled_sm], "sm_pred": [df_aligned_2["sm_pred"].values], "sm": [df_aligned_2["sm"].values]})

    return df_aligned, df_stations

#%%
# predict sm with refilled pc dataframes and get correlations

df_1_2, df_2_2 = sm_prediction_2(df_1, df_0_refilled["10% NaN"])
df_2_2["spearman"][0]






