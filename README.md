# THAAO comparison with reanalysis

> [!IMPORTANT]  
> This is a preliminary analysis.

### TODO

- [ ] extract better pixel
- [ ] revise strategy for pixel extraction
- [ ] download ``mean_sea_level_pressure``
- [ ] extract CARRA radiation
- [ ] download more ERA-L products
- [ ] focus on weather variables
- [ ] riportare valori di pressione al livello del mare?
- [ ] recuperare da Giovanni i rs che io non ho e devo aggiungere in archivio per l'analisi

## Reanalysis considered

### CARRA (Copernicus Arctic Regional Reanalysis)

- 3h resolution at 2.5 km
- 10.24381/cds.713858f6
- [cds.climate.copernicus.eu/datasets/reanalysis-carra-single-levels](https://cds.climate.copernicus.eu/datasets/reanalysis-carra-single-levels?tab=overview)

> [!NOTE]
> (from the official website) The C3S Arctic Regional Reanalysis (CARRA) dataset contains 3-hourly analyses and hourly
> short term forecasts of
> atmospheric and surface meteorological variables (surface and near-surface temperature, surface and top of atmosphere
> fluxes, precipitation, cloud, humidity, wind, pressure, snow and sea variables) at 2.5 km resolution. Additionally,
> forecasts up to 30 hours initialised from the analyses at 00 and 12 UTC are available.
> The dataset includes two domains. The West domain covers Greenland, the Labrador Sea, Davis Strait, Baffin Bay,
> Denmark Strait, Iceland, Jan Mayen, the Greenland Sea and Svalbard. The East domain covers Svalbard, Jan Mayen, Franz
> Josef Land, Novaya Zemlya, Barents Sea, and the Northern parts of the Norwegian Sea and Scandinavia.
> The dataset has been produced with the use of the HARMONIE-AROME state-of-the-art non-hydrostatic regional numerical
> weather prediction model. High resolution reanalysis for the Arctic region is particularly important because the
> climate change is more pronounced in the Arctic region than elsewhere in the Earth. This fact calls for a better
> description of this region providing additional details with respect to the global reanalyses (ERA5 for instance). The
> additional information is provided by the higher horizontal resolution, more local observations (from the Nordic
> countries and Greenland), better description of surface characteristics (high resolution satellite and physiographic
> data), high resolution non-hydrostatic dynamics and improved physical parameterisation of clouds and precipitation in
> particular.
> The inputs to CARRA reanalysis are the observations, the ERA5 global reanalysis as lateral boundary conditions and the
> physiographic datasets describing the surface characteristics of the model. The observation values and information
> about their quality are used together to constrain the reanalysis where observations are available and provide
> information for the data assimilation system in areas in where less observations are available.

### ERA5

- 1 h resolution at 0.25° x 0.25°
- 10.24381/cds.143582cf
- [cds.climate.copernicus.eu/datasets/reanalysis-era5-complete](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-complete?tab=overview)

> [!NOTE]
> (from the official website) ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate covering
> the
> period from January 1940
> to present1. It is produced by the Copernicus Climate Change Service (C3S) at ECMWF and provides hourly estimates of a
> large number of atmospheric, land and oceanic climate variables. The data cover the Earth on a 31km grid and resolve
> the atmosphere using 137 levels from the surface up to a height of 80km. ERA5 includes an ensemble component at half
> the resolution to provide information on synoptic uncertainty of its products.
> ERA5 uses a state-of-the-art numerical weather prediction model to assimilate a variety of observations,
> including satellite and ground-based measurements, and produces a comprehensive and consistent view of the Earth's
> atmosphere. These products are widely used by researchers and practitioners in various fields, including climate
> science, weather forecasting, energy production and machine learning among others, to understand and analyse past and
> current weather and climate conditions.

### ERA5-LAND

- 1 h resolution at 0.1° x 0.1°
- 10.24381/cds.e2161bac
- [cds.climate.copernicus.eu/datasets/reanalysis-era5-land](https://cds.climate.copernicus.eu/datasets/reanalysis-era5-land?tab=overview)

> [!NOTE]
> (from the official website) ERA5-Land is a reanalysis dataset providing a consistent view of the evolution of land
> variables over several decades
> at an enhanced resolution compared to ERA5. ERA5-Land has been produced by replaying the land component of the ECMWF
> ERA5 climate reanalysis. Reanalysis combines model data with observations from across the world into a globally
> complete and consistent dataset using the laws of physics. Reanalysis produces data that goes several decades back in
> time, providing an accurate description of the climate of the past.
> ERA5-Land uses as input to control the simulated land fields ERA5 atmospheric variables, such as air temperature and
> air humidity. This is called the atmospheric forcing. Without the constraint of the atmospheric forcing, the
> model-based estimates can rapidly deviate from reality. Therefore, while observations are not directly used in the
> production of ERA5-Land, they have an indirect influence through the atmospheric forcing used to run the simulation.
> In addition, the input air temperature, air humidity and pressure used to run ERA5-Land are corrected to account for
> the altitude difference between the grid of the forcing and the higher resolution grid of ERA5-Land. This correction
> is called 'lapse rate correction'.
> The ERA5-Land dataset, as any other simulation, provides estimates which have some degree of uncertainty. Numerical
> models can only provide a more or less accurate representation of the real physical processes governing different
> components of the Earth System. In general, the uncertainty of model estimates grows as we go back in time, because
> the number of observations available to create a good quality atmospheric forcing is lower. ERA5-land parameter fields
> can currently be used in combination with the uncertainty of the equivalent ERA5 fields. The temporal and spatial
> resolutions of ERA5-Land makes this dataset very useful for all kind of land surface applications such as flood or
> drought forecasting. The temporal and spatial resolution of this dataset, the period covered in time, as well as the
> fixed grid used for the data distribution at any period enables decisions makers, businesses and individuals to access
> and use more accurate information on land states.

# THAAO reference instruments

Instrument involved are HATPRO (IWV, LWP), aws_ECAPAC (temp, press, rh, windd, winds), aws_vespa (temp, press, rh,
windd, winds), VESPA (iwv), sw radiometers (up and down), radiosoundings (iwv)
The reference values are always from THAAO measurements, except for IWV (ref: VESPA) and LWP (ref:HATPRO)

> [!IMPORTANT]
> The code can be run at whichever time resolution.

> [!IMPORTANT]
> The pixel extraction is done before running this code using ``cdo remapnn`` for pixel 76.5, -68.8 which has a few
> hundreds meters offset from the exact location.

> [!TIP]
> wind component combination, conversion from dewpoint temperature to rh and similar, have been performed using the
``metpy`` package.

# Statistics

## Time series

## Scatterplots

> [!WARNING]
> The scatterplots containing all season data is a 2D density histogram, meaning that not all data are
> plotted (there are too many, and they would result overlaid), therefore the colorbar towards the warm color indicates
> higher density of values in that bin. For scatterplots with a limited number of data the result is a few number of
> points compared to the variable N (printed). Seasonal scatterplots are standard.

### Bland-Altman plot

A reference for this type of plots used for atmospheric analysis can be found
here: [Validation of the Cloud_CCI (Cloud Climate Change Initiative) cloud products in the Arctic](https://amt.copernicus.org/articles/16/2903/2023/).

## N

Total number of points analysed.

### MAE

> mae = np.nanmean(np.abs(x - y))

Excluding nan values. x(t): reference value, usually THAAO, see above; y(t): reanalysis or other

### R

Pearson correlation coefficient

### RMSE

> rmse = np.sqrt(np.nanmean((x - y) ** 2)

Excluding nan values. x(t): reference value; y(t): reanalysis or other

# Variables

## WEATHER

### Surface Pressure (``surf_pres``)

- CARRA: ``surface_pressure``
- ERA-5: ``surface_pressure``
- ERA5-L: /
- THAAO (vespa):
- THAAO (aws_ECAPAC):

### Surface temperature (``temp``)

- CARRA: ``2m_temperature``
- ERA-5: ``2m_temperature``
- ERA5-L: ``2m_temperature``
- THAAO (vespa):
- THAAO (aws_ECAPAC):

### Relative Humidity (``rh``)

- CARRA: ``2m_relative_humidity``
- ERA-5: ``2m_dewpoint_temperature`` + ``2m_temperature`` (descrivere processo per ottenere rh)
- ERA5-L: /
- THAAO (vespa):
- THAAO (aws_ECAPAC):

### Wind Direction (``windd``)

- CARRA:``10m_wind_direction``
- ERA-5: ``10m_u_component_of_wind`` + ``10m_v_component_of_wind`` (descrivere processo per ottenere velocità e
  direzione)
- ERA5-L: /
- THAAO (vespa):
- THAAO (aws_ECAPAC):

### Wind Speed (``winds``)

- CARRA: ``10m_wind_speed``
- ERA-5: ``10m_u_component_of_wind`` + ``10m_v_component_of_wind`` (descrivere processo per ottenere velocità e
  direzione)
- ERA5-L: /
- THAAO (aws_ECAPAC):

## RADIATION

> [!WARNING]
> For CARRA radiation values. These forecast variables are released t different leadtimes, with 1-hour frequency. 
> Therefore, we consider only leadtime 1, obtaining every three hours, hourly forecast valued for the following hour 
> w.r.t the chose timeframe. For example, we choose April 1, 2023 at 6:00 UTC, we analyze forecast values on April 1, 2023 
> at 7:00 UTC.

### Downward shortwave irradiance - DSI (``sw_down``)

- CARRA: ``surface_solar_radiation_downwards`` (forecast) + ``surface_net_solar_radiation`` (forecast).
- ERA-5: ``surface_net_solar_radiation`` + ``surface_solar_radiation_downwards``
- ERA5-L: /
- THAAO (pyrgeometers): ``DSI``

### Upward shortwave irradiance - USI (``sw_up``)

- CARRA: ``surface_solar_radiation_downwards`` (forecast) + ``surface_net_solar_radiation`` (forecast).
- ERA-5: ``surface_net_solar_radiation`` + ``surface_solar_radiation_downwards``
- ERA5-L: /
- THAAO (pyrgeometers): ``USI``

### Downward longwave irradiance - DLI (``lw_down``)

- CARRA: ``thermal_surface_radiation_downwards`` (forecast) + ``surface_net_thermal_radiation`` (forecast).
- ERA-5: ``surface_net_thermal_radiation`` + ``surface_thermal_radiation_downwards``
- ERA5-L: /
- THAAO (pyranometers): ``DLI``

### Upward longwave irradiance - ULI (``lw_up``)

- CARRA: ``thermal_surface_radiation_downwards`` (forecast) + ``surface_net_thermal_radiation`` (forecast).
- ERA-5: ``surface_net_thermal_radiation`` + ``surface_thermal_radiation_downwards``
- ERA5-L: /
- THAAO (pyranometers): ``ULI``

### Surface albedo (``alb``)

- CARRA: ``albedo`` (forecast). Values masked to nan for alb<0.1, since they are unrealistic.
- ERA-5: ``forecast_albedo`` (also ``snow_albedo``)
- ERA5-L: /
- THAAO (pyrgeometers): ``DSI``+``USI``

## CLOUD & ATMOSPHERE

## Precipitation (``precip``)

- CARRA: ``total_precipitation``
- ERA-5: ``total_precipitation``
- ERA5-L: /
- THAAO (rain gauge): It is calculated as cumulative value over the resampling time.

## Cloud Base Height (``cbh``)

- CARRA: ``cloud_base``
- ERA-5: ``cloud_base_height``
- ERA5-L: /
- THAAO (ceilometer): ``tcc`` CBH is calculated as the median value over 1 h form the original 15 s time resolution,
  then averaged for the comparison.

## Total Cloud Cover (``tcc``)

- CARRA: ``total_cloud_cover``
- ERA-5: ``total_cloud_cover``
- ERA5-L: /
- THAAO (ceilometer): ``cbh`` (lowermost level)

## Integrated water vapour - IWV (``iwv``)

> [!WARNING]
> - CARRA: ``total_column_integrated_water_vapour``: C'è un problema per CARRA iwv nel dicembre 2023. i valori sono
    tutti nulli. Ho provato a riscaricare ma non cambia. Alla fine ho filtrato i dati <=0.

- ERA-5: ``total_column_water_vapour``
- ERA5-L: /
- THAAO (rs): The vertical integration for rs is **missing**.
- THAAO (vespa):

> [!WARNING]
> - THAAO (hatpro): IWV HATPRO values have been masked to nan for values<0.0 and values>30.

## Liquid Water Path - LWP (``lwp``)

> [!CAUTION]  
> LWP values have issues, at least for CARRA which has been divided by 10E-06 instead of 10e-03 as expected from the
> declared uom. All LWP values have been masked to nan for LWP<0.0, and to 0 for LWP<15.

- CARRA: ``total_column_cloud_liquid_water``
- ERA-5: ``total_column_cloud_liquid_water`` (also ``total_column_water``?)
- ERA5-L: /
- THAAO (hatpro):