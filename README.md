# THAAO iwv comparison preliminary

> [!IMPORTANT]  
> This is a copy of reanalysisTHAAO specific for extended iwv from 2000 to 2024.

Install eccodes, cfgrib

https://www.ecmwf.int/sites/default/files/elibrary/2018/18727-cfgrib-easy-and-efficient-grib-file-access-xarray.pdf


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
