# pd_to_netcdf
## Code for writing pandas series to netcdf image and time series

**Time Series to netCDF**  
-Time information in index or in selected column  
-GPI must be selected  
-path must be selected  
-file meta data can be passed as dictionary (not tested) {'metainfo': 'value','metainfo2':'value2'}  
-variable meta data can be passed as dict of dicts {{varname:{MetaDict},varname:{MetaDict},...}  
-overwrite (not tested): filling up "missing" data in existing file, adding new variables, keeps existing unchanged, continues exisiting time series  
  
**GPI Series to netCDF**  
-GPIs or [lat,lon] (not tested) in index or selected column  
-GPI must NOT be selected  
-path must be selected  
-file meta data can be passed as dictionary (not tested) {'metainfo': 'value','metainfo2':'value2'}  
-variable meta data can be passed as dict of dicts {{varname:{MetaDict},varname:{MetaDict},...}  
-overwrite: filling up "missing" data in existing file, adding new variables, keeps existing unchanged  



