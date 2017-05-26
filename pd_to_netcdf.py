# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:35:49 2017

@author: wpreimes
"""

from __future__ import print_function
#Import landgrid
import types

import warnings
import pygeogrids.netcdf as nc
from pygeogrids.grids import BasicGrid
import os
import rsdata.root_path as root
import numpy as np
import pandas as pd
from datetime import datetime
from netCDF4 import Dataset,date2num
import matplotlib.pyplot as plt
from HomogeneityTesting.import_data import QDEGdata_M

def get2Dpos(gpis,globalgrid,landgrid):
    grid_points=globalgrid.get_grid_points()
    lats=np.unique(grid_points[2])[::-1]
    lons=np.unique(grid_points[1])
    lon,lat=landgrid.gpi2lonlat(gpis)
    y=datamask(lons,lon)
    x=datamask(lats,lat)
    return x,y

def globalCellgrid():

    # symettrical grid
    lon_res = 0.25
    lat_res = 0.25
    offs_h=lon_res/2.
    offs_v=lat_res/2.
    # create meshgrid
    lon, lat = np.meshgrid(
        np.arange(-180+offs_h, 180-offs_h , lon_res),
        np.arange(-90+offs_v, 90-offs_v , lat_res)
    )

    return BasicGrid(lon.flatten(), lat.flatten()).to_cell_grid(cellsize=5.)    
    
def datamask(x,y):
    index = np.argsort(x)
    sorted_x = x[index]
    sorted_index = np.searchsorted(sorted_x, y)
    
    yindex = np.take(index, sorted_index, mode="clip")
    mask = x[yindex] != y
    
    result = np.ma.array(yindex, mask=mask)
    return result

def create_cellfile_name(gpi,grid):
    #Returns filename (form:cellnumber.nc) and cell for passed gpi in passed grid
    grid_points=grid.get_grid_points()
    gpi_index = np.where(grid_points[0]==gpi)[0][0]
    cell=grid_points[3][gpi_index]
    #Create filename from cell
    file_pattern = str(cell)
    while len(file_pattern)<4:
        file_pattern=str(0)+file_pattern
        
    return cell,file_pattern  

def update_time(ncfile,time,overwrite):     
    calendar = 'standard'
    units = 'days since 1858-11-17 00:00'
    
    nctime=date2num(time.tolist(),units,calendar)  
    
    if u'time' in ncfile.variables.keys():
        
        DF_Time=pd.DataFrame()
        
    
        if np.isnan(ncfile.variables['time'][:]).all():
            overwrite=True
            #ncfile.variables['time'][:]=nctime
        
        #append new time data to old time data (if exists)
        old_values=ncfile.variables[u'time'][:]
        
        if overwrite==True:
            #File mask: Time values in file time variable that are not needed for curent dataframe --> set True
            #Frame mask: Time values that are not yet in the file and new for processsing the dataframe --> set False
            DF_Time=pd.DataFrame(data={'file_values':np.sort(old_values),'file_mask':np.in1d(old_values,nctime,invert=True),
                                       'frame_values':np.sort(nctime),'frame_mask':np.in1d(nctime,old_values,invert=False)}) #np.unique(np.concatenate((old_values,nctime))
                                
        if overwrite==False:
            DF_Time=pd.DataFrame(data={'file_values':np.sort(old_values),'file_mask':np.in1d(old_values,nctime,invert=False),
                                       'frame_values':np.sort(nctime),'frame_mask':np.in1d(nctime,old_values,invert=True)}) #np.unique(np.concatenate((old_values,nctime))
                                 
        frame_values=DF_Time['frame_values'].where(DF_Time['frame_mask']==True).dropna()
        file_values=DF_Time['file_values'].where(DF_Time['file_mask']==True).dropna()
        
        DF_Time['merged']=np.unique(np.concatenate((file_values,frame_values)))
            
        DF_Time['sort_order']=np.argsort(DF_Time['merged'].values)
        contt=DF_Time['merged'].values[DF_Time['sort_order'].values]
        if frame_values.size!=0:
            ncfile.variables['time'][:]=contt
        
        return DF_Time
        
    else:
        tvar=ncfile.createVariable(varname=u'time',
                                    datatype=nctime.dtype,
                                    dimensions=(u'time',),
                                    zlib=False)
        
        tvar[:]=np.full(nctime.size,np.nan)
        
        
        return update_time(ncfile,time,overwrite)
    
    


def update_time_var(ncfile,DF_Time,data,name,idx):

    if name in ncfile.variables.keys():
        #Old Values werden nicht gebraucht, new values mask beschreibt wo neu Zeitwerte hinzugeommen sind
        #im Prinzip sind einfach alle neuen Werte zu übernehmen, falls die Anzahl gleich der Zeitwerte ist,
        #ansonsten müssen für Zeitwerte, die nicht für die New Values gelten, Nans stehen (OLD VALUES MASK == False???)
              

        
        #Overwriting existing data:
        DF_Var=pd.DataFrame(index=ncfile.variables['time'][:],
                            data={'file_values':ncfile.variables[name][:][idx],
                                  'frame_values':data.values})  
        '''
        if new_data_dict['overwrite']==False:
            #If overwrite is not selected, keep data for existing times and just add new ones
            old_contt=np.ma.masked_array(ncfile.variables[name][:][idx],new_data_dict['old_values_mask'])
            new_contt=np.ma.masked_array(data.values,new_data_dict['add_values_mask'])
        if new_data_dict['overwrite']==True:
            pass
            #TODO
        '''
        #merged=np.concatenate((old_contt.compressed(),new_contt.compressed()))[new_data_dict['order']]
        
        contt=ncfile.variables[name][:]
        contt[idx]=merged
        ncfile.variables[name][:]=contt
        

    else:
        try:
            loc_size=ncfile.dimensions['locations'].size
            if isinstance(data, str):
                dtype=object
                contt=np.empty((loc_size,data.size,))
                contt[:]=''
            else:
                dtype=np.asarray(data).dtype
                #contt=[[np.nan]*data.size]*loc_size
                contt=np.empty((loc_size,data.size,))
                contt[:]=np.nan
            ncvar=ncfile.createVariable(varname=name,
                                        datatype=dtype,
                                        dimensions=(u'locations',u'time'),
                                        zlib=False)
            ncvar[:]=np.array(contt)
        except Exception:
            print('Cannot save data for %s to file'%name)
        
        return update_time_var(ncfile,new_data_dict,data,name,idx)
        

def update_loc_var(ncfile,data,name,grid,idx):
    
    if name in ncfile.variables.keys():
        contt=ncfile.variables[name][:]
        if idx!=None:
            if contt.ndim==2:
                x,y=get2Dpos(idx,grid[0],grid[1])
                contt[x,y]=data
            else:
                contt[idx]=data
        else:
            contt=data
        ncfile.variables[name][:]=contt
    else:
        if name in ncfile.dimensions.keys():
            dimension=[ncfile.dimensions[name]]
            dimension_size=dimension[0].size
        elif name == 'location_id':
            dimension=[ncfile.dimensions['locations']]
            dimension_size=dimension[0].size 
        elif  u'lat' and u'lon' in ncfile.dimensions:
            dimension=[ncfile.dimensions[u'lat'],ncfile.dimensions[u'lon']]
            dimension_size=dimension[0].size*dimension[1].size
        else:
            dimension=[ncfile.dimensions['locations']]
            dimension_size=dimension[0].size
        #If variable does not exist, create it with correct size and retry
        try:
            if isinstance(data, str):
                dtype=str                
                contt=np.array(['']*dimension_size,dtype=object)
            else:
                contt=np.full(dimension_size,np.nan)
                if np.asarray(data).dtype == int:
                    dtype=float
                else:
                    dtype=np.asarray(data).dtype
            ncvar=ncfile.createVariable(varname=name,
                                        datatype=dtype,
                                        dimensions=tuple([dim.name for dim in dimension]),
                                        zlib=False)
            ncvar[:]=contt
        except Exception:
            print('Cannot save data for %s to file'%name)
        
        return update_loc_var(ncfile,data,name,grid,idx) 
        
        

def fill_file(ncfile,grid,cell,gpi,dataframe,file_meta_dict,var_meta_dicts,index_col_name=None,overwrite=False,istime=True):
          
    if istime:
        grid_points=grid.grid_points_for_cell(cell)
        latitudes,longitudes=grid_points[2][::-1],grid_points[1]
        locations=grid_points[0]
        if gpi and not index_col_name:
            time=dataframe.index
        elif gpi and index_col_name:
            time=dataframe[index_col_name]
        else: 
            raise Exception, 'To save time series data define GPI and time column (if not DF index)'
        pos=None
        glob_pos=None
        
    else:
        globgrid=globalCellgrid()
        grid_points=grid.get_grid_points()
        global_grid_points=globgrid.get_grid_points()
        
        latitudes,longitudes=np.unique(global_grid_points[2])[::-1],np.unique(global_grid_points[1])
        locations=grid_points[0]
        if index_col_name:
            locs=dataframe[index_col_name]
        else:
            locs=dataframe.index
        #glob_pos contains the indices of points to process in the overall grid
        pos=datamask(np.array(locations),np.array(locs))
            
    n_gpis=locations.size
    
    #Create data dimensions for Time series and global image
    if not ncfile.dimensions:
        ncfile.createDimension(dimname='locations',size=n_gpis)
        if istime:
            ncfile.createDimension(dimname=u'time',size=None) 
        else:
            ncfile.createDimension(dimname='lat',size=latitudes.size)
            ncfile.createDimension(dimname='lon',size=longitudes.size)
    #TODO: Add Metadata for netcdf file to dict
    if not ncfile.ncattrs():
        meta_dict={'geospatial_lon_min':longitudes[0],
                   'geosptial_lon_max':longitudes[-1],
                   'geospatial_lat_min':latitudes[-1],
                   'geospatial_lat_max':latitudes[0],
                   'id':cell,
                   'date_created':datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        if file_meta_dict:
            meta_dict.update(file_meta_dict)
        ncfile.setncatts(meta_dict) 


    #Create variable for locations and add value
    #GPI, LAT, LON werden beim erstellen immer gefüllt je nach grid unabhängig vom GPI
    #Statt None: gpi_index: Nur für den prozessierten gpi werden idx,lat,lon ins file gespeichert
    meta={'long_name':'Location Index','standard_name':'GPI','valid_range':'[0 Grid Dependant'}
    update_loc_var(ncfile,locations,u'location_id',grid,pos)
    meta={'units':'degrees_east','long_name':'location longitude','standard_name':'longitude','valid_range':'[-180. 180.]'}
    update_loc_var(ncfile,longitudes,u'lon',grid,None)
    ncfile.variables[u'lon'].setncatts(meta)
    meta={'units':'degrees_north','long_name':'location latitude','standard_name':'latitude','valid_range':'[-90. 90.]'}
    update_loc_var(ncfile,latitudes,u'lat',grid,None)
    ncfile.variables[u'lat'].setncatts(meta)

               
    if istime==True:
        gpi_index = np.where(grid_points[0]==gpi)[0][0]
        DF_Time=update_time(ncfile,time,overwrite)
        ncfile.variables['time'].setncatts({'units':'days since 1858-11-17 00:00:00',
                                            'long_name':'time of measurement',
                                            'standard_name':'time'})
        for i,var in enumerate(dataframe.columns.values):
            
        #Create vaiable for time and other variables in dataframe
            update_time_var(ncfile,DF_Time,dataframe[var],var,gpi_index)
            try:
                ncfile.variables[var].setncatts({'units':'days since 1858-11-17 00:00:00',
                                                 'long_name':'time of measurement',
                                                 'standard_name':'time'})
                ncfile.variables[var].setncatts(var_meta_dicts[var])
            except KeyError:
                #Catch if no metadata was specified and create some trivial one
                #TODO: Make more useful auto meta data
                var_meta_auto={'name':var,'info':'Automatically generated meta data'}
                ncfile.variables[var].setncatts(var_meta_auto)  
    else:
        for i,var in enumerate(dataframe.columns.values):
            glob_pos=datamask(global_grid_points[0],locs.values)
            update_loc_var(ncfile,dataframe[var].values,var,[globgrid,grid],glob_pos)
            try:
                ncfile.variables[var].setncatts(var_meta_dicts[var])
            except KeyError:
                ##TODO: Make more useful auto meta data
                var_meta_auto={'name':var,'info':'Automatically generated meta data'}
                ncfile.variables[var].setncatts(var_meta_auto)
                


def df_to_netcdf(dataframe,
                 path,
                 index_col_name=None,
                 gpi=None,
                 filename=None,
                 file_meta_dict=None,
                 var_meta_dicts=None,
                 overwrite=False):
    
    '''
    Writes Pandas Dataframe to netCDF file using lat/lon for the given gpi
    2 ways of using:
        1) Write time dependant data (data series, data frame) to file:
            -pandas object must contain datetime information as index or in 
            selected column (index_col_name)
            -GPI must be selected
        2) Write spatial data (data series, data frame) to file:
            -pandas object must contain GPIs as index or in the selected 
            column (index_col_name)
    Parameters
    ----------
    dataframe (mandatory): pandas data frame or data series
        pandas object with data for writing to file
        for time series data: date time as index
        for spatial data: gpi as index
    path (mandatory): string
        path where netcdf file is saved to
    index_col_name (optional): string
        name of the column with time/location data in the pandas object
    gpi (optional): integer or tuple
        for time series data: index of the ground point which is processed or
        (lat,lon) for finding nearest ground point
    filename (optional): string
        for time series data: filename is automatically "*cell*.nc"
        for spatial data: select file name
    file_meta_dict (optional): dictionary
        additional meta information on the netcdf file
    var_meta_dict (optional): dictionary of dictionaries
        additional meta information on the written variables
        for each column in the dataframe, there is 1 dictionary in this list
    overwrite (optional): boolean
        If a (cell)file already exits at the chosen location, existing ground 
        point data is overwritten
    '''    

    if not index_col_name:
        index_col=dataframe.index
    else:
        index_col=dataframe[index_col_name]
    
    if isinstance(index_col, pd.DatetimeIndex):
        if gpi:
            istime=True
        else:
            raise Exception, 'Found time data in pandas object, select a GPI!'
    elif isinstance(index_col, pd.Int64Index):
        if gpi:
            raise Exception, 'Found GPI data in dataframe, select GPI=None and provide GP index as dataframe index'
        istime=False
    else:
        raise Exception, 'Index Column must either be datetime format for time series or int64 to process multiple GPIs'
    
    if istime==False:
        gpi=index_col.values

    grid=nc.load_grid(os.path.join(root.r,'Datapool_processed','GLDAS','GLDAS_NOAH025_3H.020',
                                       'ancillary','GLDASv2_025_land_grid.nc'))
    
    #Find cell and create filename for given gpi
    if not filename:
        cell,filename=create_cellfile_name(gpi,grid)
    else: 
        cell,filename='global',filename
    
    #Create or open netcdf cell file
    if os.path.isfile(os.path.join(path,filename+'.nc')):
        ncfile=Dataset(os.path.join(path,filename+'.nc'), "a", format="NETCDF4")
    else:
        ncfile=Dataset(os.path.join(path,filename+'.nc'), "w", format="NETCDF4")
    try:
        fill_file(ncfile,grid,cell,gpi,dataframe,file_meta_dict,var_meta_dicts,index_col_name,overwrite,istime)
    except Exception:
        #TODO: handle the case that no metadata was passed
        #print('Error during filling file %s'%filename)
        pass
    ncfile.close() 
            
    

        
  
def gotest(testtype):
    #TODO: The variable-meta-dict is not good

    if testtype=='time':
        gpi_file=r"H:\workspace\HomogeneityTesting\csv\pointlist_United_457_quarter.csv"
        df=pd.read_csv(gpi_file,index_col=0)
        ttime=['2005-07-01','2011-10-01','2012-07-01']
        data=QDEGdata_M(products=['merra2','cci_22'])
        for i, gpi in enumerate(df.index.values):
            print('Writing gpi %i of %i to netcdf'%(i,df.index.values.size))
            dataframe_time=data.read_gpi(gpi,ttime[0],ttime[1])
        
            var1_meta={'longname':'Soil Moisture','units':'kg/m^2'}
            var2_meta={'longname':'Soil Moisture','units':'kg/m^2'}

                         
            df_to_netcdf(dataframe=dataframe_time,
                         path=r'D:\users\wpreimes\datasets\ncwriter\newfiles',
                         index_col_name=None,
                         var_meta_dicts={'merra2':var1_meta,'cci_22':var2_meta},
                         gpi=gpi,
                         overwrite=False)
                         
    elif testtype=='points':
        
        dataframe_gpi=pd.read_csv(r"H:\workspace\HomogeneityTesting\output\global_merra2_cci22D\DF_Points_merra2_1998-01-01_2002-07-01_2007-01-01.csv",index_col=0)
        
        df_to_netcdf(dataframe=dataframe_gpi[['h_all','h_FK','h_WK']],
                     path=r'D:\users\wpreimes\datasets\ncwriter\newfiles',
                     filename='testfile',
                     gpi=None,
                     index_col_name=None,
                     file_meta_dict={'metainfo': 'value','metainfo2':'value2'},
                     var_meta_dicts={'lat':{'longname':'Latitude','units':'degree'},
                                     'h_all':{'longname':'Break Test Class','units':'class','additional info':'awesome!'}},
                     overwrite=False)
    else:
        print("use gotest('time') or gotest(points) for testing this script")
            
        
                         
gotest('time')
    
    