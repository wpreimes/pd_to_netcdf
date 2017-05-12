# -*- coding: utf-8 -*-
"""
Created on Thu Apr 13 14:35:49 2017

@author: wpreimes
"""


#Import landgrid
import types

import warnings
import pygeogrids.netcdf as nc
import os
import rsdata.root_path as root
import numpy as np
import pandas as pd
from datetime import datetime
from netCDF4 import Dataset,date2num
import matplotlib.pyplot as plt
from HomogeneityTesting.import_data import QDEGdata_M

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


#TODO: merge these 2 functions to 1
def update_time_var(ncfile,time,data,name,idx,meta=None):
    
    calendar = 'standard'
    units = 'days since 1858-11-17 00:00'
    loc_size=ncfile.dimensions['locations'].size
    
    nctime=date2num(time.tolist(),units,calendar)
    if data.size != time.size:
        raise Exception, 'Number of time values and variable values do not match'
    
    #TODO: Shorten this
    if u'time' in ncfile.variables.keys():
        contt=ncfile.variables[u'time'][:]
        if any(contt==False):
            raise Exception, 'Time period for variable in file does not match time period from dataframe'
    else:
        tvar=ncfile.createVariable(varname=u'time',
                                    datatype=nctime.dtype,
                                    dimensions=(u'time',),
                                    zlib=False)
        #ncvar.setncatts({'units':})
        #contt=np.full(size,np.nan)
        contt=nctime
        tvar[:]=np.array(contt)
        tvar.setncatts({'units':'days since 1858-11-17 00:00:00','long_name':'time of measurement','standard_name':'time'})
        
    if name in ncfile.variables.keys():
        contt=ncfile.variables[name][:]
        contt[idx]=data
        ncfile.variables[name][:]=contt

    else:
        try:
            if isinstance(data, str):
                dtype=object
                contt=[['']*data.size]*loc_size
            else:
                dtype=np.asarray(data).dtype
                contt=[[np.nan]*data.size]*loc_size

            ncvar=ncfile.createVariable(varname=name,
                                        datatype=dtype,
                                        dimensions=(u'locations',u'time'),
                                        zlib=False)
    
        except Exception:
            print('Cannot save data for %s to file'%name)
        
        contt[idx]=data.tolist()
        dataarray=np.array(contt,dtype=dtype)
        ncvar[:]=dataarray
        
        if meta: ncvar.setncatts(meta)



def update_loc_var(ncfile,data,name,idx,meta=None):
    
    if name in ncfile.dimensions.keys():
        dimension=[ncfile.dimensions[name]]
        size=dimension[0].size
    elif name == 'location_id':
        dimension=[ncfile.dimensions['locations']]
        size=dimension[0].size 
    elif  u'lat' and u'lon' in ncfile.dimensions:
        dimension=[ncfile.dimensions[u'lat'],ncfile.dimensions[u'lon']]
        size=dimension[0].size*dimension[1].size
    else:
        dimension=[ncfile.dimensions['locations']]
        size=dimension[0].size

    if name in ncfile.variables.keys():
        if idx!=None:
            contt=ncfile.variables[name][:]
            if data.size==1:
                contt[idx]=data
            else:
                for i in idx:
                    contt[i]=data[i]
        else:
            contt=data
        ncfile.variables[name][:]=contt
    else:
        try:
            if isinstance(data, str):
                dtype=str                
                contt=np.array(['']*size,dtype=object)
            else:
                contt=np.full(size,np.nan)
                dtype=np.asarray(data).dtype
            ncvar=ncfile.createVariable(varname=name,
                                    datatype=dtype,
                                    dimensions=tuple([dim.name for dim in dimension]),
                                    zlib=False)
            ncvar[:]=contt
        except Exception:
            print('Cannot save data for %s to file'%name)
        
        if meta:
            ncvar.setncatts(meta)
        
        update_loc_var(ncfile,data,name,idx,meta) 
        
        

def fill_file(ncfile,grid,cell,gpi,dataframe,file_meta_dict,var_meta_dicts,index_col_name=None,istime=True):
          
    if istime:
        grid_points=grid.grid_points_for_cell(cell)
        latitudes=grid_points[2]
        longitudes=grid_points[1]
        if gpi and not index_col_name:
            time=dataframe.index
        elif gpi and index_col_name:
            time=dataframe[index_col_name]
        else: 
            raise Exception, 'To save time series data define GPI and time column (if not DF index)'
        glob_pos=None
        
    else:

        grid_points=grid.get_grid_points()
        latitudes=np.unique(grid_points[2])
        longitudes=np.unique(grid_points[1])
        locations=grid_points[0]
        if index_col_name:
            locs=dataframe[index_col_name]
        else:
            locs=dataframe.index
            
    n_gpis=grid_points[0].size
    
    #Create 2 data dimensions 
    if not ncfile.dimensions:
        ncfile.createDimension(dimname='locations',size=n_gpis)
        ncfile.createDimension(dimname=u'time',size=None) 
        if not istime:
            ncfile.createDimension(dimname='lat',size=latitudes.size)
            ncfile.createDimension(dimname='lon',size=longitudes.size)
    #TODO: Add Metadata for netcdf file to dict
    if not ncfile.ncattrs():
        meta_dict={'geospatial_lon_min':longitudes[0],
                   'geosptial_lon_max':longitudes[-1],
                   'geospatial_lat_min':latitudes[0],
                   'geospatial_lat_max':latitudes[-1],
                   'id':cell,
                   'date_created':datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
        if file_meta_dict:
            meta_dict.update(file_meta_dict)
        ncfile.setncatts(meta_dict) 


    #Create variable for locations and add value
    #GPI, LAT, LON werden beim erstellen immer gefüllt je nach grid unabhängig vom GPI
    #Statt None: gpi_index: Nur für den prozessierten gpi werden idx,lat,lon ins file gespeichert
    glob_pos=datamask(np.array(locations),np.array(locs))
    update_loc_var(ncfile,locations,u'location_id',glob_pos,meta=None)
    meta={'units':'degrees_east','long_name':'location longitude','standard_name':'longitude','valid_range':'[-180. 180.]'}
    update_loc_var(ncfile,longitudes,u'lon',None,meta)
    meta={'units':'degrees_north','long_name':'location latitude','standard_name':'latitude','valid_range':'[-90. 90.]'}
    update_loc_var(ncfile,latitudes,u'lat',None,meta)

               
    if istime==True:
        gpi_index = np.where(grid_points[0]==gpi)[0][0]
        for i,var in enumerate(dataframe.columns.values):
        #Create vaiable for time and other variables in dataframe
            update_time_var(ncfile,time,dataframe[var],var,gpi_index)
            try:
                ncfile.variables[var].setncatts(var_meta_dicts[var])
            except KeyError:
                #Catch if no metadata was specified and create some trivial one
                #TODO: Make more useful auto meta data
                var_meta_auto={'name':var,'info':'Automatically generated meta data'}
                ncfile.variables[var].setncatts(var_meta_auto)
                
    else:
        for i,var in enumerate(dataframe.columns.values):
            for gpi in locs.values:
                gpi_index = np.where(grid_points[0]==gpi)[0][0]
                glob_pos=datamask(np.array(locations),gpi)
                update_loc_var(ncfile,dataframe.ix[gpi][var],var)
            try:
                ncfile.variables[var].setncatts(var_meta_dicts[var])
            except KeyError:
                #TODO: Make more useful auto meta data
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
            raise Exception, 'Found time data in pandas object, select GPI!'
    elif isinstance(index_col, pd.Int64Index):
        istime=False
    else:
        raise Exception, 'Index Column must either be datetime format or int64 (GPIs)'
    
    grid=nc.load_grid(os.path.join(root.r,'Datapool_processed','GLDAS','GLDAS_NOAH025_3H.020',
    'ancillary','GLDASv2_025_land_grid.nc'))

        
    if gpi!=None and (istime==True or gpi[0]==0):
        try:
            if gpi.__len__()==2:
                gpi, _= grid.find_nearest_gpi(gpi[0],gpi[1])
            if gpi[0]==0:
                np.delete(gpi,0)
        except Exception:
            pass
        
        #Find cell and create filename for given gpi
        if not filename:
            cell,filename=create_cellfile_name(gpi,grid)
        else: 
            cell,filename='global',filename
        
        #Create or open netcdf cell file
        if os.path.isfile(os.path.join(path,filename+'.nc')) and overwrite==False:
            ncfile=Dataset(os.path.join(path,filename+'.nc'), "a", format="NETCDF4")
        else:
            ncfile=Dataset(os.path.join(path,filename+'.nc'), "w", format="NETCDF4")
        try:
            fill_file(ncfile,grid,cell,gpi,dataframe,file_meta_dict,var_meta_dicts,index_col_name,istime)
        except Exception:
            #TODO: handle the case that no metadata was passed
            #print('Error during filling file %s'%filename)
            pass
        ncfile.close() 
    elif gpi==None and istime==False:
        gpi=np.append(0,index_col.values)
        df_to_netcdf(dataframe,path,index_col_name,gpi,filename,file_meta_dict,var_meta_dicts,overwrite)
    else:
        raise Exception, 'DataFrame cannot be saved to NetCDF, check index, gpi and filename input'
    
    

        
  
def gotest(testtype):

    #TODO: The variable-meta-dict is not good
    gpi_file=r"H:\workspace\HomogeneityTesting\csv\pointlist_United_457_quarter.csv"
    df=pd.read_csv(gpi_file,index_col=0)
    ttime=['2007-07-01','2011-10-01','2012-07-01']
    if testtype=='time':
        data=QDEGdata_M(products=['merra2','cci_22'])
        for gpi in df.index.values:
            
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
        dataframe_gpi=pd.read_csv(r"H:\workspace\HomogeneityTesting\output\v22\DF_Points_merra2_2011-10-01_2012-07-01_2015-05-01.csv",index_col=0)
        dataframe_gpi=dataframe_gpi[1:1000]
        
        var1_meta={'longname':'Latitude','units':'degree'}        
        var2_meta={'longname':'Break Test Class','units':'class','additional info':'awesome!'}
        
        df_to_netcdf(dataframe=dataframe_gpi[['h_all','h_WK','h_FK']],
                     path=r'D:\users\wpreimes\datasets\ncwriter\newfiles',
                     filename='testfile',
                     gpi=None,
                     index_col_name=None,
                     file_meta_dict={'metainfo': 'value','metainfo2':'value2'},
                     var_meta_dicts={'lat':var1_meta,'lon':var2_meta},
                     overwrite=False)
    else:
        print("use gotest('time') or gotest(points) for testing this script")
            
        
                         
gotest('points')
    
    