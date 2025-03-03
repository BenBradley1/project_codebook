import numpy as np
import xarray as xr
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import hvplot.xarray

#ds = xr.open_dataset("./ptom_t021_2010010418.nc") # 1 month
ds = xr.open_dataset("./ptom_t021_2010120412.nc") # 1 year ?


plot = ds.CO.hvplot.quadmesh(x='lon', y='lat', projection=ccrs.PlateCarree(), coastline=True, widget_location='bottom', width=1500, clim=(0,1500))
#plot = ds.CO.hvplot.image(x='lon', y='lat', projection=ccrs.PlateCarree(), coastline=True, widget_location='bottom', width=1500, clim=(0,vmax))
hvplot.show(plot)

#--- add slider (did not work) ---
#import panel as pn
#vmax = pn.widgets.IntSlider(name="vmax", start=0, end=2000, step=100, value=1000)
#graph = pn.Column(vmax, plot).show()