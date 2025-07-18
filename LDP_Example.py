#
#
#
import pandas as pd
import geopandas as gpd
import pyproj
import numpy as np
import pygeodesy as pgd
from io import StringIO
gpd.options.io_engine = "pyogrio"

def dd2DMS( dd, PREC=7, POS=''  ):
    '''convert degree to DMS string'''
    return pgd.dms.toDMS( dd, prec=PREC,pos=POS )

##################################################
LDP = '''+proj=tmerc +lat_0=0.0 +lon_0=101.63333333333334 +k_0=1.000036  +x_0=5000  +y_0=-1641000
+a=6378137.0 +b=6356752.314245179 +units=m +no_defs +type=crs'''
LDP = pyproj.CRS( LDP )

TAB = '''
  BLDG-1   5200.000   5300.000      260.000   783541.304   1647258.039      260.000
  BLDG-2   5225.000   5300.000      260.000   783566.316   1647258.335      260.000
  BLDG-3   5225.000   5400.000      260.000   783565.134   1647358.384      260.000
  BLDG-4   5200.000   5400.000      260.000   783540.122   1647358.088      260.000
'''

COLS = ['Pnt','LDP_E','LDP_N','LDP_Elev', 'UTM_E', 'UTM_N',  'UTM_Elev'] 
df = pd.read_csv( StringIO(TAB), header=None, delim_whitespace=True,names=COLS )
gdf = gpd.GeoDataFrame( df, crs=LDP, geometry=gpd.points_from_xy(df.LDP_E,df.LDP_N) )

gdfUTM47 = gdf.to_crs( 'EPSG:32647' )
gdfWGS84 = gdf.to_crs( 'EPSG:4326' )

def TransLatLng(row):
    Lat_DD,  Lng_DD = row.geometry.y, row.geometry.x
    Lat_DMS, Lng_DMS = dd2DMS(Lat_DD,PREC=5), dd2DMS(Lng_DD,PREC=5)
    return Lat_DMS,Lng_DMS ,Lat_DD,Lng_DD
gdfWGS84[['Lat_DMS','Lng_DMS' , 'Lat_DD', 'Lng_DD']] = \
           gdfWGS84.apply( TransLatLng, axis=1, result_type='expand' ) 


print( gdf.iloc[:,[0,1,2,4,5]].to_markdown(floatfmt=[None,None, '.3f','.3f', '.3f','.3f'] ) )
FMT = [None, None, '.3f','.3f', '.0f' , '.3f','.3f', '.0f', '.3f'  ]
print( gdf.to_markdown(floatfmt=FMT) )
print( gdfUTM47.iloc[:,[0,1,2,3,7]].to_markdown( floatfmt=[None, None, '.3f','.3f', '.0f', '.3f' ]  )   )
print( gdfWGS84.iloc[:,[0,1,2,8,9,10,11]].to_markdown( floatfmt=[None, None, '.3f','.3f', None, None , '.9f', '.9f' ]  )  )
import pdb; pdb.set_trace()



