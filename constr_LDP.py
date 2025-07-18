# -*- coding: utf-8 -*-
"""  
PROGRAM : Constr_LDP
***Design of a Low Distortion Projection for a Construction Project***</br>  
 Program to design  low distortion projections (LDPs) using conformal map projections (TM) for minimizing linear distortion between projected coordinates eg. UTM grid and the true distance at the surface of the engineering project.</br>  
Phisan Santitamonont,</br>  
Faculty of Engineering, Chulalongkorn University © 2022
*Phisan.Chula@gmail.com*</br>

History:  
  Update :  30 May 2024
"""
import tomllib
import sys,re
import numpy as np
import pandas as pd
import geopandas as gpd
from shapely.affinity import scale
from shapely.geometry import Point,LineString,MultiPoint
import pygeodesy as pgd
import pyproj
import utm
from pathlib import Path
import argparse

TM = '+proj=tmerc +lat_0=0.0 +lon_0={lon_0} +k_0={k_0} +x_0={x_0} +y_0={y_0}'\
     ' +a=6378137.0 +b=6356752.3142 +units=m +no_defs'
LCC = '+proj=lcc +lat_1={lat_0} +lat_0={lat_0} +lon_0=0  +k_0={k_0}'\
      ' +x_0={x_0} +y_0={y_0} +a=6378137.0 +b=6356752.3142 +units=m +no_defs'
COL_LDP = ['UNDUL', 'h','HSF','PSF','CSF', 'CSF_ppm', 'LDP_E', 'LDP_N']

def dd2DMS( dd, PREC=7, POS=''  ):
    '''convert degree to DMS string'''
    return pgd.dms.toDMS( dd, prec=PREC,pos=POS )

def dd2DM( dd ): 
    ''' truncate degree to degree-minute DDD:MM'''
    return pgd.dms.toDMS( dd, prec=0 )[:-3]

def parse_dms( dms ):
    return pgd.parseDMS( dms, sep=':')

class LDP_Design:
    def __init__(self, args ):
        self.ARGS = args
        self.STEM = Path( args.LDP_TOML ).stem
        Path('./CACHE').mkdir( parents=True, exist_ok=True )
        self.GetOFFSET_PP()
        #import pdb ; pdb.set_trace()
        ###########################################################
        self.ELLPS  = pgd.datums.Ellipsoids.WGS84
        TGM_2017 = "/usr/share/GeographicLib/geoids/tgm2017-1.pgm"
        if Path( TGM_2017 ).is_file():
            self.GEOID = pgd.geoids.GeoidKarney( TGM_2017 )
        else:
            self.GEOID = pgd.geoids.GeoidKarney( './tgm2017-1.pgm' )
        self.dfPP = self.LoadTestPoint()
        UNDUL = self.GEOID.height( self.dfPP.lat.mean(),self.dfPP.lng.mean() )
        self.MSL_PP = self.dfPP.MSL.mean() + self.DATA.OFFSET_PP[0]
        self.HAE_PP = UNDUL + self.MSL_PP                 # h = N + H
        RG = self.ELLPS.rocGauss( self.dfPP.lat.mean() )  # RG = sqrt(MN)
        self.k0 = np.round(1 + self.HAE_PP/RG, 6)    #  M.Dennis 2016 : Ground Truth ... (...5 to 6 digits)
        self.CTR_PAR = self.dfPP.lat.mean() 
        ##########################################################
        if (self.DATA['FALSE_EN'])=='AUTO':
            self.CreateLDP( 0, 0 )
            self.dfPP[COL_LDP] = self.dfPP.apply( self.CalcLDP, axis=1, result_type='expand' )
            FalseE,FalseN = self._FindFalse(self.dfPP.LDP_E),self._FindFalse(self.dfPP.LDP_N)
            self.dfPP[['LDP_E','LDP_N']]= self.dfPP[['LDP_E','LDP_N']] + [FalseE,FalseN]
            self.CreateLDP( FalseE,FalseN )
        else:
            FalseE, FalseN = self.DATA['FALSE_EN'] # user defined
            self.CreateLDP( FalseE, FalseN )
            self.dfPP[COL_LDP] = self.dfPP.apply( self.CalcLDP, axis=1, 
                                                        result_type='expand' )

    def GetOFFSET_PP(self):
        with open( self.ARGS.LDP_TOML,"rb") as f :
            self.DATA = pd.Series( tomllib.load( f ) )
        #import pdb; pdb.set_trace()
        if 'OFFSET_PP' in self.DATA.keys():
            self.DATA['OFFSET_PP'] = [ self.DATA['OFFSET_PP'],'defined in TOML']
        else:
            self.DATA['OFFSET_PP'] = [0.0, 'default hPP from average topo or pnts']
        if self.ARGS.OFFSET_PP is not None:  # most prioritized
            self.DATA['OFFSET_PP'] = [ self.ARGS.OFFSET_PP,'defined by CLI args' ]
        print( f'{self.DATA}' )

    def CalcLDP(self, row):
        UNDUL = self.GEOID.height( row.lat,row.lng )
        RG  = self.ELLPS.rocGauss( row.lat )
        h   = UNDUL + row.MSL
        HSF = RG/(RG+h)
        ldp_crs = pyproj.Proj(self.DATA.LDP_CRS)
        PSF = ldp_crs.get_factors( row.lng, row.lat ).meridional_scale
        CSF = PSF*HSF
        CSF_ppm = (CSF-1)*1E6
        TR = pyproj.Transformer.from_crs( 'epsg:4326', self.DATA.LDP_CRS ) 
        LDP_E,LDP_N = TR.transform( row.lat, row.lng )
        return [UNDUL, h, HSF, PSF, CSF, CSF_ppm, LDP_E, LDP_N]
        ### COL_LDP = ['UNDUL', 'h','HSF','PSF','CSF', 'CSF_ppm', 'LDP_E', 'LDP_N']

    def CreateLDP(self, FalseE, FalseN):
        if self.DATA.LDP[0] == 'TM':
            LDP_STR = TM.format( lon_0=parse_dms(self.DATA.LDP[1]),k_0=self.k0,x_0=FalseE,y_0=FalseN )
        elif self.DATA.LDP[0] == 'LCC':
            LDP_STR = LCC.format(lat_0=parse_dms(self.DATA.LDP[1]),k_0=self.k0,x_0=FalseE,y_0=FalseN )
        else:
            print( f'UNKNOWN***LDP_TYPE = {self.LDP} ...'); raise Exception('***ERROR***')
        #import pdb; pdb.set_trace()
        self.DATA['LDP_CRS'] = pyproj.CRS( LDP_STR )

    def _FindFalse(self, dfEN):
        digits = int( np.log10(np.ptp(dfEN))) 
        lo = np.round( dfEN.min(), -(digits-1) ) ; up = np.round( dfEN.max(), -(digits-1) )
        mid = (lo+up)/2
        target = 5*(10**digits)
        return  int( np.round( target-mid,-digits) )

    def LoadTestPoint( self ):
        if 'GPKG' in self.DATA['TEST_POINT'].keys():
            gdfs = list()
            for lay in self.DATA['TEST_POINT']['GPKG'][1].split('|'):
                gdf = gpd.read_file( self.DATA['TEST_POINT']['GPKG'][0], layer=lay )
                gdfs.append( gdf )
            gdf = pd.concat( gdfs, ignore_index=True, copy=True )
            dfPP = gdf.to_crs('EPSG:4326')
            dfPP['lng'] = dfPP.geometry.x; dfPP['lat'] = dfPP.geometry.y
        elif 'XLSX' in self.DATA['TEST_POINT'].keys():
            df = pd.read_excel( self.DATA['TEST_POINT']['XLSX'][0], engine='openpyxl' )
            gdf = gpd.GeoDataFrame( df, crs=self.DATA['TEST_POINT']['XLSX'][1], 
                                        geometry=gpd.points_from_xy( df.UTM_E, df.UTM_N ) )
            dfPP = gdf.to_crs('EPSG:4326')
            dfPP['lng'] = dfPP.geometry.x; dfPP['lat'] = dfPP.geometry.y
        else:  # make 3 planes from buffering a single point
            import pdb ; pdb.set_trace()
            POS  = self.DATA['TEST_POINT']['POS_LATLNG']
            MSL  = self.DATA['TEST_POINT']['MSL']
            HOR_BUF, VER_BUF = self.DATA['TEST_POINT']['BUFFER']
            EWNS =Point( POS[1],POS[0] ).buffer( HOR_BUF/111_000, cap_style=3 ).exterior.coords.xy
            PP_EWNS = np.vstack( (np.array( [POS[1],POS[0]] ) , np.array(EWNS).T) )[:-1]
            dfPP = pd.DataFrame( {'Point':['P0','P1','P2','P3','P4'],
                                  'lng':PP_EWNS[:,0], 'lat':PP_EWNS[:,1]  } )
            dfPP = pd.concat( 3*[dfPP] ,  ignore_index=True) # create 3 planes ...
            dfPP['MSL'] = 5*[MSL+VER_BUF]+  5*[MSL] +  5*[MSL-VER_BUF]
            dfPP = gpd.GeoDataFrame( dfPP, crs='EPSG:4326', geometry=gpd.points_from_xy( dfPP.lng, dfPP.lat ) )
        return dfPP

    def Print_Defintion(self):
        parm=self.DATA.LDP_CRS.to_dict()
        ll_dm = dd2DM(parse_dms(self.DATA.LDP[1]))
        LDP_DEF = self.DATA.LDP_CRS.to_string()
        if self.DATA.LDP[0]=='TM':
            LDP_DEF_= re.sub(r'\+lon_0=[^\s]+', f'+lon_0={ll_dm}', LDP_DEF)
        elif self.DATA.LDP[0]=='LCC':
            LDP_DEF_= re.sub(r'\+lat_0=[^\s]+', f'+lat_0={ll_dm}', LDP_DEF)
            LDP_DEF_= re.sub(r'\+lat_1=[^\s]+', f'+lat_1={ll_dm}', LDP_DEF_)
        def PrintTwoLine( P4STR ):
            LEN = P4STR.find('+a='); print(f'{P4STR[:LEN]}', end=""); print(f'{P4STR[LEN:]}' )
        print(f'================================ LDP Defintion ===============================')
        PrintTwoLine( LDP_DEF)
        PrintTwoLine( LDP_DEF_)
        WKT = self.DATA.LDP_CRS.to_wkt( pretty=True )
        with open( f'CACHE/{self.STEM}_CRS.WKT', 'w' ) as f:
            f.write( WKT+'\n' )
        if self.DATA.LDP[0]=='TM':
            cm_sp = LineString( [ [ parm['lon_0'], self.dfPP.lat.min() ],
                                  [ parm['lon_0'], self.dfPP.lat.max() ] ] )
        elif self.DATA.LDP[0]=='LCC':
            cm_sp = LineString( [ [ self.dfPP.lng.min(), parm['lat_0'] ],
                                  [ self.dfPP.lng.max(), parm['lat_0'] ] ] )
        dfCM = gpd.GeoDataFrame( {'geometry':[ scale(cm_sp,xfact=1.25,yfact=1.25,origin='centroid'),] } )
        dfCM.to_file(      f'./CACHE/{self.STEM}.gpkg', driver='GPKG', layer='CM_SP' )

        #import pdb; pdb.set_trace()
        self.dfPP.to_file( f'./CACHE/{self.STEM}.gpkg', driver='GPKG', layer='LDP_Analysis' )

    def Print_Summary(self):
        print(f'=========================== LDP {self.DATA.LDP[0]} ===========================')
        centr = MultiPoint(self.dfPP.geometry).centroid
        print( f'@ Centroid @     <lng={dd2DM(centr.x)}>    <lat={dd2DM(centr.y)}>' )
        msl = self.dfPP.MSL.describe()
        csf = self.dfPP.CSF_ppm.describe()
        print( f'Offset Projection Plane = {self.DATA.OFFSET_PP[0]:+.1f} m. <{self.DATA.OFFSET_PP[1]}>' )
        print( f'Designed Project Plane : hPP={self.HAE_PP:+.2f} m. / MSL={self.MSL_PP:+.2f} m.'\
               f' tied to k0 = {self.k0:.6f}' )
        print( f'Points CSF min/mean/max [ppm] : {csf["min"]:+.1f} / {csf["mean"]:+.1f} / {csf["max"]:+.1f}'  ) 
        print( f'Points MSL min/mean/max [m.]  : {msl["min"]:+.1f} / {msl["mean"]:+.1f} / {msl["max"]:+.1f}' ) 
        min_max = self.dfPP[['LDP_E', 'LDP_N']].agg(['min','max'])
        min_max.loc['diff'] = min_max.iloc[1]-min_max.iloc[0]
        print( min_max.to_markdown(floatfmt=',.3f' ) )

    def Print_CSFppm(self):
        print(f'===================== CSF_ppm Distribution =====================')
        COLS = ['Point', 'lng', 'lat', 'MSL', 'CSF_ppm', 'LDP_E', 'LDP_N']
        FMT = [ None, None, ',.6f', ',.6f', '+,.1f', '+,.1f', ',.3f', ',.3f', ] 
        print( self.dfPP[COLS].to_markdown( floatfmt=FMT ) )

    def Print_UTM(self):
        epsg = self.dfPP.estimate_utm_crs().to_epsg()
        zn = str(epsg)[-2:]
        dfUTM = self.dfPP.to_crs( epsg )
        #import pdb; pdb.set_trace()
        ####utm_crs = pyproj.Proj(epsg)
        ####PSF = ldp_crs.get_factors( row.lng, row.lat ).meridional_scale
        def makeUTM(row):
            utm_e = row.geometry.x ; utm_n = row.geometry.y 
            return utm_e,utm_n
        dfUTM[[f'UTM{zn}_E',f'UTM{zn}_N']] = dfUTM.apply( makeUTM, axis=1, result_type='expand' )
        FMT = [ None, ',.3f', ',.3f', '.3f', ',.3f', ',.3f' ]
        print( dfUTM[['Point', f'UTM{zn}_E', f'UTM{zn}_N', 'MSL', 'LDP_E', 'LDP_N' ]].to_markdown(
                     index=False , floatfmt=FMT ) )

    def DoTransformation(self,TOML_SECT ):
        FR_PRJ = pyproj.CRS( self.DATA[TOML_SECT]['PROJ'] )
        TO_PRJ = self.DATA.LDP_CRS
        FR_COL = ['UTM_E','UTM_N','UTM_Elev']
        TO_COL = ['LDP_E','LDP_N','LDP_Elev']
        if TOML_SECT=='UTM_LDP':
            pass
        elif TOML_SECT=='LDP_UTM':
            FR_PRJ,TO_PRJ = TO_PRJ,FR_PRJ
            FR_COL,TO_COL = TO_COL,FR_COL
        del self.DATA[TOML_SECT]['PROJ']
        pnts = self.DATA[TOML_SECT]
        df = pd.DataFrame.from_dict( pnts, orient='index', columns=FR_COL )
        TR = pyproj.Transformer.from_crs( FR_PRJ, TO_PRJ )
        def Transf( row,TR, FR_COL,TO_COL ):
            E,N = TR.transform( row[FR_COL[0]], row[FR_COL[1]] )
            return E,N, row[FR_COL[2] ]
        df[TO_COL] = df.apply( Transf, axis=1, result_type='expand',
                                args=(TR,FR_COL,TO_COL) )
        return df

###########################################################################
###########################################################################
###########################################################################
parser = argparse.ArgumentParser( prog='constr_LDP',
                    description='Calculate CSF from points defined by GPKG, CSV or single point',
                    epilog='P.Santitamnont ( phisan.chula@gmail.com ) July,2024) ')
parser.add_argument('LDP_TOML', help="TOML file, LDP parameters e.g. LDP=['TM','99:38']")
parser.add_argument('-c', '--csf', action='store_true', help='show CSF table') 
parser.add_argument('-u', '--utm', action='store_true', help='show UTM-LDP table')
parser.add_argument('-o', '--OFFSET_PP', type=int, help='offset for project plane')
args = parser.parse_args()
print(args)

ldp = LDP_Design(args)
ldp.Print_Summary()
ldp.Print_Defintion()

if args.csf: ldp.Print_CSFppm()
if args.utm: ldp.Print_UTM()
if 1:
    dfLDP = ldp.DoTransformation( 'UTM_LDP' )
    print( dfLDP.to_markdown( floatfmt=",.3f" ) )
    dfUTM = ldp.DoTransformation( 'LDP_UTM' )
    print( dfUTM.to_markdown( floatfmt=",.3f" ) )
############################################################################

