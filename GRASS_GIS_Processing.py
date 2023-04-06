# -*- coding: utf-8 -*-
"""
@author: ailing

"""

#!/usr/bin/env python

import os
import sys
import subprocess
import shutil
import binascii
import tempfile
import time

################################################   
####  Setup GRASS GIS Environment
def Setup_GRASS():
    #### Configure GRASS Environment
    # path to the GRASS GIS launch script
    # WinGRASS installer
    grass7bin_win = r'C:\Program Files\GRASS GIS 7.8\grass78.bat'
    # Linux
    grass7bin_lin = 'grass78'
    # Mac OS X
    grass7bin_mac = '/Applications/GRASS/GRASS-7.8.app/'
    ########### SOFTWARE
    if sys.platform.startswith('linux'):
        # we assume that the GRASS GIS start script is available and in the PATH
        # query GRASS 7 itself for its GISBASE
        grass7bin = grass7bin_lin
    elif sys.platform.startswith('win'):
        grass7bin = grass7bin_win
    else:
        raise OSError('Platform not configured.')
        
    # Create a tempory GRASS Database
    gisdb = os.path.join(tempfile.gettempdir(), 'grassdata')
    try:
        os.stat(gisdb)
    except:
        os.mkdir(gisdb)
    
    # Query GRASS 7 itself for its GISBASE
    startcmd = [grass7bin, '--config', 'path']
    
    p = subprocess.Popen(startcmd, shell=False,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    if p.returncode != 0:
        print("ERROR: Cannot find GRASS GIS 7 start script (%s)" % startcmd, file=sys.stderr)
        sys.exit(-1)
    gisbase = out.decode('ascii').strip('\n\r')
    # Set GISBASE environment variable
    os.environ['GISBASE'] = gisbase
    
    osPath = os.environ['PATH']
    # the following not needed with trunk
    if osPath.find(os.path.join(gisbase, 'extrabin')) == -1:
        os.environ['PATH'] += os.pathsep + os.path.join(gisbase, 'extrabin')
        
    # add path to GRASS addons
    home = os.path.expanduser("~")
    if osPath.find(os.path.join(home, '.grass7', 'addons', 'scripts')) == -1:
        os.environ['PATH'] += os.pathsep + os.path.join(home, '.grass7', 'addons', 'scripts')
        
    # Set GISDBASE environment variable
    os.environ['GISDBASE'] = gisdb    
     
    # define GRASS-Python environment
    gpydir = os.path.join(gisbase, "etc", "python")
    if gpydir not in sys.path:
        sys.path.append(gpydir)
    
    # For windows, probably need to explicitly set PYTHONHOME
    PyPath = os.getenv('PYTHONHOME')
    Pydirr = os.path.join(gisbase, 'Python37')
    if not PyPath:
        PyPath = Pydirr
    elif PyPath.find(Pydirr) == -1:
        PyPath = Pydirr + os.pathsep + PyPath
        
    os.environ['PYTHONHOME'] = PyPath
    return grass7bin, gisdb, gisbase

    
################################################   
### GRASS GIS r.watershed analysis
def GRASS_WatershedAnalysis(demraw_Path, threshold, out_dir):
    strT = time.time()
    # Set up GIS Environment
    grass7bin, gisdb, gisbase = Setup_GRASS()

    ### Paths to output files
    flwdir_Path = os.path.join(out_dir,'flwdir.tif')
    subbasin_Path = os.path.join(out_dir,'subbasins.tif')
    strseg_Path = os.path.join(out_dir,'strseg.tif')
    flwacc_Path = os.path.join(out_dir,'flwacc.tif')
    strvector_Path = os.path.join(out_dir,'streams.shp')
    
    ### Start a GRASS GIS Session    
    # location/mapset: use random names for batch jobs
    string_length = 16
    location = binascii.hexlify(os.urandom(string_length)).decode('ascii')
    mapset   = 'PERMANENT'
    location_path = os.path.join(gisdb, location)
    
    ### Create new location
    #  from EPSG code:
    #  startcmd = [grass7bin,'-c', myepsg, '-e', location_path]
    #  from SHAPE or GeoTIFF file
    startcmd = [grass7bin, '-c', demraw_Path ,'-e', location_path]
      
    p = subprocess.Popen(startcmd, shell=True, 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    
    if p.returncode != 0:
        print('%s' %err.decode('ascii'))
        print('ERROR: %s' % err,file=sys.stderr)
        print('Cannot generate location')
        sys.exit(-1)
    else:
        print('%s' %err.decode('ascii'))
        print('Successifully created location') 
    
    ### Now we can use PyGRASS or GRASS Scripting library etc. after 
    # having started the session with gsetup.init() etc
    # import GRASS Python bindings (see also pygrass)
    import grass.script as grass
    import grass.script.setup as gsetup
    
    
    # launch session
    gsetup.init(gisbase,gisdb, location, mapset)
     
    # say hello
    grass.message('--- GRASS GIS 7.8: Current GRASS GIS 7 environment:')
    grassenv = grass.gisenv()
    print("--- Current GRASS GIS 7 environment: %s'" % grassenv)

    
    # print full proj parameter printing
    grass.message('--- GRASS GIS 7.8: Checking projection info:')
    in_proj = grass.read_command('g.proj', flags = 'jf')
    in_proj = in_proj.strip()
    grass.message("--- Found projection parameters: '%s'" % in_proj)
    print("--- Found projection parameters: '%s'" % in_proj)
    
    # set and show current region:
    grass.message('--- GRASS GIS 7.8: Checking computational region info:')
    demraw='demraw'
    grass.run_command('r.in.gdal', input=demraw_Path, output=demraw, flags='e',overwrite = True)
    grass.run_command('g.region',raster=demraw)
    in_region = grass.region()
    grass.message("--- Computational region: '%s'" % in_region)
    print("--- Computational region: '%s'" % in_region)
    
    # Perform watershed anaysis and output the results
    print("--- Perform r.watershed Analysis")   
    flwdir = 'flwdir'
    strseg = 'strseg'
    subbasin = 'subbasin'
    flwacc = 'flwacc'
    # grass.run_command('r.watershed', flags='4b', elevation=demraw, threshold=threshold, 
    #                   drainage=flwdir, basin=subbasin, stream=strseg, accumulation=flwacc, overwrite = True)
    grass.run_command('r.watershed', flags='4b', elevation=demraw, threshold=threshold, 
                      drainage=flwdir, basin=subbasin, stream=strseg, accumulation=flwacc, overwrite = True)

    print("--- Create stream vector")
    strsegThin = 'strseg_thin'
    strvector = "streams"
    grass.run_command('r.thin',input=strseg,output=strsegThin,overwrite = True)
    grass.run_command('r.to.vect', flags='v',input=strsegThin,  type='line',output=strvector, overwrite = True)

   
    print("--- Export results")       
    grass.run_command('r.out.gdal', flags='c', input=flwdir, format='GTiff', output=flwdir_Path, overwrite = True)
    grass.run_command('r.out.gdal', flags='c', input=subbasin, format='GTiff', output=subbasin_Path, overwrite = True)
    grass.run_command('r.out.gdal', flags='c', input=strseg, format='GTiff', output=strseg_Path, overwrite = True)
    grass.run_command('r.out.gdal', flags='c', input=flwacc, format='GTiff', output=flwacc_Path, overwrite = True)
    grass.run_command('v.out.ogr', input=strvector, format='ESRI_Shapefile', output=strvector_Path, 
                      type='line', flags='e2c', overwrite = True)
    
    # Finally remove the temporary batch location from disk
    print('Removing location %s' % location_path)
    shutil.rmtree(location_path)
    endT = time.time()
    print('time used: %s' %(endT-strT))
    
    return flwdir_Path, subbasin_Path, flwacc_Path, strseg_Path, strvector_Path
        
################################################   
### GRASS GIS r.stream.order  analysis
def GRASS_StreamOrder(strsegPath, flwdirPath, flwaccPath, demPath, out_dir):
    strT = time.time()
    # Set up GIS Environment
    grass7bin, gisdb, gisbase = Setup_GRASS()

    ### Paths to output files
    strOrder_Path = os.path.join(out_dir,'str_Strahler.tif')
    #strOrder_Path = os.path.join(out_dir,'str_Horton.tif')
    
    ### Start a GRASS GIS Session    
    # location/mapset: use random names for batch jobs
    string_length = 16
    location = binascii.hexlify(os.urandom(string_length)).decode('ascii')
    mapset   = 'PERMANENT'
    location_path = os.path.join(gisdb, location)
    
    ### Create new location
    #  from EPSG code:
    #  startcmd = [grass7bin,'-c', myepsg, '-e', location_path]
    #  from SHAPE or GeoTIFF file
    startcmd = [grass7bin, '-c', demPath ,'-e', location_path]
      
    p = subprocess.Popen(startcmd, shell=True, 
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    
    if p.returncode != 0:
        print('%s' %err.decode('ascii'))
        print('ERROR: %s' % err,file=sys.stderr)
        print('Cannot generate location')
        sys.exit(-1)
    else:
        print('%s' %err.decode('ascii'))
        print('Successifully created location') 
    
    ### Now we can use PyGRASS or GRASS Scripting library etc. after 
    # having started the session with gsetup.init() etc
    # import GRASS Python bindings (see also pygrass)
    import grass.script as grass
    import grass.script.setup as gsetup
    
    
    # launch session
    gsetup.init(gisbase,gisdb, location, mapset)
     
    # say hello
    grass.message('--- GRASS GIS 7.8: Current GRASS GIS 7 environment:')
    grassenv = grass.gisenv()
    print("--- Current GRASS GIS 7 environment: %s'" % grassenv)

    
    # print full proj parameter printing
    grass.message('--- GRASS GIS 7.8: Checking projection info:')
    in_proj = grass.read_command('g.proj', flags = 'jf')
    in_proj = in_proj.strip()
    grass.message("--- Found projection parameters: '%s'" % in_proj)
    print("--- Found projection parameters: '%s'" % in_proj)
    
    # load input raster
    demElv='dem'
    strseg = 'strseg'
    flwdir = 'flwdir'
    flwacc = 'flwacc'
    grass.run_command('r.in.gdal', input=demPath, output=demElv, flags='e',overwrite = True)        
    grass.run_command('r.in.gdal', input=strsegPath, output=strseg, flags='e',overwrite = True)
    grass.run_command('r.in.gdal', input=flwdirPath, output=flwdir, flags='e',overwrite = True)
    grass.run_command('r.in.gdal', input=flwaccPath, output=flwacc, flags='e',overwrite = True)
    
    
    # set and show current region:
    grass.message('--- GRASS GIS 7.8: Checking computational region info:') 
    grass.run_command('g.region',raster=demElv)
    in_region = grass.region()
    grass.message("--- Computational region: '%s'" % in_region)
    print("--- Computational region: '%s'" % in_region)
    
    # Perform stream order analysis anaysis and output the results
    print("--- Perform r.stream.order Analysis")
    outName = 'str_Strahler'
    grass.run_command('r.stream.order', flags='z', stream_rast=strseg, direction=flwdir,
                      accumulation=flwacc, elevation=demElv, strahler=outName, overwrite = True)
    # outName = 'str_Horton'
    # grass.run_command('r.stream.order', flags='z', stream_rast=strseg, direction=flwdir,
    #                   accumulation=flwacc, elevation=dem, horton=outName, overwrite = True)
    print("--- Export results")       
    grass.run_command('r.out.gdal', flags='c', input=outName, format='GTiff', output=strOrder_Path, overwrite = True)

    # Finally remove the temporary batch location from disk
    print('Removing location %s' % location_path)
    shutil.rmtree(location_path)
    endT = time.time()
    print('time used: %s' %(endT-strT))
    
    return strOrder_Path      