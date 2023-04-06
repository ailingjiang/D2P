# -*- coding: utf-8 -*-
"""
Modified from lidar python package (https://github.com/opengeos/lidar)
to better characterize depressions in terms of the sceening method

"""

import os
import math
import time
import shutil
import numpy as np
import richdem as rd
from scipy import ndimage
from skimage import measure
from osgeo import gdal, ogr, osr

class Depression:
    """The class for storing depression info.
    """      
    def __init__(self, id, level, count, size, volume, meanDepth, maxDepth, minElev, bndElev, spillElv, inNbrId, regionId,
                 perimeter, major_axis, minor_axis, elongatedness, eccentricity, orientation, area_bbox_ratio):
        self.id = id
        self.level = level
        self.count = count
        self.size = size
        self.volume = volume
        self.meanDepth = meanDepth
        self.maxDepth = maxDepth
        self.minElev = minElev
        self.bndElev = bndElev
        self.spillElv = spillElv
        self.inNbrId = inNbrId
        self.regionId = regionId
        self.perimeter = perimeter
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.elongatedness = elongatedness
        self.eccentricity = eccentricity
        self.orientation = orientation
        self.area_bbox_ratio = area_bbox_ratio


    
def get_min_max_nodata(image):
    """Gets the minimum, maximum, and no_data value of a numpy array.

    Args:
        image (np.array): The numpy array containing the image. 

    Returns:
        tuple: The minimum, maximum, and no_data value.
    """
    max_elev = np.max(image)
    nodata = pow(10, math.floor(math.log10(np.max(image))) + 2) - 1  # assign no data value
    image[image <= 0] = nodata  # change no data value
    min_elev = np.min(image)
    return min_elev, max_elev, nodata


# set input image parameters for level set method
def set_image_paras(no_data, min_size, min_size_sub, min_depth, resolution):
    """Sets the input image parameters for level-set method.

    Args:
        no_data (float): The no_data value of the input DEM.
        min_size (int): The minimum nuber of pixels to be considered as a depression.
        min_size_sub (int): The minimum nuber of pixels to be considered as a leaf within a compound depression.
        min_depth (float): The minimum depth to be considered as a depression.
        resolution (float): The spatial resolution of the DEM.

    Returns:
        dict: A dictionary containing image parameters.
    """    
    image_paras = {}
    image_paras["no_data"] = no_data
    image_paras["min_size"] = min_size
    image_paras["min_size_sub"] = min_size_sub
    image_paras["min_depth"] = min_depth
    image_paras["resolution"] = resolution
    return image_paras


def get_image_paras(image_paras):
    """Gets image parameters.

    Args:
        image_paras (dict): The dictionary containing image parameters.

    Returns:
        tuple: A tuple containing no_data, min_size, min_size_sub, min_depth, resolution.
    """    
    no_data = image_paras["no_data"]
    min_size = image_paras["min_size"]
    min_size_sub = image_paras["min_size_sub"]
    min_depth = image_paras["min_depth"]
    resolution = image_paras["resolution"]
    return no_data, min_size, min_size_sub, min_depth, resolution


def regionGroup(img_array, min_size, no_data):
    """IdentifIies regions based on region growing method

    Args:
        img_array (np.array): The numpy array containing the image.
        min_size (int): The minimum number of pixels to be considered as a depression.
        no_data (float): The no_data value of the image.

    Returns:
        tuple: The labelled objects and total number of labels. 
    """       
    img_array[img_array == no_data] = 0
    label_objects, nb_labels = ndimage.label(img_array)
    sizes = np.bincount(label_objects.ravel())
    mask_sizes = sizes > min_size
    mask_sizes[0] = 0
    image_cleaned = mask_sizes[label_objects]
    label_objects, nb_labels = ndimage.label(image_cleaned)
    # nb_labels is the total number of objects. 0 represents background object.
    return label_objects, nb_labels


def writeObject(img_array, obj_array, bbox):
    """Writes depression objects to the original image.

    Args:
        img_array (np.array): The output image array.
        obj_array (np.array): The numpy array containing depression objects.
        bbox (list): The bounding box of the depression object.

    Returns:
        np.array: The numpy array containing the depression objects.
    """    
    min_row, min_col, max_row, max_col = bbox
    roi = img_array[min_row:max_row, min_col:max_col]
    roi[obj_array > 0] = obj_array[obj_array > 0]
    return img_array


def writeRaster(arr, out_path, template):
    """Saves an numpy array as a GeoTIFF.

    Args:
        arr (np.array): The numpy array containing the image.
        out_path (str): The file path to the output GeoTIFF.
        template (str): The file path to the template image containing projection info. 

    Returns:
        np.array: The numpy array containing the image.
    """
    no_data = 0
    # First of all, gather some information from the template file
    data = gdal.Open(template)
    [cols, rows] = arr.shape
    trans = data.GetGeoTransform()
    proj = data.GetProjection()
    # nodatav = 0 #data.GetNoDataValue()
    # Create the file, using the information from the template file
    outdriver = gdal.GetDriverByName("GTiff")
    # http://www.gdal.org/gdal_8h.html
    # GDT_Byte = 1, GDT_UInt16 = 2, GDT_UInt32 = 4, GDT_Int32 = 5, GDT_Float32 = 6,
    outdata   = outdriver.Create(str(out_path), rows, cols, 1, gdal.GDT_UInt32)
    # Write the array to the file, which is the original array in this example
    outdata.GetRasterBand(1).WriteArray(arr)
    # Set a no data value if required
    outdata.GetRasterBand(1).SetNoDataValue(no_data)
    # Georeference the image
    outdata.SetGeoTransform(trans)
    # Write projection information
    outdata.SetProjection(proj)
    return arr


def polygonize(img, shp_path):
    """Converts a raster image to vector.

    Args:
        img (str): File path to the input image.
        shp_path (str): File path to the output shapefile.
    """      
    # mapping between gdal type and ogr field type
    type_mapping = {gdal.GDT_Byte: ogr.OFTInteger,
                    gdal.GDT_UInt16: ogr.OFTInteger,
                    gdal.GDT_Int16: ogr.OFTInteger,
                    gdal.GDT_UInt32: ogr.OFTInteger,
                    gdal.GDT_Int32: ogr.OFTInteger,
                    gdal.GDT_Float32: ogr.OFTReal,
                    gdal.GDT_Float64: ogr.OFTReal,
                    gdal.GDT_CInt16: ogr.OFTInteger,
                    gdal.GDT_CInt32: ogr.OFTInteger,
                    gdal.GDT_CFloat32: ogr.OFTReal,
                    gdal.GDT_CFloat64: ogr.OFTReal}

    ds = gdal.Open(img)
    prj = ds.GetProjection()
    srcband = ds.GetRasterBand(1)
    dst_layername = "Shape"
    drv = ogr.GetDriverByName("ESRI Shapefile")
    dst_ds = drv.CreateDataSource(shp_path)
    srs = osr.SpatialReference(wkt=prj)

    dst_layer = dst_ds.CreateLayer(dst_layername, srs=srs)
    raster_field = ogr.FieldDefn('id', type_mapping[srcband.DataType])
    dst_layer.CreateField(raster_field)
    gdal.Polygonize(srcband, srcband, dst_layer, 0, [], callback=None)
    del img, ds, srcband, dst_ds, dst_layer


def img_to_shp(in_img_dir, out_shp_dir):
    """Converts images in a selected folder to shapefiles

    Args:
        in_img_dir (str): The input iimage directory.
        out_shp_dir (str): The output shapefile directory. 
    """    
    img_files = os.listdir(in_img_dir)
    for img_file in img_files:
        if img_file.endswith(".tif"):
            img_filename = os.path.join(in_img_dir, img_file)
            shp_filename = os.path.join(out_shp_dir, img_file.replace("tif", "shp"))
            polygonize(img_filename, shp_filename)


# # parallel processing
# def task(region, out_image, no_data, min_size, min_depth, interval, resolution):

#     label_id = region.label
#     img = region.intensity_image
#     # img[img == 0] = no_data
#     bbox = region.bbox
#     # out_obj = identifyDepression(img,label_id,no_data,min_size,min_depth)
#     # writeObject(out_image,out_obj,bbox)
#     out_obj = levelSet(img, label_id, no_data, min_size, min_depth, interval, resolution)
#     writeObject(out_image, out_obj, bbox)


def levelSet(img, region_id, obj_uid, image_paras,regionOutlet):
    """Identifies nested depressions using level-set method.

    Args:
        img (np.array): The numpy array containing the image.
        region_id (int): The unique id of the region.
        obj_uid (int): The object id of the region.
        image_paras (dict): The dictionary containing image parameters.

    Returns:
        tuple: (level image, depression list)
    """
    # unzip input parameters from dict
    no_data, min_size, min_size_sub, min_depth, resolution = get_image_paras(image_paras)

    level_img = np.zeros(img.shape)     # init output level image
    SmoothElev_img = np.zeros(img.shape)     # init output level image
    # flood_img = np.zeros(img.shape)     # init output flood time image

    (rows, cols) = img.shape
    max_elev = np.max(img)
    img[img == 0] = no_data
    min_elev = np.min(img)
    
    
    levelArray = -np.sort(-np.unique(img))
    if rows == 1 or cols == 1: 
        levelArray = np.insert(levelArray,0,no_data,axis = 0)


    # print("Processing Region # {} ...".format(region_id))
    # print("=========================================================================== Region: {}".format(region_id))
    unique_id = obj_uid
    parent_ids = {}  # store current parent depressions
    nbr_ids = {}  # store the inner-neighbor ids of current parent depressions
    dep_list = []  # list for storing depressions

    #for elev in np.arange(max_elev, min_elev, interval):  
    ### Slicing operation using top-down approach
    for elev in levelArray:
        img[img >= elev] = 0  # set elevation higher than xy-plane to zero
        label_objects, nb_labels = regionGroup(img, 0, no_data)

        if nb_labels == 0:   # if slicing results in no objects, quit
            break

        # objects = measure.regionprops(label_objects, img, coordinates='xy')
        ### find the root depression
        objects = measure.regionprops(label_objects, img)
        for i, object in enumerate(objects):
            (row, col) = object.coords[0]  # get a boundary cell
            bbox = object.bbox

            if len(parent_ids) == 0:  # This is the first depression, maximum depression
                # print("This is the maximum depression extent.")
                cells = object.area
                spill_elev = regionOutlet
                mean_depth = (regionOutlet * cells - np.sum(object.intensity_image)) / cells  # depression mean depth
                if (mean_depth >= min_depth) & (cells > min_size_sub):
                    size = cells * pow(resolution, 2)  # depression size
                    max_depth = spill_elev- object.min_intensity  # depression max depth
                    volume = mean_depth * cells * pow(resolution, 2)  # depression volume
                    min_elev = object.min_intensity   # depression min elevation
                    max_elev = object.max_intensity     # depression max elevation
                    # print("size = {}, max depth = {:.2f}, mean depth = {:.2f}, volume = {:.2f}, spill elev = {:.2f}".format(
                    #     size, max_depth, mean_depth, volume, spill_elev))
                    unique_id += 1
                    level = 1
                    if rows == 1 or cols == 1: 
                        perimeter = cells * resolution
                        major_axis = cells * resolution
                        minor_axis = resolution
                    else:
                        perimeter = object.perimeter * resolution                    
                        major_axis = object.major_axis_length * resolution
                        minor_axis = object.minor_axis_length * resolution
                    if minor_axis == 0:
                        minor_axis = resolution
                    elongatedness = major_axis * 1.0 / minor_axis
                    eccentricity = object.eccentricity
                    orientation = object.orientation / np.pi * 180
                    area_bbox_ratio = object.extent
                    dep_list.append(Depression(unique_id,level,cells,size,volume,mean_depth,max_depth,min_elev,max_elev,spill_elev,[],
                                               region_id, perimeter, major_axis, minor_axis, elongatedness, eccentricity,
                                               orientation, area_bbox_ratio))
                    parent_ids[unique_id] = 0  # number of inner neighbors
                    nbr_ids[unique_id] = []   # ids of inner neighbors
                    tmp_img = np.zeros(object.image.shape)
                    tmp_img[object.image] = unique_id
                    writeObject(level_img, tmp_img, bbox)  # write the object to the final image
                else:
                    tmp_img = np.zeros(object.image.shape)
                    tmp_img[object.image] = 0
                    writeObject(level_img, tmp_img, bbox)
                    return level_img, dep_list, SmoothElev_img

            else:  # identify inner neighbors of parent depressions
                # print("current id: {}".format(parent_ids.keys()))
                # (row, col) = object.coords[0]
                parent_id = level_img[row,col]
                if parent_id in parent_ids:
                    parent_ids[parent_id] += 1
                    nbr_ids[parent_id].append(i)

        ### find the leaves depression
        isChild = {}
        for key in parent_ids.copy():  # check how many inner neighbors each upper level depression has
            isChild[key] = 0
            # print('key: '+str(key))
            if parent_ids[key] > 1:  # if the parent has two or more children
                # print("Object id: {} has split into {} objects".format(key, parent_ids[key]))
                new_parent_keys = nbr_ids[key].copy()
                for new_key in new_parent_keys:
                    object = objects[new_key]
                    cells = object.area
                    spill_elev = elev
                    if cells > min_size_sub:
                        size = cells * pow(resolution, 2)
                        max_depth = spill_elev - object.min_intensity
                        mean_depth = (elev * cells - np.sum(object.intensity_image)) / cells
                        volume = mean_depth * cells * pow(resolution, 2)
                        min_elev = object.min_intensity
                        max_elev = object.max_intensity
                        # print("  --  size = {}, max depth = {:.2f}, mean depth = {:.2f}, volume = {:.2f}, spill elev = {:.2f}".format(
                        #         size, max_depth, mean_depth, volume, spill_elev))
                        unique_id += 1
                        level = 1
                        if rows == 1 or cols == 1: 
                            perimeter = cells * resolution
                            major_axis = cells * resolution
                            minor_axis = resolution
                        else:
                            perimeter = object.perimeter * resolution                    
                            major_axis = object.major_axis_length * resolution
                            minor_axis = object.minor_axis_length * resolution
                        
                        if minor_axis == 0:
                            minor_axis = resolution
                        elongatedness = major_axis * 1.0 / minor_axis
                        eccentricity = object.eccentricity
                        orientation = object.orientation / np.pi * 180
                        area_bbox_ratio = object.extent
                        dep_list.append(
                            Depression(unique_id, level, cells, size, volume, mean_depth, max_depth, min_elev, max_elev, spill_elev, [],
                                       region_id, perimeter, major_axis, minor_axis, elongatedness, eccentricity,
                                       orientation, area_bbox_ratio))
                        dep_list[key-1-obj_uid].inNbrId.append(unique_id)
                        parent_ids[unique_id] = 0
                        nbr_ids[unique_id] = []
                        bbox = object.bbox
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = unique_id
                        writeObject(level_img, tmp_img, bbox)
                        isChild[key] = isChild[key]+1
                    else:
                        nbr_ids[key].remove(new_key)
                        bbox = object.bbox
                        tmp_img = np.zeros(object.image.shape)
                        tmp_img[object.image] = spill_elev
                        writeObject(SmoothElev_img, tmp_img, bbox) 
                        writeObject(img, tmp_img, bbox)  

                if isChild[key] == 1:                    
                    parent_ids[key] = 0
                    nbr_ids[key] = []
                    dep_list.pop()
                    dep_list[key-1-obj_uid].inNbrId = []
                    parent_ids.pop(unique_id)
                    nbr_ids.pop(unique_id)
                    # print(dep_list)
                    # print(unique_id)
                    # print(key)                 
                    level_img[level_img==unique_id] = key
                    unique_id = unique_id-1                   
                elif key in parent_ids.keys():    # remove parent id that has split
                    parent_ids.pop(key)
                   #print('Afterall, pop from parent_ids',parent_ids)
            else:
                parent_ids[key] = 0     # if a parent depression has not split, keep it
                nbr_ids[key] = []
                # print('unspplit--key:'+str(key))
                # print('unsplot parent ' + str(parent_ids))
                # print('unsplit nbr_ids ' + str(nbr_ids))

    #     print(obj_uid)
    dep_list = updateLevel(dep_list, obj_uid)   # update the inner neighbors of each depression
    # for dep in dep_list:
    #     print("id: {} is level {}".format(dep.id, dep.level))
    del img

    return level_img, dep_list, SmoothElev_img


def updateLevel(dep_list, obj_uid):
    """Updates the inner neighbors of each depression.

    Args:
        dep_list (list): A list containing depression info.
        obj_uid (int): The unique id of an object.

    Returns:
        list: A list containing depression info.
    """    
    for dep in reversed(dep_list):
        if len(dep.inNbrId) == 0:
            dep.level = 1
        else:
            max_children_level = 0
            for id in dep.inNbrId:
               # print('idx: '+str(id-1-obj_uid))
                if dep_list[id-1-obj_uid].level > max_children_level:
                    max_children_level = dep_list[id-1-obj_uid].level
            dep.level = max_children_level + 1
    return dep_list


def obj_to_level(obj_img, dep_list):
    """Derives depression level image based on the depression id image and depression list.

    Args:
        obj_img (np.array): The numpy array containing the object image.
        dep_list (list): A list containing depression info.

    Returns:
        np.array: The numpy array containing the object level image.
    """    
    level_img = np.copy(obj_img)

    max_id = int(np.max(level_img))
    # print("max id = " + str(max_id))
    if max_id > 0:
        min_id = int(np.min(level_img[np.nonzero(level_img)]))
        # print("min_id = " + str(min_id))
        for i in range(min_id, max_id+1):
            level_img[level_img == i] = dep_list[i-1].level + max_id
    level_img = level_img - max_id

    return level_img


def write_dep_csv(dep_list, csv_file):
    """Saves the depression list to a CSV file.


    Args:
        dep_list (list): A list containing depression info. 
        csv_file (str): File path to the output CSV file.
    """    
    csv = open(csv_file, "w")
    header = "id" +","+"level"+","+"count"+","+"area"+","+"volume"+","+"avg-depth"+","+"max-depth"+","+\
             "min-elev"+","+"max-elev"+","+"spill-elev"+","+"children-id"+","+"region-id" + "," + "perimeter" + "," + "major-axis" + \
             "," + "minor-axis" + "," + "elongatedness" + "," + "eccentricity" + "," + "orientation" + "," + \
             "area-bbox-ratio"
    csv.write(header + "\n")
    for dep in dep_list:
        # id, level, size, volume, meanDepth, maxDepth, minElev, bndElev, inNbrId, nbrId = 0
        line = "{},{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{},{},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f},{:.3f}," \
               "{:.3f}".format(dep.id, dep.level, dep.count, dep.size, dep.volume, dep.meanDepth, dep.maxDepth,
                               dep.minElev,dep.bndElev, dep.spillElv, str(dep.inNbrId).replace(",",":"), dep.regionId, dep.perimeter,
                               dep.major_axis, dep.minor_axis, dep.elongatedness, dep.eccentricity, dep.orientation,
                               dep.area_bbox_ratio)
        csv.write(line + "\n")
    csv.close()


def extract_levels(level_img, obj_img, min_size, no_data, out_img_dir, out_shp_dir, template, bool_comb=False):
    """Extracts individual level image.

    Args:
        level_img (np.array): The numpy array containing the level image.
        obj_img (np.array): The numpy array containing the object image.
        min_size (int): The minimum number of pixels to be considered as a depression.
        no_data (float): The no_data value of the image.
        out_img_dir (str): The output image directory.
        out_shp_dir (str): The output shapefile directory.
        template (str): The file path to the template image.
        bool_comb (bool, optional): Whether to extract combined level image. Defaults to False.

    Returns:
        tuple: The single level image, properties of region grouped level image, properties of region grouped object image.
    """    
    max_level = int(np.max(level_img))
    combined_images = []
    single_images = []
    img = np.copy(level_img)

    digits = int(math.log10(max_level)) + 1  # determine the level number of output file name
    for i in range(1, max_level + 1):
        img[(img > 0) & (img <= i) ] = i
        tmp_img = np.copy(img)
        tmp_img[tmp_img > i] = 0
        if bool_comb == True:  # whether to extract combined level image
            combined_images.append(np.copy(tmp_img))
            filename_combined = "Combined_level_" + str(i).zfill(digits) + ".tif"
            out_file = os.path.join(out_shp_dir, filename_combined)
            writeRaster(tmp_img,out_file,template)

        lbl_objects, n_labels = regionGroup(tmp_img, 0, no_data)
        # regs = measure.regionprops(lbl_objects, level_img, coordinates='xy')
        regs = measure.regionprops(lbl_objects, level_img)
        # regs2 = measure.regionprops(lbl_objects, obj_img, coordinates='xy')
        regs2 = measure.regionprops(lbl_objects, obj_img)

        sin_img = np.zeros(img.shape)

        for index, reg in enumerate(regs):
            uid = regs2[index].min_intensity
            if reg.max_intensity >= i:
                bbox = reg.bbox
                tmp_img = np.zeros(reg.image.shape)
                tmp_img[reg.image] = uid
                writeObject(sin_img, tmp_img, bbox)

        # for reg in regs:
        #     if reg.max_intensity >= i:
        #         bbox = reg.bbox
        #         tmp_img = np.zeros(reg.image.shape)
        #         tmp_img[reg.image] = i
        #         writeObject(sin_img, tmp_img, bbox)
        del tmp_img
        # single_images.append(np.copy(sin_img))
        filename_single = "Single_level_" + str(i).zfill(digits) + ".shp"
        out_shp_file = os.path.join(out_shp_dir, filename_single)

        out_img_file = os.path.join(out_img_dir, "tmp.tif")
        writeRaster(sin_img, out_img_file, template)
        polygonize(out_img_file, out_shp_file)
        # writeRaster(sin_img,out_file,template)
        del sin_img, regs, regs2

    del img
    return True


def getMetadata(img):
    """Gets rdarray metadata.

    Args:
        img (rdarray): The richDEM array containing the image.

    Returns:
        tuple: no_data, projection, geotransform, cell_size
    """    
    no_data = img.no_data
    projection = img.projection
    geotransform = img.geotransform
    cell_size = np.round(geotransform[1], decimals=2)
    return no_data, projection, geotransform, cell_size


def np2rdarray(in_array, no_data, projection, geotransform):
    """Converts numpy array to rdarray.

    Args:
        in_array (np.array): The input numpy array containing the image.
        no_data (float): The no_data value of the image.
        projection (str): The projection coordinate system of the image.
        geotransform (str): The geotransform of the image.

    Returns:
        rdarray: The richDEM array containing the image.
    """    
    out_array = rd.rdarray(in_array, no_data=no_data)
    out_array.projection = projection
    out_array.geotransform = geotransform
    return out_array


def DelineateDepressions_Mod(in_fill, in_orig, min_size, min_size_sub, min_depth, out_dir, bool_level_shp=False):
    """Delineates nested depressions.

    Args:
        in_sink (str): The file path to the sink image.
        min_size (int): The minimum number of pixels to be considered as a depression.
        min_size_sub (int): The minimum nuber of pixels to be considered as a leaf within a compound depression.
        min_depth (float): The minimum depth to be considered as a depression.
        out_dir (str): The file path to the output directory.
        bool_level_shp (bool, optional): Whether to generate shapefiles for each individual level. Defaults to False.

    Returns:
        tuple: The output level image, and the output object image.
    """
    # The following parameters can be used by default
    out_img_dir = os.path.join(out_dir, "img-level")
    out_shp_dir = os.path.join(out_dir, "shp-level")
    out_obj_file = os.path.join(out_dir, "depression_id.tif")
    out_level_file = os.path.join(out_dir, "depression_level.tif")
    out_vec_file = os.path.join(out_dir, "depressions.shp")
    out_csv_file = os.path.join(out_dir, "depressions_info.csv")
    smoothElev_obj_file = os.path.join(out_dir, "smoothed_elev_within_compound_depressions.tif")

    init_time = time.time()

    # delete contents in output folder if existing
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    if os.path.exists(out_img_dir):
        shutil.rmtree(out_img_dir)
    os.mkdir(out_img_dir)
    if os.path.exists(out_shp_dir):
        shutil.rmtree(out_shp_dir)
    os.mkdir(out_shp_dir)

    #print("Reading filled and original data ...")
    read_time = time.time()
    dem_fill = rd.LoadGDAL(in_fill)
    image= rd.LoadGDAL(in_orig)
    #curvature = rd.TerrainAttribute(image, attrib='curvature')
    dem_diff =  dem_fill-image
    image[dem_diff == 0] = 0
    #curvature[dem_diff == 0] = 0
    

    #image = rd.LoadGDAL(in_sink)
    no_data_raw, projection, geotransform, resolution = getMetadata(image)
    rows_cols = image.shape
    #print("rows, cols: " + str(rows_cols))
    #print("Pixel resolution: " + str(resolution))
    print("Read data time: {:.4f} seconds".format(time.time() - read_time))

    min_elev, max_elev, no_data = get_min_max_nodata(image)  # set nodata value to a large value, e.g., 9999
    # initialize output image
    obj_image = np.zeros(image.shape)  # output depression image with unique id for each nested depression
    level_image = np.zeros(image.shape)  # output depression level image
    smoothElev_image = np.zeros(image.shape) # output smooth elevation within compound depression

    # nb_labels is the total number of objects. 0 represents background object.
    label_objects, nb_labels = regionGroup(image, min_size, no_data)
    
    # regions = measure.regionprops(label_objects, image, coordinates='xy')
    regions = measure.regionprops(label_objects, image)
    del image  # delete the original image to save memory
    prep_time = time.time()
    #print("Data preparation time: {:.4f} seconds".format(prep_time - init_time))
    #print("Total number of regions: {}".format(nb_labels))

    identify_time = time.time()

    obj_uid = 0
    global_dep_list = []

    # loop through regions and identify nested depressions in each region using level-set method
    for region in regions:  # iterate through each depression region
        region_id = region.label
        img = region.intensity_image  # dem subset for each region
       # print(img.shape)
        bbox = region.bbox
        regionOutlet = dem_fill[region.coords[0][0],region.coords[0][1]]
        #curv_img = curvature[bbox[0]:bbox[2],bbox[1]:bbox[3]]
        #print(curv_sub.shape)
        # save all input parameters needed for level set methods as a dict
        image_paras = set_image_paras(no_data, min_size, min_size_sub, min_depth, resolution)

        # execute level set methods
        out_obj, dep_list, smoothElev_obj = levelSet(img, region_id, obj_uid, image_paras, regionOutlet)

        for dep in dep_list:
            global_dep_list.append(dep)

        obj_uid += len(dep_list)

        level_obj = obj_to_level(out_obj, global_dep_list)
        obj_image = writeObject(obj_image, out_obj, bbox)       # write region to whole image        
        level_image = writeObject(level_image, level_obj, bbox)
        smoothElev_image =  writeObject(smoothElev_image, smoothElev_obj, bbox)
        


        del out_obj, level_obj, region, smoothElev_obj 

    del regions, label_objects

    print("=========== Run time statistics =========== ")
    print("(rows, cols):\t\t\t {0}".format(str(rows_cols)))
    print("Pixel resolution:\t\t {0} m".format(str(resolution)))
    print("Number of regions:\t\t {0}".format(str(nb_labels)))
    print("Data preparation time:\t\t {:.4f} s".format(prep_time - init_time))
    print("Identify level time:\t\t {:.4f} s".format(time.time() - identify_time))

    write_time = time.time()
    # writeRaster(obj_image, out_obj_file, in_sink)
    # writeRaster(level_image, out_level_file, in_sink)
    # SaveGDAL function can only save data as floating point
    level_image = np2rdarray(np.int32(level_image), no_data_raw, projection, geotransform)
    rd.SaveGDAL(out_level_file, level_image)
    obj_image = np2rdarray(np.int32(obj_image), no_data_raw, projection, geotransform)
    rd.SaveGDAL(out_obj_file, obj_image)
    smoothElev_image = np2rdarray(np.float32(smoothElev_image), no_data_raw, projection, geotransform)
    rd.SaveGDAL(smoothElev_obj_file, smoothElev_image)
    print("Write image time:\t\t {:.4f} s".format(time.time() - write_time))

    # converting object image to polygon
    level_time = time.time()
    polygonize(out_obj_file, out_vec_file)
    write_dep_csv(global_dep_list, out_csv_file)
    print("Polygonize time:\t\t {:.4f} s".format(time.time() - level_time))

    # extracting polygons for each individual level
    if bool_level_shp:
        level_time = time.time()
        extract_levels(level_image, obj_image, min_size, no_data, out_img_dir, out_shp_dir, in_orig, False)
        print("Extract level time:\t\t {:.4f} s".format(time.time() - level_time))
        shutil.rmtree(out_img_dir)
    else:
        shutil.rmtree(out_shp_dir)
        shutil.rmtree(out_img_dir)
    del level_image
    del obj_image

    end_time = time.time()
    print("Total run time:\t\t\t {:.4f} s".format(end_time - init_time))
    return out_obj_file, out_level_file, smoothElev_obj_file