
#### CHECK ####
#arcpy.AddMessage(f"Data type of {out_feature_class} is {type(out_feature_class)}")
#fields = arcpy.ListFields(out_feature_class) ### REMOVE
#varlist = [field.name for field in fields[1:]] ### REMOVE
#arcpy.AddMessage(f"Variable list: {varlist} FOR {out_feature_class}") ### REMOVE
#sum_p_NISP = 0
#with arcpy.da.SearchCursor(out_feature_class, ['p_NISP']) as cursor:
#    for row in cursor:
#        sum_p_NISP += row[0]
#arcpy.AddMessage(f"The sum of p_NISP is: {sum_p_NISP}")    
#rc = arcpy.management.GetCount(out_feature_class) ### REMOVE
#arcpy.AddMessage(f"Row count: {rc}") ### REMOVE
#arcpy.AddMessage("Breakpoint reached") ### REMOVE
#pdb.set_trace() ### REMOVE



        #rc = arcpy.management.GetCount(hadash_hexes) ### REMOVE
        #arcpy.AddMessage(f"Row count: {rc}") ### REMOVE
        #arcpy.AddMessage("Breakpoint reached") ### REMOVE
        #pdb.set_trace() ### REMOVE   


#### END CHECK ####



import arcpy
import os
import numpy as np
import pdb
import locale
import json
import time
import gc



def getColorScheme(number_of_breaks, color_pallette):
    # color pallette generator
    try:
        first_break_rgba = [255, 255, 255, 0]
        second_break_rgba = [255, 255, 255, 100]
        last_colors = {
            "green": [0, 100, 0, 100],
            "blue": [0, 0, 139, 100],
            "black": [0, 0, 0, 100]
        }

        last_break_rgba = last_colors.get(color_pallette, last_colors["black"])

        gradient_colors = []
        for i in range(1, number_of_breaks-2):
            alpha = i / (number_of_breaks) 
            color = np.array(second_break_rgba) * (1 - alpha) + np.array(last_break_rgba) * alpha
            gradient_colors.append(list(color.astype(int)))

        color_scheme = [first_break_rgba] + [second_break_rgba] + gradient_colors + [last_break_rgba]
        return color_scheme

    except Exception as e:
        arcpy.AddError(f"Error during validation: {str(e)}")
        raise

def getMinMax(input, list_of_variables):
    # gets the values range,- KEEP NOTICE, that it gets the range out of varlist !!!
    try:
        overall_max = float('-inf')
        overall_min = float('inf')

        with arcpy.da.SearchCursor(input, list_of_variables) as cursor:
            for row in cursor:
                for value in row:
                    if value is not None:
                        if value > overall_max:
                            overall_max = value
                        if value < overall_min:
                            overall_min = value

        arcpy.AddMessage(f"Minimal output value: {overall_min:.2e}\nMaximal output value: {overall_max:.2e} ")
        return overall_min, overall_max

    except Exception as e:
        arcpy.AddError(f"Error during validation: {str(e)}")
        return False

def varlist_generator(layer, mode):
    try:
        # Define technical and excluded variables
        technical_vars = ["OBJECTID", "GRID_ID", "Shape", "Shape_Length", "Shape_Area", "CU", "Unit", "Collection_Unit", "Number_of_surveyors", "Surface_Visibility", "Elevation", "Slope"]
        excluded_vars = ["date_f_LP", "date_p_LB_I"]
        
        # Get all fields from the layer
        whole_list_of_fields = arcpy.ListFields(layer)

        # Filter fields based on exclusions
        fields = [field.name for field in whole_list_of_fields if field.name not in technical_vars and field.name not in excluded_vars]

        var_sets = {
            'condition': ['f_Total_Patinated', 'f_Total_Fresh', 'f_Total_Burnt', 'f_Total_not_burnt', 'f_Broken', 'f_Complete'],
            'dating': [f for f in fields if f.startswith("date")],
            'general': ["p_NISP", "f_NISP_Total"],
            'types': [f for f in fields if f.startswith("f_MNI")],
            'specific': ["date_f_UP_EPI", "date_f_NEO_CHAL", "f_MNI_Total"]
        }

        varlist = var_sets.get(mode, fields)
        arcpy.AddMessage(f'List of variables used: {varlist}')
        return varlist
    
    except Exception as e:
        arcpy.AddError(f"Error during varlist generation: {str(e)}")
        raise

def calculateClassBreaksAndColors(project_name, map_name, feature_name, list_of_variables, breaks_count, decimals_number, color_pallette):
    try:
        aprx = arcpy.mp.ArcGISProject(project_name)
        mp = aprx.listMaps(map_name)[0]
        lyr = mp.listLayers(feature_name)[0]
        
        min_val, max_val = getMinMax(lyr, list_of_variables)
        classBreakValues = [round(min_val + (max_val - min_val) * i / (breaks_count - 1), decimals_number) for i in range(breaks_count)]
        classBreakLabels = ['\u2264 ' + str(value) for value in classBreakValues]
        color_scheme = getColorScheme(len(classBreakValues), color_pallette)

        return classBreakValues, classBreakLabels, color_scheme
        
    except Exception as e:
        arcpy.AddError(f"Error during varlist generation: {str(e)}")
        raise

def setSymbology(project_name, map_name, feature_name, classBreakValues, classBreakLabels, color_scheme, field_name, mode, color_pallette):
    try:
        # Load the necessary layer
        aprx = arcpy.mp.ArcGISProject(project_name)
        mp = aprx.listMaps(map_name)[0]
        lyr = mp.listLayers(feature_name)[0]

        # Update the symbology of the layer
        sym = lyr.symbology
        
        sym.updateRenderer('GraduatedColorsRenderer')
        
        if mode == "natural breaks":
            sym.renderer.classificationMethod = "NaturalBreaks"
            classBreakValues = [break_.upperBound for break_ in sym.renderer.classBreaks]
            classBreakLabels = ['\u2264 ' + str(round(value, 1)) for value in classBreakValues]
        elif mode == "standard deviation":
            sym.renderer.classificationMethod = "StandardDeviation"
            num_breaks = sym.renderer.breakCount
            color_scheme = getColorScheme(num_breaks, color_pallette)
            classBreakLabels = ['\u2264 ' + str(round(brk.upperBound, 1)) + ' (' + brk.label.replace("Std. Dev.", "SD") + ')' for brk in sym.renderer.classBreaks]
        else:
            sym.renderer.classificationMethod = "EqualInterval"
        
        
        sym.renderer.classificationField = field_name
        
        if mode != "standard deviation":
            sym.renderer.breakCount = len(classBreakValues)

        # Update the breaks
        
        for count, brk in enumerate(sym.renderer.classBreaks):
            if mode == "min to max":
                brk.upperBound = classBreakValues[count]
            brk.symbol.color = {'RGB': color_scheme[count]}
            if count == 0:
                brk.symbol.outlineColor = {'RGB': [155, 155, 155, 0]}
            else:
                brk.symbol.outlineColor = {'RGB': [155, 155, 155, 50]}
        
        # Save layer's symbology
        lyr.symbology = sym
        
        # Save the JSON style file and removing the layer from the map
        output_path = os.path.dirname(arcpy.env.workspace)
        tmp_lyrx = f"{output_path}/temp.lyrx"
        lyr.saveACopy(tmp_lyrx)
        mp.removeLayer(lyr)

        # Load the JSON file
        with open(tmp_lyrx, 'r') as file:
            data = json.load(file)

        # Locate the renderer and update the breaks' labels
        renderer = data["layerDefinitions"][0]["renderer"]
        if "breaks" in renderer and len(renderer["breaks"]) == len(classBreakLabels):
            for i, break_item in enumerate(renderer["breaks"]):
                break_item["label"] = classBreakLabels[i]

        # Save the updated JSON file
        updated_file_path = tmp_lyrx
        with open(updated_file_path, 'w') as updated_file:
            json.dump(data, updated_file, indent=4)

        # Add the layer back with edited JSON and remove the lyrx file
        mp.addLayer(arcpy.mp.LayerFile(tmp_lyrx))
        os.remove(tmp_lyrx)

    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages())


def generate_outline(project_name, main_map, feature_name, bufsize):
    """
    Description:
        Generates an 'Area of Investigation' polygon for the given input features using a buffer approach.
    Args:
        input (str): Path to the input feature layer.
        bufsize (float): Buffer distance in meters.
    Returns:
        None
    """
    try:
        arcpy.AddWarning("You haven't specified an 'Area of Investigation'.")
        arcpy.AddMessage(f"Generating 'Area of Investigation' for '{feature_name}' with smoothing factor {bufsize}")
        
        aprx = arcpy.mp.ArcGISProject(project_name)
        mp = aprx.listMaps(main_map)[0]
        lyr = mp.listLayers(feature_name)[0]
        output = "research_area"

        # Buffer the input features to create a buffered polygon
        buffered_polygon = "in_memory/buffered_polygon"
        arcpy.analysis.Buffer(
            in_features = lyr,
            out_feature_class = buffered_polygon,
            buffer_distance_or_field = f"{bufsize} Meters",
            line_side = "FULL",
            line_end_type = "ROUND",
            dissolve_option = "ALL",
            dissolve_field = None,
            method = "GEODESIC"
        )

        # Debuffer the polygon back
        debuffered_polygon = "in_memory/debuffered_polygon"
        arcpy.analysis.Buffer(
            in_features = buffered_polygon,
            out_feature_class = debuffered_polygon,
            buffer_distance_or_field = f"-{bufsize-bufsize/10} Meters",
            line_side = "FULL",
            line_end_type = "ROUND",
            dissolve_option = "NONE",
            dissolve_field = None,
            method = "GEODESIC"
        )

        # Smooth the polygon outline
        smoothed_polygon = "in_memory/smoothed_polygon"
        arcpy.cartography.SmoothPolygon(
            in_features = debuffered_polygon,
            out_feature_class = smoothed_polygon,
            algorithm = "PAEK",
            tolerance = f"{bufsize*3} Meters",
            endpoint_option = "FIXED_ENDPOINT",
            error_option = "NO_CHECK",
            in_barriers = None
        )
                
        # Save the result to the output
        arcpy.management.CopyFeatures(smoothed_polygon, output)
        arcpy.AddMessage(f"Area of Investigation '{output}' has been generated.")

        mp.addDataFromPath(os.path.join(arcpy.env.workspace, output))
        

        res_area = mp.listLayers(output)[0]
        sym = res_area.symbology
    
        sym.renderer.symbol.color = {'RGB': [255, 255, 190, 25]}
        sym.renderer.symbol.outlineColor = {'RGB': [204, 204, 204, 100]}
        sym.renderer.symbol.outlineWidth = 1

        res_area.symbology = sym



    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise


def heatmapSymbology(project_name, map_name, feature_name, points_layer_name, field_name, cell_size, search_dist, number_of_breaks, barrier_fc=None):
    try:
        arcpy.env.overwriteOutput = True

        # Initiate the necessary layers
        aprx = arcpy.mp.ArcGISProject(project_name)
        mp = aprx.listMaps(map_name)[0]
        lyr_points = mp.listLayers(points_layer_name)[0]
        prefix = "HeatmapKD"
        heatmap_layer_name = f"{prefix}_{field_name}"
        heatmap_raster_layer_path = os.path.join(arcpy.env.workspace, heatmap_layer_name)
        layers_to_remove = ("Kernel", "Hotspot", "Heatmap")
        excluded_layers = "research_area"
        rasters = arcpy.ListRasters("HeatmapKD*")

        for layer in mp.listLayers():
            if layer.name.startswith(layers_to_remove):
                arcpy.AddMessage(f"deleting layers: {layer}")
                mp.removeLayer(layer)
            if layer.isFeatureLayer and layer.name not in excluded_layers:
                layer.visible = False

        arcpy.management.ClearWorkspaceCache(arcpy.env.workspace)
        aprx.save()

        if rasters:
            for raster in rasters:
                can_lock = arcpy.TestSchemaLock(raster)
                if can_lock:
                    arcpy.AddMessage("Raster is not locked. Schema lock acquired successfully.")
                else:
                    arcpy.AddMessage("Raster is locked or cannot acquire schema lock.")
                arcpy.management.Delete(raster)
                arcpy.AddMessage(f"Deleting raster: {raster}")
                arcpy.management.Delete(raster)
                arcpy.ClearEnvironment("workspace")
            arcpy.AddMessage(f"Deleted {len(rasters)} raster(s) starting with '{prefix}'.")
                


        if arcpy.Exists(heatmap_raster_layer_path):
            arcpy.AddMessage(f"Attempting to delete existing raster: {heatmap_raster_layer_path}")
            arcpy.Delete_management(heatmap_raster_layer_path)
            arcpy.AddMessage(f"Deleted existing raster: {heatmap_raster_layer_path}")
            aprx.save()

        kernel_raster = None
        if kernel_raster is not None and arcpy.Exists(kernel_raster):
            arcpy.AddMessage("it exists!")
        else:
            arcpy.AddMessage("it DOESNT EXIST!!")
        # Creating a Kernel Density raster
        kernel_raster = arcpy.sa.KernelDensity(
            in_features = lyr_points,
            population_field = field_name,
            cell_size = cell_size,
            search_radius = search_dist,
            area_unit_scale_factor = "SQUARE_METERS",
            out_cell_values = "EXPECTED_COUNTS",
            method = "GEODESIC",
            in_barriers = barrier_fc                                      # REDO
        )

        kernel_raster.save(heatmap_raster_layer_path)
        mp.addDataFromPath(heatmap_raster_layer_path)
        arcpy.AddMessage(f"Heatmap layer '{heatmap_layer_name}' added successfully.")
        lyr = mp.listLayers(heatmap_layer_name)[0]

        # Move the heatmap layer to the top (before the first layer)
     #   if lyr:
    #        top_layer = mp.listLayers()[0]  # Get the topmost layer
    #        mp.moveLayer(top_layer, lyr, "BEFORE")

#        arcpy.AddMessage(f"Kernel has been added to the map.")
     #   aprx.save()   
        sym = lyr.symbology

      #  sym.renderer.classificationMethod = "EqualInterval"
        sym.colorizer.breakCount = number_of_breaks
        #class_breaks = sym.colorizer.classBreaks
        classBreakLabels = ['\u2264 ' + str(round(brk.upperBound, 1)) for brk in sym.colorizer.classBreaks]
        sym.colorizer.classBreaks[0].color = {'RGB': [0,0,0,0]}

        for brk, label in zip(sym.colorizer.classBreaks, classBreakLabels):
            brk.label = label

        lyr.symbology = sym
#        
#        arcpy.AddMessage("Breakpoint reached") ### REMOVE
#        pdb.set_trace() ### REMOVE
#        arcpy.AddMessage(f"labels: {classBreakLabels}")



        # Save the JSON style file and removing the layer from the map
        output_path = os.path.dirname(arcpy.env.workspace)
        tmp_lyrx = os.path.join(output_path, "temp.lyrx")
        #arcpy.management.SaveToLayerFile(in_layer = lyr, out_layer = tmp_lyrx)  #is_relative_path="RELATIVE"
        
        lyr.saveACopy(tmp_lyrx)
  #      mp.removeLayer(lyr)

        # Load the JSON file
        with open(tmp_lyrx, 'r') as file:
            data = json.load(file)
#       
#        if "layerDefinitions" in data and data["layerDefinitions"]:
#            for layer_def in data["layerDefinitions"]:
#                if "name" in layer_def:
#                    layer_def["name"] += "_edited"
                    
        arcpy.AddMessage(f"Original class break labels: {classBreakLabels}")

        # Locate the renderer and update the breaks' labels
        colorizer = data["layerDefinitions"][0]["colorizer"]
        if "classBreaks" in colorizer and len(colorizer["classBreaks"]) == len(classBreakLabels):
           for i, break_item in enumerate(colorizer["classBreaks"]):
               break_item["label"] = classBreakLabels[i]




        # Save the updated JSON file
     #   updated_file_path = tmp_lyrx
        with open(tmp_lyrx, 'w') as updated_file:
            json.dump(data, updated_file, indent=4)
     #   arcpy.AddMessage("Breakpoint reached") ### REMOVE
 #       pdb.set_trace() ### REMOVE
 #       mp.removeLayer(lyr)
        
  #      out_layer = arcpy.management.ApplySymbologyFromLayer(in_layer=heatmap_layer_name, in_symbology_layer=tmp_lyrx)
  #      mp.addLayer(arcpy.mp.LayerFile(out_layer))    
  #  aprx.save()
#        mp.removeLayer(lyr)



 #       arcpy.AddMessage("Breakpoint reached") ### REMOVE
#        pdb.set_trace() ### REMOVE


       # updated_layer = arcpy.mp.LayerFile(tmp_lyrx)
      #  for l in updated_layer.listLayers():
      #      l.name = "Kernel Density"  # Set the desired name
   #     mp.addLayer(updated_layer)
    #    lyr.updateLayer(updated_layer)
   #     lyr = mp.listLayers(heatmap_layer_name)[0]
  #      arcpy.AddMessage(f"LAYER: {tmp_lyrx}")
  #      arcpy.AddMessage(f"updated_layer.name: {updated_layer}")
#       arcpy.AddMessage(f"heatmap_layer_name: {heatmap_layer_name}")

   #     lyr.symbology = tmp_lyrx

 #       mp.addLayer(updated_layer) 
    #    arcpy.SetParameterSymbology(in_layer=heatmap_layer_name, in_symbology_layer=tmp_lyrx)
     #   reclass = arcpy.sa.Reclassify(
     #       heatmap_layer_name,
     #       "Value",

      #  lyr = arcpy.ApplySymbologyFromLayer_management(str(heatmap_layer_name), tmp_lyrx, update_symbology="MAINTAIN")[0]
   #     lyr = arcpy.management.ApplySymbologyFromLayer(in_layer=str(heatmap_layer_name), in_symbology_layer=tmp_lyrx, update_symbology="MAINTAIN")[0]

    #    aprx.save()
     #   lyr_two = arcpy.mp.LayerFile(updated_file_path)
 #       for_adding = arcpy.management.MakeRasterLayer(lyr_two, "TEST RASTER LAYER")[0]
        # Add the layer back with edited JSON and remove the lyrx file
   #     mp.addLayer(arcpy.mp.LayerFile(tmp_lyrx)) # THIS WORKS BUT CORRUPTS THE RASTER
   #     os.remove(tmp_lyrx)

    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages())
    except Exception as e:
        arcpy.AddError(f"Unexpected error: {e}")
    finally:
        # Ensure cleanup of any in-memory datasets
        arcpy.management.Delete("in_memory")
        arcpy.ClearEnvironment("workspace")

def main():
    try:
        gc.collect()
        # Parameters
        symbology_mode = arcpy.GetParameterAsText(0)
        layout_name = arcpy.GetParameterAsText(1)
        polygons_layer_name = arcpy.GetParameterAsText(2) # RENAME
        points_layer_name = f"{polygons_layer_name}_points"
        breaks = int(arcpy.GetParameterAsText(3)) # RENAME
        pallette = arcpy.GetParameterAsText(4)
        var_filter = arcpy.GetParameterAsText(5) # THINK OF MAKING UNIVERSAL
        cell_size = int(arcpy.GetParameterAsText(6))
        search_dist = int(arcpy.GetParameterAsText(7))
        output_dir = arcpy.GetParameterAsText(8) # IF NULL - USE ENV WORKDIR

        # Constants
        decimals = 1
        project = "CURRENT"
        main_map = "Map"
        aprx = arcpy.mp.ArcGISProject(project)
        mp = aprx.listMaps(main_map)[0]
        lyr = mp.listLayers(polygons_layer_name)[0]

        lyt = None
        for l in aprx.listLayouts():
            if l.name == layout_name:
                lyt = l
                break

        map_frame_name = "Map Frame"
        map_frame = lyt.listElements("MAPFRAME_ELEMENT", map_frame_name)[0]
        map_layer = map_frame.map.listLayers(polygons_layer_name)[0]
        pathway = os.path.join(output_dir, var_filter, symbology_mode)
        os.makedirs(pathway, exist_ok = True)

        # dictionary for variable labels
        substitutions = {
            'p_NISP': 'Pottery Total (NISP)',
            'f_MNI_Total' : 'Flints Total (MNI)',
            'f_NISP_Total': 'Flints Total (NISP)',
            'f_MNI_CTEs': 'CTE-s (MNI)',
            'f_MNI_Cores': 'Cores (MNI)',
            'f_MNI_Tools': 'Tools (MNI)',
            'f_MNI_OtherTypes': 'Other Flint Types (MNI)',
            'f_Total_Patinated': 'Patinated Flints (NISP)',
            'f_Total_Fresh': 'Fresh Flints (NISP)',
            'f_Total_Burnt': 'Burnt Flints (NISP)',
            'f_Total_not_burnt': 'Not Burnt Flints (NISP)',
            'f_Broken': 'Broken Flints (NISP)',
            'f_Complete': 'Complete Flints (NISP)',
            'date_f_LP': 'Lower Paleolithic (NISP)',
            'date_f_MP': 'Middle Paleolithic (NISP)',
            'date_f_UP_EPI': 'Blades and Bladelets (NISP)',
            'date_f_NEO_CHAL': 'Neolithic - Chalcolithic (NISP)',
            'date_p_LB_I': 'Late Bronze - Iron (NISP)',
            'date_p_Classic': 'Classic Period (NISP)',
            'date_p_OTT': 'Ottoman Period (NISP)'
            }
        
        # get the list of variables based on mode selected
        varlist = varlist_generator(lyr, var_filter)

        #if symbology_mode == 'min to max':
        classBreakValues = None
        classBreakLabels = None
        color_scheme = None
        if symbology_mode in ["min to max", "natural breaks"]:
            classBreakValues, classBreakLabels, color_scheme = calculateClassBreaksAndColors(project, main_map, polygons_layer_name, varlist, breaks, decimals, pallette)
        elif symbology_mode == "heatmap":
            generate_outline(project, main_map, polygons_layer_name, 50)

        ## LEGEND ADJUSTING ##
        # Loop through the variables in the list, update the symbology, and export each layout
#        varlist = ["f_NISP_Total"] # REMOVE
        for var in varlist:
            arcpy.AddMessage(f"Updating symbology for: {var}")
            if symbology_mode == "heatmap":
                heatmapSymbology(project, main_map, polygons_layer_name, points_layer_name, var, cell_size, search_dist, breaks, "research_area")
            else:
                setSymbology(project, main_map, polygons_layer_name, classBreakValues, classBreakLabels, color_scheme, var, symbology_mode, pallette)
            #aprx.save()


    ## EXPORTING ##
            # Export the layout to a TIFF file

            safe_var = var.replace(" ", "_").replace("/", "_").replace("\\", "_")
            tiff_output = os.path.join(pathway, f"Maharal_{safe_var}_{symbology_mode.replace(' ', '')}.tif")

            lyt.exportToTIFF(
                tiff_output,
                resolution = 300,
                tiff_compression = "LZW"
                )
            arcpy.AddMessage(f"Exported layout to {tiff_output}")

        arcpy.AddMessage("All layouts exported successfully.")

        ## Clean up in-memory
        arcpy.management.Delete("in_memory")
        collected = gc.collect()
        arcpy.AddMessage(f"Garbage collector: collected {collected} objects." )
    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages())

if __name__ == "__main__":
    main()