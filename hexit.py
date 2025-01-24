"""
Script Name: hexit.py
Authors: Mitya Kletzerman, Yitzchak Jaffe
Version: 1.0
Date: 2025-01-11
Description:
    This script performs spatial analysis and tessellation processing for a given set of input data. 
    It automates the generation of a research area outline, creates hexagonal tessellations, and calculates 
    spatial statistics such as elevation and slope for each tessellation tile. Additionally, the script generates 
    summary statistics tables and allows exporting data to Excel format. Users can optionally add the output layers 
    to the current ArcGIS Pro map session.

    The key functionalities include:
    - Validation of input feature layers to ensure spatial integrity.
    - Creation of a buffered "Area of Investigation" outline if not provided.
    - Generation of variable lists for analysis while excluding untouchable fields.
    - Processing tessellation (e.g., clipping, removing irregular shapes, adding calculated attributes).
    - Extracting elevation and slope data from a Digital Elevation Model (DEM).
    - Updating or creating summary tables with statistical measures (e.g., min, max, mean, median).
    - Exporting summary and spatial data to Excel format.
    - Optional addition of output layers to the ArcGIS Pro map.

Dependencies:
    - arcpy: For ArcGIS geoprocessing and analysis.
    - os: For handling file paths and directories.

Usage:
    This script is designed to be run in an ArcGIS Pro environment with the following required parameters:
    1. Input feature class containing spatial data (e.g., survey points or density data).
    2. Research area feature class (optional; will be generated if not provided).
    3. Size of tessellation tiles (in square meters).
    4. Path to the Digital Elevation Model (optional).
    5. Option to convert tessellation hexagons to points ("true" or "false").
    6. Option to add output layers to the ArcGIS Pro map ("true" or "false").
    7. Option to save summary and spatial tables to an Excel file ("true" or "false").
    8. Output directory for saving Excel files (if enabled).

Notes:
    - This script is designed to process geospatial data in ArcGIS Pro and assumes the user has set up a proper workspace.
    - Outputs are saved as in-memory layers unless explicitly copied or exported.
    - Ensure that the input data and spatial references are consistent to avoid processing errors.

"""

import arcpy
import os

def validate_feature_layer(feature_layer):
    """
    Description:
        Validates a feature layer by checking its existence, feature count, and spatial reference.
    Args:
        feature_layer (str): Path to the feature layer to validate.
    Returns:
        bool: True if the layer is valid, False otherwise.
    """
    try:
        # Existence check
        if not arcpy.Exists(feature_layer):
            arcpy.AddError(f"{feature_layer} does not exist.")
            return False
        
        # Feature count check
        rc = arcpy.management.GetCount(feature_layer)
        if rc == 0:
            arcpy.AddError(f"{feature_layer} has no features.")
            return False
        
        # Spatial reference check
        spatial_ref = arcpy.Describe(feature_layer).spatialReference
        if spatial_ref.name == "Unknown":
            arcpy.AddError(f"Spatial reference for {feature_layer} is not defined.")
            return False

        return True

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def generate_outline(input, bufsize):
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
        arcpy.AddMessage(f"Generating 'Area of Investigation' for '{input}' with smoothing factor {bufsize}")
        output = "research_area"

        # Buffer the input features to create a buffered polygon
        buffered_polygon = "in_memory/buffered_polygon"
        arcpy.analysis.Buffer(
            in_features = input,
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
            buffer_distance_or_field = f"-{bufsize} Meters",
            line_side = "FULL",
            line_end_type = "ROUND",
            dissolve_option = "NONE",
            dissolve_field = None,
            method = "GEODESIC"
        )

        # Save the result to the output
        arcpy.management.CopyFeatures(debuffered_polygon, output)
        arcpy.AddMessage(f"Area of Investigation '{output}' has been generated.")

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def area_counting(extent, input, output, unit):
    try:
        def get_total_area(data):
            total_area = 0
            with arcpy.da.SearchCursor(data, ["SHAPE@"]) as cursor:
                for row in cursor:
                    # Add the area of each polygon
                    total_area += row[0].area
            if unit == "km2":
                return (f"{round(total_area/1000000, 2)} km2")
            else:
                return (f"{round(total_area, 0)} m2")
    
        arcpy.AddMessage(f"\n Area of investigation: {get_total_area(extent)}\n Input collection units total area: {get_total_area(input)}\n Hexagons total area: {get_total_area(output)}\n\n")

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def varlist_generator(input, surv_field, vis_field):
    """
    Description:
        Creates a list of variable field names from the input feature layer, excluding technical fields.
    Args:
        input (str): Path to the input feature layer.
        surv_field (str): Name of the surveyor number field.
        vis_field (str): Name of the surface visibility field.
    Returns:
        list: List of variable field names.
    """
    try:
        arcpy.AddMessage(f"Creating a list of variables for {input}")
        
        # Exclude untouchable fields
        untouchable = ["OBJECTID", "Shape", "Shape_Length", "Shape_Area", "CU", f"{surv_field}", f"{vis_field}"]
        fields = arcpy.ListFields(input)
        varlist = [field.name for field in fields if field.name not in untouchable]
        
        arcpy.AddMessage(f"Variable list: {varlist}.")
        return varlist

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def remove_prefixes(field_name):
    """
    Description:
    Removes predefined prefixes from the given field name.
    Args:
        field_name (str): The original field name.
    Returns:
        str: Field name without the prefix.
    """
    try:
        prefixes = ["SUM_", "hadash_hexes_", "elevations_table_", "slopes_table_"]
        for prefix in prefixes:
            if field_name.startswith(prefix):
                return field_name[len(prefix):]

        return field_name

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def tesselation(input, spatial_extent, tesselation_shape_type, tile_area, variable_list):
    """
    Description:
        Generates a tessellation within the given spatial extent and processes variable values.
    Args:
        input (str): Path to the input feature polypolygon layer.
        spatial_extent (str): Path to the spatial extent feature class.
        tesselation_shape_type (str): Shape type for tessellation.
        tile_area (float): Area of each tessellation tile in square meters.
        variable_list (list): List of variable fields to process.
    Returns:
        str: Path to the processed tessellation feature class.
    """
    try:
        arcpy.AddMessage(f"Creating and processing tessellation with tile size {tile_area} square meters")
        
        # Generate the tessellation
        spatial_ref = arcpy.Describe(spatial_extent).spatialReference
        temp_tessellation = "in_memory/temp_tessellation"
        arcpy.management.GenerateTessellation(
            Output_Feature_Class = temp_tessellation,
            Extent = spatial_extent,
            Shape_Type = tesselation_shape_type,
            Size = tile_area,
            Spatial_Reference = spatial_ref
        )

        # Clip tessellation to the research area
        clipped_tessellation = "in_memory/clipped_tessellation"
        arcpy.analysis.Clip(
            temp_tessellation,
            spatial_extent,
            clipped_tessellation
            )

        # Add a field to store the shape area value and make a feature layer
        arcpy.management.AddField(clipped_tessellation, "AREA", "DOUBLE")
        arcpy.management.CalculateGeometryAttributes(clipped_tessellation, [["AREA", "AREA"]], "", "", spatial_ref)
        arcpy.management.MakeFeatureLayer(clipped_tessellation, "clipped_lyr")

        # Remove ragged hexagons
        arcpy.management.SelectLayerByAttribute("clipped_lyr", "NEW_SELECTION", f"\"AREA\" < {tile_area - 1}")
        arcpy.management.DeleteFeatures("clipped_lyr")
        arcpy.management.DeleteField("clipped_lyr", "AREA")

        # create FL with parts of hexagons overlapping the collection units
        intersected_hexes = "in_memory/intersected_hexes"
        arcpy.analysis.Intersect(
            in_features = ["clipped_lyr", input],
            out_feature_class = intersected_hexes,
            join_attributes = "ALL",
            cluster_tolerance = None,
            output_type = "INPUT"
        )

        # transorm back to number of items
        intersected_ind_hexes = "in_memory/intersected_ind_hexes"
        arcpy.management.AddField(intersected_hexes, "AREA", "DOUBLE")
        arcpy.management.CalculateGeometryAttributes(intersected_hexes, [["AREA", "AREA"]], "", "", spatial_ref)
        arcpy.management.MakeFeatureLayer(intersected_hexes, intersected_ind_hexes)
        with arcpy.da.UpdateCursor(intersected_ind_hexes, ['AREA'] + variable_list) as cursor:
            for row in cursor:
                shape_area = row[0]
                if shape_area != 0:
                    for i in range(1, len(row)):
                        row[i] = row[i] * shape_area
                    cursor.updateRow(row)
                else:
                    pass

        arcpy.management.DeleteField(intersected_ind_hexes, "AREA")
        
        # create FL with parts of hexagons NOT overlapping the collection units
        not_intersected_hexes = "in_memory/not_intersected_hexes"
        arcpy.analysis.Erase(
            in_features = "clipped_lyr",
            erase_features = intersected_ind_hexes,
            out_feature_class  = not_intersected_hexes,
            cluster_tolerance = None
        )
        
        # add variables and set them to zero
        for field in variable_list:
            arcpy.management.AddField(not_intersected_hexes, field, "DOUBLE")

        with arcpy.da.UpdateCursor(not_intersected_hexes, variable_list) as cursor:
            for row in cursor:
                # Set all fields to 0
                for i in range(len(row)):
                    row[i] = 0
                cursor.updateRow(row)

        merged_hexes = "in_memory/merged_hexes"
        arcpy.management.Merge([not_intersected_hexes, intersected_ind_hexes], merged_hexes)

        stats_fields = [(field, "SUM") for field in variable_list]

        hadash_hexes = "in_memory/hadash_hexes"
        arcpy.management.Dissolve(
            in_features = merged_hexes,
            out_feature_class = hadash_hexes,
            dissolve_field = "GRID_ID",
            statistics_fields = stats_fields
        )
  
        hex_field_list = arcpy.ListFields(hadash_hexes)

        for field in hex_field_list:
            if field.name.startswith("SUM_"):
                new_name = remove_prefixes(field.name)
                if new_name != field.name:
                    arcpy.management.AlterField(hadash_hexes, field.name, new_name, new_name)
        
        # switch back to dencities
        arcpy.management.AddField(hadash_hexes, "AREA", "DOUBLE")
        arcpy.management.CalculateGeometryAttributes(hadash_hexes, [["AREA", "AREA"]], "", "", spatial_ref)
        with arcpy.da.UpdateCursor(hadash_hexes, ['AREA'] + variable_list) as cursor:
            for row in cursor:
                shape_area = row[0]
                if shape_area != 0:
                    for i in range(1, len(row)):
                        row[i] = row[i] / shape_area
                    cursor.updateRow(row)
                else:
                    pass
        arcpy.management.DeleteField(hadash_hexes, "AREA")
        arcpy.AddMessage("Tessellation processing finished.")

        return hadash_hexes

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def elevation_and_slope(digital_elevation_model, hexagon_layer):
    """
    Description:
        Extracts elevation and slope values from a digital elevation model (DEM) and calculates mean values per hexagon tile.
    Args:
        digital_elevation_model (str): Path to the DEM raster.
        hexagon_layer (str): Path to the polypolygon feature class.
    Returns:
        str: Path to the polypolygon feature class with elevation and slope attributes.
    """
    try:
        arcpy.AddMessage(f"Fetching elevation and slope data from {digital_elevation_model}")
        
        # Get elevation & slope mean values per feature
        if digital_elevation_model is not None:
            cell_width = digital_elevation_model.meanCellWidth
            cell_height = digital_elevation_model.meanCellHeight
            new_cell_size = str(cell_width / 10) + " " + str(cell_height / 10)
            
            resampled_elevation = "in_memory/resampled_elevation"
            arcpy.management.Resample(
                in_raster = digital_elevation_model,
                out_raster = resampled_elevation,
                cell_size = new_cell_size,
                resampling_type = "NEAREST"
            )
            
            slope_raster = "in_memory/slope_raster"
            arcpy.ddd.Slope(
                in_raster = digital_elevation_model,
                out_raster = slope_raster,
                output_measurement = "DEGREE",
                z_factor = 1,
                method = "GEODESIC",
                z_unit = "METER",
                analysis_target_device = "GPU_THEN_CPU"
            )
            
            resampled_slope = "in_memory/resampled_slope"
            arcpy.management.Resample(
                in_raster = "slope_raster",
                out_raster = resampled_slope,
                cell_size = new_cell_size,
                resampling_type = "NEAREST"
            )
            
            arcpy.sa.ZonalStatisticsAsTable(
                in_zone_data = hexagon_layer,
                zone_field = "GRID_ID",
                in_value_raster = resampled_elevation,
                out_table = "elevations_table",
                ignore_nodata = "DATA",
                statistics_type = "MEAN",
                process_as_multidimensional = "CURRENT_SLICE",
                percentile_values = 90,
                percentile_interpolation_type = "AUTO_DETECT",
                circular_calculation = "ARITHMETIC",
                circular_wrap_value = 360
            )
            
            arcpy.sa.ZonalStatisticsAsTable(
                in_zone_data = hexagon_layer,
                zone_field="GRID_ID",
                in_value_raster = resampled_slope,
                out_table = "slopes_table",
                ignore_nodata="DATA",
                statistics_type="MEAN",
                process_as_multidimensional="CURRENT_SLICE",
                percentile_values=90,
                percentile_interpolation_type="AUTO_DETECT",
                circular_calculation="ARITHMETIC",
                circular_wrap_value=360,
                out_join_layer="HEX_SLOPE"
            )
            
            # Rename elevation & slope variables
            arcpy.management.AlterField("elevations_table", "MEAN", "Elevation", "Elevation")
            arcpy.management.AlterField("slopes_table", "MEAN", "Slope", "Slope")
            
            # Join elevation & slope with the feature layer
            arcpy.management.MakeFeatureLayer(hexagon_layer, "clipped_lyr")
            arcpy.management.AddJoin("clipped_lyr", "GRID_ID", "elevations_table", "GRID_ID", "KEEP_COMMON")
            arcpy.management.AddJoin("clipped_lyr", "GRID_ID", "slopes_table", "GRID_ID", "KEEP_COMMON")
            
            # Clean feature layer field names
            equipped_hexes = "in_memory/equipped_hexes"
            arcpy.management.CopyFeatures("clipped_lyr", equipped_hexes)
            remove_fields = [
                "elevations_table_OBJECTID", "elevations_table_GRID_ID", "elevations_table_ZONE_CODE",
                "elevations_table_COUNT", "elevations_table_AREA", "slopes_table_OBJECTID", "slopes_table_GRID_ID",
                "slopes_table_ZONE_CODE", "slopes_table_COUNT", "slopes_table_AREA"
                ]
            
            fields = arcpy.ListFields(equipped_hexes)
            for field in fields:
                if field.name in remove_fields:
                    arcpy.management.DeleteField(equipped_hexes, field.name)
            
            fields = arcpy.ListFields(equipped_hexes)
            for field in fields:
                new_name = remove_prefixes(field.name)
                if new_name != field.name:
                    arcpy.management.AlterField(equipped_hexes, field.name, new_name, new_name)
            
            # Output
            arcpy.AddMessage(f"Slope and elevation have been added to the feature layer '{equipped_hexes}'")
            return equipped_hexes

        else:
            arcpy.AddError(f"There is some problem with '{digital_elevation_model}'.")
            return None

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def table_generator(sum_table, input, output, varlist):
    """
    Description:
        Updates or creates a summary table containing statistics (e.g., min, max, mean, median) for variables.

    Args:
        sum_table (str): Path to the summary table (Supposed to be created due to subset tool first).
        input (str): Path to the input feature class.
        output (str): Path to the output feature class.
        varlist (list): List of variables for which statistics are generated.

    Returns:
        None
    """
    try:
        arcpy.AddMessage(f"Updating the summary table '{sum_table}'")
        
        # Check if SUM table exists and create it if it doesn't
        if not arcpy.Exists(sum_table):
            arcpy.AddWarning(f"'{sum_table}' was not found. Please run the 'Subset' script to get the 'SUM' table.")
            arcpy.AddMessage(f"Create the summary table '{sum_table}'")
            
            # Create the table and add fields if it doesn't exist
            table = arcpy.management.CreateTable("in_memory", "var_sums")
            arcpy.management.AddField(table, "Field_Name", "TEXT")
            arcpy.management.AddField(table, "SUM_Adjusted", "DOUBLE", 2)

            # Insert data into the summary table
            with arcpy.da.InsertCursor(table, ["Field_Name", "SUM_Adjusted"]) as cursor:
                for field_name in varlist:
                    sum_value = 0
                    with arcpy.da.SearchCursor(input, [field_name]) as search_cursor:
                        for row in search_cursor:
                            sum_value += row[0] if row[0] is not None else 0
                    cursor.insertRow([field_name, sum_value])
            arcpy.management.CopyRows(table, "SUM")
        
        # Updating SUM table
        fields_to_add = [
            ("MIN", "DOUBLE", 1),
            ("MAX", "DOUBLE", 1),
            ("MEAN", "DOUBLE", 2),
            ("MEDIAN", "DOUBLE", 2)
        ]
        existing_fields = [field.name for field in arcpy.ListFields(sum_table)]
        
        # Check and add fields if they don't exist
        for field_name, field_type, field_precision in fields_to_add:
            if field_name not in existing_fields:
                arcpy.management.AddField(sum_table, field_name, field_type, field_precision=field_precision)

        with arcpy.da.UpdateCursor(sum_table, ["Field_Name", "MIN", "MAX", "MEAN", "MEDIAN"]) as update_cursor:
            for row in update_cursor:
                field_name = row[0]  # Get the field name
                values = []

                # Calculate statistics by searching through the original feature class
                with arcpy.da.SearchCursor(output, [field_name]) as search_cursor:
                    for search_row in search_cursor:
                        if search_row[0] is not None:
                            values.append(search_row[0])
                
                if values:
                    # Calculate min, max, mean, and median
                    min_value = min(values)
                    max_value = max(values)
                    mean_value = sum(values) / len(values)
                    median_value = sorted(values)[len(values) // 2] if len(values) % 2 == 1 else \
                                (sorted(values)[len(values) // 2 - 1] + sorted(values)[len(values) // 2]) / 2.0
                else:
                    min_value = max_value = mean_value = median_value = None

                # Update the statistics in the table
                row[1] = min_value
                row[2] = max_value
                row[3] = mean_value
                row[4] = median_value

                update_cursor.updateRow(row)
        
        arcpy.AddMessage(f"Table '{sum_table}' has been updated.")

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def path_checker(path):

    try:
        # Check if the path is a directory
        if os.path.isdir(path):
            test_file = os.path.join(path, 'test_write_permission.txt')
            with open(test_file, 'w') as f:
                f.write('Testing write permission.')
            os.remove(test_file)  # Remove the test file after the check
        else:
            return False
        return True
    except (PermissionError, IOError):
        return False

def table_exporter(summary_table, input_table, output_table, save_path, filename):
    """
    Description:
        Exports specified tables to an Excel file.
    Args:
        summary_table (str): Path to the summary table.
        input_table (str): Path to the input table.
        output_table (str): Path to the output table.
        save_path (str): Directory path to save the Excel file.
        filename (str): Name of the Excel file.
    Returns:
        None
    """
    try:
        arcpy.AddMessage(f"Saving the tables {summary_table}, {input_table}, and {output_table}")
        if path_checker(save_path):
            ## export tables ##
            arcpy.conversion.TableToExcel(                                                  
                Input_Table = f"{input_table}; {output_table}; {summary_table}",
                Output_Excel_File = os.path.join(save_path, f"{filename}.xlsx"),
                Use_field_alias_as_column_header = "ALIAS",
                Use_domain_and_subtype_description = "CODE"
            )
            arcpy.AddMessage(f"Data and summary tables were saved to:\n{save_path}\{filename}.xlsx")
        else:
            arcpy.AddWarning(f"The specified folder is protected from writing or doesn't exist:\n{save_path}\{filename}.xlsx")
    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def add_to_the_map(feature_layer, workspace, project):
    """
    Description:
        Adds a feature layer to the active map in the current ArcGIS Pro project.
    Args:
        feature_layer (str): Path to the feature layer to add.
    Returns:
        None
    """
    try:
        arcpy.AddMessage(f"Adding '{feature_layer}' to the active map of the current project")
        # Adding feature layers to the map
        project.save()
        map_doc = project.activeMap
        map_doc.addDataFromPath(os.path.join(workspace, feature_layer))
        arcpy.AddMessage(f"'{feature_layer}' has been added to the map.")     

    except ValueError as e:
        arcpy.AddError(str(e))
        raise

def main():
    """
    Main function to control the execution workflow of the script.
    Args:
        Parameter inputs are passed from user through ArcGIS tool.
    Returns:
        None
    """
    try:
        # Parameter inputs from user 
        in_dencities = arcpy.GetParameterAsText(0)  # Multi-polygon feature class storing geometry and data on colelction units
        research_area = arcpy.GetParameterAsText(1) # Single-polygon feature class defining the research extent
        size = float(arcpy.GetParameterAsText(2)) # Tile size in m2 for the tesselation
        dem_path = arcpy.GetParameterAsText(3) # Path to the raster digital elevation model that covers the research area
        points_mode = arcpy.GetParameterAsText(4) # True if user want the final output to be transformed into poly-point feature class 
        add_to_map = arcpy.GetParameterAsText(5) # True if user wants to add the output to the active map in current project
        save_output_tables = arcpy.GetParameterAsText(6) # True if user wants to save the output tables
        output_dir = arcpy.GetParameterAsText(7) # Path to the output folder
        
        # Constants
        environment = arcpy.env.workspace
        aprx = arcpy.mp.ArcGISProject("CURRENT")
        output_fc = "Hexagons" # The name of output multi-polygon tesselation feature
        shape_type = "HEXAGON" # Shape type for the tessalation
        surveyor_number_field = "Number_of_surveyors"  # Field name for the number of surveyors participated in items collection
        visibility_field = "Surface_Visibility" # Field name for surface visibility
        sum_table_name = "SUM" # The name of the summary table that is supposed to be produced with Subset tool
        project_name = "Maharal-2025" # The name for the output Escel file with tables
        buffer_size = 100 # The buffer size used to create an Area of Investigation for case when not provided by user
        units = "km2" # Km2 gives information message in km2
        # Validate the input feature layer
        if not validate_feature_layer(in_dencities):
            arcpy.AddError(f"Validation of input feature '{in_dencities}' layer failed. Stopping execution")
            return
    
        # Generate research area if not provided
        if not research_area:
            generate_outline(in_dencities, buffer_size)
            research_area = "research_area"

        # Generate a list of variables
        varlist = "in_memory/varlist"
        varlist = varlist_generator(in_dencities, surveyor_number_field, visibility_field)

        # Create tessellation
        hadash_hexes = "in_memory/hadash_hexes"
        hadash_hexes = tesselation(in_dencities, research_area, shape_type, size, varlist)
     
        # Add elevation and slope data if DEM is provided
        if dem_path:
            dem_raster = arcpy.sa.Raster(dem_path)
            equipped_hexes = elevation_and_slope(dem_raster, hadash_hexes)
        else:
            equipped_hexes = hadash_hexes
            arcpy.AddWarning("DEM not provided. Skipping elevation and slope calculations")
        
        # Validate output features
        if not validate_feature_layer(equipped_hexes):
            arcpy.AddError(f"Validation of output feature layer '{equipped_hexes}' failed. Stopping execution")
            return
        
        # Copy features to final output
        arcpy.management.CopyFeatures(equipped_hexes, output_fc)        
        arcpy.AddMessage(f"Feature layer '{output_fc}' has been created successfully.")

        # Convert hexagons to points if required
        if points_mode == "true":
            points_hexagons = "in_memory/points_hexagons"
            arcpy.management.FeatureToPoint(
                in_features = output_fc,
                out_feature_class = points_hexagons,
                point_location = "CENTROID"
            )

            arcpy.management.CopyFeatures(points_hexagons, f"{output_fc}_points")   
            arcpy.AddMessage(f"Feature layer '{output_fc}_points' has been created successfully.")
        
        # Count areas
        area_counting(research_area, in_dencities, output_fc, units)

        # Generate summary table
        table_generator(sum_table_name, in_dencities, output_fc, varlist)

        # Export tables if required
        if save_output_tables == "true":
            if len(output_dir) == 0:
                arcpy.AddWarning("Output directory not specified. Tables were not saved.")
            else:
                table_exporter(sum_table_name, in_dencities, output_fc, output_dir, project_name)
        
        # Add layers to the map if required
        if add_to_map == "true":
            add_to_the_map(output_fc, environment, aprx)

            if points_mode == "true":
                add_to_the_map(f"{output_fc}_points", environment, aprx)

        # Clean up
        arcpy.management.Delete("equipped_hexes")
        arcpy.management.Delete("slopes_table")
        arcpy.management.Delete("elevations_table")
        arcpy.management.Delete("in_memory")
        arcpy.management.Delete("slope_raster")

    except arcpy.ExecuteError:
        arcpy.AddError(arcpy.GetMessages())

if __name__ == "__main__":
    main()