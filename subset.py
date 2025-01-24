"""
Script Name: subset.py
Authors: Mitya Kletzerman, Yitzchak Jaffe
Version: 1.1
Date: 2025-01-11
Description:
    This script is designed to process and analyze spatial data related to archaeological collection units.
    It integrates polygon geometries, surveyor information, surface visibility, and artifact counts from 
    Excel data. The script performs several key functions:

    - Validates and processes polygon features, including clipping to research extents and resolving overlaps.
    - Integrates artifact data by joining Excel-based counts to polygon features.
    - Applies adjustments to artifact counts based on surface visibility and the number of surveyors.
    - Cleans output data by removing unnecessary prefixes and renaming fields.
    - Generates a summary table containing sums of variables, with both original and adjusted values.
    - Checks for collection units with zero pottery or flint counts and generates warnings.
    - Optionally adds the processed output layers to the active ArcGIS Pro map.

    This script is suitable for automating the analysis of archaeological survey data and improving data accuracy
    through visibility and surveyor-based adjustments.

Dependencies:
    - arcpy: For geoprocessing and spatial analysis.
    - os: For file and directory path operations.
    - re: For regular expression-based field name cleaning.

Usage:
    This script is designed to be run in an ArcGIS Pro environment with the following parameters:
    1. Input polygon feature class (collection units).
    2. Research area feature class (optional; defines the clipping extent).
    3. Path to an Excel file containing artifact counts.
    4. Surface visibility multiplier (a value ≤ 2).
    5. Option to add the output layer to the ArcGIS Pro map ("true" or "false").

Outputs:
    - A cleaned and processed feature class named "Densities" with adjusted artifact counts.
    - A summary table "SUM" stored in memory, showing original and adjusted variable sums.
    - Optionally, the processed feature class is added to the ArcGIS Pro map.

Notes:
    - Ensure that the input data, spatial references, and field names are consistent to avoid errors.
    - The script generates intermediate datasets in memory to optimize performance and minimize disk usage.
    - Surface visibility adjustments rely on a user-defined multiplier, which should not exceed 2.
    - Zero counts in pottery or flint fields are checked and reported for further investigation.
"""

import arcpy
import os
import re

def multiplier_checker(value):
    """
    Description:
        Checks and returns the surface visibility multiplier if it's valid (<=2).
    Args:
        value (str): The surface visibility multiplier input as a string.
    Returns:
        float: The validated multiplier value.
    """
    try:
        num_value = float(value)
        if abs(num_value) <= 2:
            return num_value
        else:
            raise ValueError("Surface visibility multiplier should be less than or equal to 2")
    
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

def fetch_excel_data(path):
    """
    Description:
        Fetches data from an Excel file.
    Args:
        path (str): The file path to the Excel file.
    Returns:
        A list of variable names extracted from the table and the in-memory table of the items.
    """
    try:
        arcpy.AddMessage(f"Fetching excel data from {path}")
        ## Creating a table out of excel data file
        temp_data = "in_memory/temp_data"
        varlist = "in_memory/varlist"
        arcpy.conversion.ExcelToTable(Input_Excel_File = path, Output_Table = temp_data)
        
        # Check if temp_data has rows
        row_count = int(arcpy.GetCount_management(temp_data)[0])
        if row_count == 0 or temp_data is None or varlist is None:
            arcpy.AddError("Error: The data could not be fetched, or the table is empty.")
        
        fields = arcpy.ListFields(temp_data)
        varlist = [field.name for field in fields[2:]]

        arcpy.AddMessage("Data has been successfully fetched from Excel file.")
        arcpy.AddMessage(f"Variable list: {varlist}.")
        
        # Return the varlist and temp_data if everything is fine
        return varlist, temp_data
    
    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def process_polygons(polygons, extention, id_field_for_cu, surveyours_field_name, visibility_field_name):
    """
    Description:
        Processes a set of input polygons to refine and prepare them for further analysis by performing clipping, splitting, and de-overlaying operations.
    Args:
        polygons: The input feature class or layer of polygons to be processed.
        extention: The feature class or layer defining the research area boundary for clipping.
        id_field_for_cu: The ID field.
        surveyours_field_name: Field name for the number of surveyors participated in items collection
        visibility_field_name: Field name for surface visibility
    Returns:
        The feature class of processed polygons.
    """
    try:
        arcpy.AddMessage(f"Processing '{polygons}'")
        ## Make a copy of the input polygons to work with
        temp_polygons = "in_memory/temp_polygons"
        arcpy.management.CopyFeatures(polygons, temp_polygons)

        ## Clip the polygons to the research area
        clipped_polygons = "in_memory/clipped_polygons"
        arcpy.analysis.Clip(temp_polygons, extention, clipped_polygons)

        ## Split any polygons that cross the boundary of the research area
        split_polygons = "in_memory/split_polygons"
        arcpy.analysis.Intersect([clipped_polygons, extention], split_polygons, "ALL")
        
        ## Remove overlaying of initial polygons over each other
        overlays_polygons = "in_memory/overlays_polygons"
        arcpy.analysis.Intersect(split_polygons, overlays_polygons)
        arcpy.management.DeleteIdentical(overlays_polygons, ["Shape"])

        revealed_polygons = "in_memory/revealed_polygons" 
        arcpy.analysis.Erase(split_polygons, overlays_polygons, revealed_polygons)

        separated_polygons = "in_memory/separated_polygons"
        arcpy.management.Merge(
            [revealed_polygons, overlays_polygons],
            separated_polygons,
            add_source = "NO_SOURCE_INFO"
            )

        good_polygons = "in_memory/good_polygons"
        arcpy.management.Dissolve(
            in_features = separated_polygons,
            out_feature_class = good_polygons,
            dissolve_field = id_field_for_cu,
            statistics_fields = f"{surveyours_field_name} MEAN; {visibility_field_name} MEAN",
            multi_part = "MULTI_PART",
            unsplit_lines = "DISSOLVE_LINES",
            concatenation_separator = ""
            )

        arcpy.AddMessage(f"Subsetting finished. Resultant '{good_polygons}' were de-overlapped.")
        return good_polygons
    
    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def remove_prefixes(field_name):
    """
    Description:
        Removes specified prefixes from a field name.
    Args:
        The name of the field to process.
    Returns:
        The field name with specified prefixes removed, if present.
    """
    try:
        prefixes = ["good_polygons_MEAN_", "good_polygons_", "temp_data_"]
        prefix_regex = re.compile(f"^({'|'.join(map(re.escape, prefixes))})")
        
        return prefix_regex.sub("", field_name)
    
    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def clean_output(output, id_field_for_cu):
    """
    Description:
        Cleans and renames fields in the output feature class by removing prefixes and renaming specified fields.
    Args:
        output (str): The path to the output feature class to be cleaned.
        id_field_for_cu (str): The field name used as the identifier for collection units.
    Returns:
        str: The cleaned and updated feature class with renamed fields.
    """
    try:
        data = output
        arcpy.management.AlterField(data, f"good_polygons_{id_field_for_cu}", id_field_for_cu, id_field_for_cu)
        reserved = ["OBJECTID", "Shape", "Shape_Length", "Shape_Area", id_field_for_cu]
        fields = arcpy.ListFields(data)

        new_fields = []
        for field in fields:
            if field.name not in reserved:
                new_name = remove_prefixes(field.name)
                if new_name != field.name:
                    arcpy.management.AddField(data, new_name, "DOUBLE")
                    new_fields.append((field.name, new_name))
        
        for old_name, new_name in new_fields:
            arcpy.management.CalculateField(data, new_name, f"!{old_name}!", "PYTHON3")
            arcpy.management.DeleteField(data, old_name)

        arcpy.AddMessage("Cleaning completed.")
        return data

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def create_sum_table(digits_after_dot, list_of_variables, output, adjusted_output):
    """
    Description:
        Creates a summary table of variable sums, including the original and adjusted values.
    Args:
        digits_after_dot (int): The number of decimal places for the sum fields.
        list_of_variables (list): A list of variable (field) names to include in the summary.
        output (str): The feature class or table containing the original data.
        adjusted_output (str): The feature class or table containing the adjusted data.
    Returns:
        None: The function creates an in-memory summary table and saves it as "SUM".
    """
    try:
        arcpy.AddMessage("Creating a summary table")
        # Creating a summary table #
        sum_table = arcpy.management.CreateTable("in_memory", "var_sums")

        arcpy.management.AddField(sum_table, "Field_Name", "TEXT")
        arcpy.management.AddField(sum_table, "SUM_Original", "DOUBLE", digits_after_dot)
        arcpy.management.AddField(sum_table, "SUM_Adjusted", "DOUBLE", digits_after_dot)

        # Initialize an insert cursor to add rows to the in-memory table
        with arcpy.da.InsertCursor(sum_table, ["Field_Name", "SUM_Original"]) as cursor:
            for field_name in list_of_variables:
                sum_value = 0
                with arcpy.da.SearchCursor(output, [field_name]) as search_cursor:
                    for row in search_cursor:
                        sum_value += row[0] if row[0] is not None else 0
                cursor.insertRow([field_name, sum_value])

        # Count the adjusted sums of items by variable
        with arcpy.da.UpdateCursor(sum_table, ["Field_Name", "SUM_Adjusted"]) as update_cursor:
           for row in update_cursor:
                field_name = row[0]
                sum_value = 0
                with arcpy.da.SearchCursor(adjusted_output, [field_name]) as search_cursor:
                    for search_row in search_cursor:
                        sum_value += round(search_row[0],1) if search_row[0] is not None else 0
                
                row[1] = sum_value  # Update the Sum_Adjusted field
                update_cursor.updateRow(row)
        
        arcpy.management.CopyRows(sum_table, "SUM")
        arcpy.AddMessage("Table 'SUM' has been created.\nYou may inspect the initial and adjusted SUMs there.")

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def transform_vector(vector, multiplier):
    """
    Description:
        Transforms a vector by applying a linear transformation based on a multiplier.
    Args:
        vector (list): A list of values to be transformed.
        multiplier (float): The multiplier applied to the transformation formula.
    Returns:
        list: A transformed list where each value is calculated as `1 + multiplier * index`.
    """
    try:
        transformed = []
        for i in range(len(vector)):
            transformed_value = 1 + multiplier * i
            transformed.append(transformed_value)
        
        return transformed
    
    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def adjuster(visibility_multiplier, output, surv, vis, list_of_variables):
    """
    Description:
        Adjusts values in the dataset based on surface visibility and the number of surveyors.
    Args:
        visibility_multiplier (float): A multiplier used to adjust values based on surface visibility.
        output (str): The path to the dataset to adjust.
        surv (str): The field name representing the number of surveyors.
        vis (str): The field name representing surface visibility.
        list_of_variables (list): A list of variable (field) names to adjust.
    Returns:
        str: The path to the adjusted dataset.
    """
    try:
        data = output

        # Remove the NULL values
        null_count = 0

        # Create an update cursor to iterate over rows
        with arcpy.da.UpdateCursor(data, [f"{surv}", f"{vis}"]) as cursor:
            for row in cursor:
                # Check if either Var1 or Var2 is NULL
                if row[0] is None or row[1] is None:
                    # Increase the counter for NULL values
                    null_count += 1
                    # Delete the row with NULL values
                    cursor.deleteRow()

        # Print the message with the number of rows deleted
        if null_count > 0:
            arcpy.AddWarning(f"{null_count} fearures has been removed because of the NULL values in either '{surv}' or '{vis}'.")

        # ADJUSTING THE VALUES REGARDING SURF AND SURV #
        if visibility_multiplier is not None and visibility_multiplier != 0:
            stat_table = "in_memory/stat_table"
            arcpy.analysis.Statistics(
                in_table = data,
                out_table = stat_table,
                statistics_fields = f"{surv} MIN;{surv} MAX;{vis} MIN;{vis} MAX"
            )

            # Use SearchCursor to access the min/max values from the data table
            with arcpy.da.SearchCursor(stat_table, [f"MIN_{surv}", f"MAX_{surv}", 
                                                    f"MIN_{vis}", f"MAX_{vis}"]) as cursor:
                for row in cursor:
                    min_surveyors, max_surveyors, min_visibility, max_visibility = row
                    
                    # Create a numeric vector (list of integers) from min to max for both fields
                    if visibility_multiplier is not None:
                        surf_range = list(range(int(min_visibility), int(max_visibility) + 1))
                        surf_multipliers = transform_vector(surf_range, abs(visibility_multiplier))
                        
                        if visibility_multiplier < 0:
                            surf_multipliers = surf_multipliers[::-1]
                        
                        # Create a dictionary to map surf values to multipliers
                        surf_multiplier_dict = dict(zip(surf_range, surf_multipliers))
                        
                        # Divide values
                        with arcpy.da.UpdateCursor(data, ['Surface_Visibility'] + list_of_variables) as cursor:
                            for row in cursor:
                                surf_value = row[0]
                                if surf_value in surf_multiplier_dict:
                                    multiplier = surf_multiplier_dict[surf_value]
                                    
                                    # Apply the multiplier to each variable in varlist
                                    for i, var in enumerate(list_of_variables):
                                        row[i + 1] /= multiplier
                                    cursor.updateRow(row)
                                else:
                                    arcpy.AddWarning(f"Surface visibility value '{surf_value}' isn't within a surf_range.")

                        arcpy.AddMessage(f"Surface Visibility Range: {surf_range}.")
                        arcpy.AddMessage(f"Surface Visibility Multipliers: {surf_multipliers}.")
                
        # Divide values by surveyors
        with arcpy.da.UpdateCursor(data, [f"{surv}"] + list_of_variables) as cursor:
            for row in cursor:
                number_of_surveyors = row[0] 
                
                if number_of_surveyors != 0:
                    for i in range(1, len(row)):
                        row[i] = row[i] / number_of_surveyors
                    cursor.updateRow(row)
                else:
                    arcpy.AddWarning("Number_of_surveyors is 0 for one or more rows.")

        arcpy.AddMessage("Values have been adjusted for surveyor number and surface visibility.")
        return data

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def zeroes_checker(output, id_field_for_cu, pottery_field_name, flints_field_name):
    """
    Description:
        Checks for collection units with zero pottery, zero flints, or both, and generates warnings.
    Args:
        output (str): The path to the dataset to check.
        id_field_for_cu (str): The field name representing the collection unit identifier.
        pottery_field_name (str): The field name representing pottery count.
        flints_field_name (str): The field name representing flint count.
    Returns:
        None: The function generates warnings for collection units with zero items.
    """
    try:
        arcpy.AddMessage(f"Checking for zero counts in '{output}'")
        ## Run check for zero items polygons
        zero_pottery_count = 0
        zero_flints_count = 0
        zero_pottery_cu_list = []
        zero_flints_cu_list = []
        zero_both_count = 0
        zero_both_cu_list = []
        
        final_row_count = arcpy.management.GetCount(output)
        arcpy.AddMessage(f"Number of Collection Units in the subset: {final_row_count[0]}.")

        # Use a SearchCursor to iterate through the rows
        with arcpy.da.SearchCursor(output, [id_field_for_cu, pottery_field_name, flints_field_name]) as cursor:
            for row in cursor:
                cu_id = row[0]
                pottery_count = row[1]
                flints_count = row[2]
                
                # Check for zero pottery
                if pottery_count == 0:
                    zero_pottery_count += 1
                    zero_pottery_cu_list.append(cu_id)
                
                # Check for zero flints
                if flints_count == 0:
                    zero_flints_count += 1
                    zero_flints_cu_list.append(cu_id)
                
                # Check for both zero pottery and zero flints
                if pottery_count == 0 and flints_count == 0:
                    zero_both_count += 1
                    zero_both_cu_list.append(cu_id)
        
        # Generate warning messages
        if zero_pottery_count > 0:
            zero_pottery_cu_list_str = ", ".join(str(cu) for cu in zero_pottery_cu_list)
            arcpy.AddWarning(f"There are {zero_pottery_count} Collection Units with zero Pottery:\n{zero_pottery_cu_list_str}")

        if zero_flints_count > 0:
            zero_flints_cu_list_str = ", ".join(str(cu) for cu in zero_flints_cu_list)
            arcpy.AddWarning(f"There are {zero_flints_count} Collection Units with zero Flints:\n{zero_flints_cu_list_str}")
        
        if zero_both_count > 0:
            zero_both_cu_list_str = ", ".join(str(cu) for cu in zero_both_cu_list)
            arcpy.AddWarning(f"There are {zero_both_count} Collection Units with neither Pottery nor Flint items:\n{zero_both_cu_list_str}")
    
    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

def main():
    """
    Description:
        Main function that integrates and executes all processing steps for analyzing collection units.
        It processes input polygon features and Excel data, applies adjustments, generates a summary table, and validates the outputs.
    """
    try:
        # PARAMETERS
        in_polygons = arcpy.GetParameterAsText(0)  # Multi-polygon feature class storing geometry, surveyors, and visibility
        research_area = arcpy.GetParameterAsText(1)  # Single-polygon feature class defining the research extent
        excel_file_path = arcpy.GetParameterAsText(2)  # Path to the Excel file with item counts
        make_adjustments = arcpy.GetParameterAsText(3)
        multiplier_SURF = multiplier_checker(arcpy.GetParameterAsText(4))  # Validate the surface visibility multiplier
        add_to_map = arcpy.GetParameterAsText(5)  # Option to add the output layer to the current map

        # CONSTANTS
        aprx = arcpy.mp.ArcGISProject("CURRENT")  # Setting the current ArcGIS project
        map_doc = aprx.activeMap  # Getting the active map in the project
        out_feature_class = "Densities"  # Name for the final output feature class
        joinFieldCU = "CU"  # Field name for collection unit ID in polygons
        joinFieldData = "Unit"  # Field name for collection unit ID in the Excel table
        surveyor_number_field = "Number_of_surveyors"  # Field name for the number of surveyors participated in items collection
        visibility_field = "Surface_Visibility"  # Field name for surface visibility
        pottery_total_field = "p_NISP"  # Field name for total pottery count
        flints_total_field = "f_NISP_Total"  # Field name for total flint count
        precision = 1  # Precision for decimal places in the summary table
        buffer_size = 100 # The buffer size used to create an Area of Investigation for case when not provided by user

        # Validate the input feature layer
        if not validate_feature_layer(in_polygons):
            arcpy.AddError(f"Validation of input feature '{in_polygons}' layer failed. Stopping execution")
            return

        # Generate research area if not provided
        if not research_area:
            generate_outline(in_polygons, buffer_size)
            research_area = "research_area"

        # FETCHING EXCEL DATA
        varlist, temp_data = fetch_excel_data(excel_file_path)

        # PROCESSING POLYGONS
        processed_polygons = process_polygons(in_polygons, research_area, joinFieldCU, surveyor_number_field, visibility_field)

        # CREATING A FEATURE LAYER AND JOINING EXCEL DATA TO POLYGONS
        rc_features = arcpy.management.GetCount(in_polygons)
        rc_data = arcpy.management.GetCount(temp_data)
        arcpy.AddMessage(f"Joining geometry layer '{processed_polygons}' with the the count data '{temp_data}'")
        arcpy.management.MakeFeatureLayer(processed_polygons, "processed_polygons")
        arcpy.management.AddJoin("processed_polygons", joinFieldCU, temp_data, joinFieldData, join_type = 'KEEP_COMMON')

        # SAVING THE JOINED FEATURE LAYER

        joint_polygons = "in_memory/joint_polygons"
        arcpy.management.CopyFeatures("processed_polygons", joint_polygons)
        rc_final = arcpy.management.GetCount(joint_polygons)
        arcpy.AddMessage(f"Feature layer '{joint_polygons}' has been saved.\n Count of unique Collection Units in geometry feature layer: {rc_features}.\n Count of unique Collection Units in data: {rc_data}.\n Count of unique Collection Units in {joint_polygons}: {rc_final}.")
        
        # CLEANING UP JOINT FEATURE LAYER
        arcpy.AddMessage(f"Сleaning '{joint_polygons}'")
        arcpy.management.DeleteField(joint_polygons, ["temp_data_ObjectID", f"temp_data_{joinFieldData}"])
        pretty_polygons = clean_output(joint_polygons, joinFieldCU)

######################################################################

        arcpy.conversion.FeatureClassToShapefile(
            Input_Features=pretty_polygons,
            Output_Folder=r"C:\Users\PUMIGATOR\Desktop\temp"
            )

        # APPLYING ADJUSTMENTS FOR VISIBILITY AND SURVEYOR NUMBER
        adjusted_polygons = "in_memory/adjusted_polygons"
        arcpy.management.CopyFeatures(pretty_polygons, adjusted_polygons)
        if make_adjustments == 'true':
            arcpy.AddMessage(f"Adjusting '{pretty_polygons}' for visibility and surveyor number")
            adjusted_polygons = adjuster(multiplier_SURF, adjusted_polygons, surveyor_number_field, visibility_field, varlist)

        # CREATING A SUMMARY TABLE
        create_sum_table(precision, varlist, pretty_polygons, adjusted_polygons)

        # SAVING THE FINAL OUTPUT
        arcpy.management.CopyFeatures(adjusted_polygons, out_feature_class)
        arcpy.AddMessage(f"'{out_feature_class}' has been saved to as a layer in current database.")

        # CHECKING FOR ZERO VALUES IN THE FINAL OUTPUT
        zeroes_checker(out_feature_class, joinFieldCU, pottery_total_field, flints_total_field)

        # SAVING THE PROJECT
        aprx.save()

        # ADDING THE OUTPUT LAYER TO THE MAP IF REQUESTED
        if add_to_map.lower() == "true":
            arcpy.AddMessage(f"Adding the output '{out_feature_class}' to the current map")
            map_doc.addDataFromPath(os.path.join(arcpy.env.workspace, out_feature_class))
            arcpy.AddMessage(f"Layer '{out_feature_class}' has been added to the map.")

        # CLEANING UP IN-MEMORY WORKSPACE
        arcpy.management.Delete("in_memory")

    except Exception as e:
        arcpy.AddError(f"An error occurred: {str(e)}")
        raise

if __name__ == "__main__":
    main()