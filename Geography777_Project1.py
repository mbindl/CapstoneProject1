## Geography777_Project1.py

## Created as part of the requreiments for the Master's in Cartography, GIS, and Web Map Programming.
## University of Wisconsin - Masdison - Geography 777 - Capstone - Project 1
## By: Mason Bindl
## Spring 2020

## This script investigates the relationship between meassured nitrate values in wells and cancer rates.
## The analysis interpolates values of nitrates using the Inverse Distance Weighted method.
## The analysis workflow then aggregates cancer rates meassured at the census tract or county level and
## the interpolated nitrate values to a generated tesselation feature class using Zonal Statistics.
## Oridnary Least Squares is used to meassure the correlation between the dependent variable (cancer rates)
## and the independent variable (nitrate prediction). The script then moves the residual and estimated values to
## attibute fields setup in the tesselation feature class.

## --------------------------------------------------------------------------------------------##

# Import system modules
import os
import arcpy
from arcpy import env
from arcpy.sa import *
# Set Workspace environment
env.workspace = "//ARCPROD//c$//GISData//Geog777_Project_1.gdb"
# Set Overwrite option
env.OverwriteOutput = True

## --------------------------------------------------------------------------------------------##

# function to transfer attribute values from table to feature class
def fieldJoinCalc(updateFC, updateFieldsList, sourceFC, sourceFieldsList):
    from time import strftime
    print("Started data transfer: " + strftime("%Y-%m-%d %H:%M:%S"))

    # Use list comprehension to build a dictionary from a da SearchCursor
    valueDict = {r[0]: (r[1:]) for r in arcpy.da.SearchCursor(sourceFC, sourceFieldsList)}

    with arcpy.da.UpdateCursor(updateFC, updateFieldsList) as updateRows:
        for updateRow in updateRows:
            # store the Join value of the row being updated in a keyValue variable
            keyValue = updateRow[0]
            # verify that the keyValue is in the Dictionary
            if keyValue in valueDict:
                # transfer the value stored under the keyValue from the dictionary to the updated field.
                updateRow[1] = valueDict[keyValue][0]
                updateRows.updateRow(updateRow)
    del valueDict
    print("Finished data transfer: " + strftime("%Y-%m-%d %H:%M:%S"))

## --------------------------------------------------------------------------------------------##

## Generate Tesselation

# delete feature class if it already exists
if arcpy.Exists("hexbin"):
    arcpy.Delete_management("hexbin")

tessellation_extent = arcpy.Extent(-10340405, 5234953, -9656979, 5992793)
spatial_ref = arcpy.SpatialReference(3857)
arcpy.GenerateTessellation_management("hexbin",
                                      tessellation_extent, "HEXAGON",
                                      "500 SquareMiles", spatial_ref)

# Set local variables
inFeatures = "hexbin"
fieldName1 = "ID"
fieldName2 = "CancerRate"
fieldName3 = "NitrateRate"
fieldName4 = "Residual"
fieldName5 = "Estimated"
fieldName6 = "StdResidual"

# Execute AddField twice for two new fields
arcpy.AddField_management(inFeatures, fieldName1, "SHORT")
arcpy.AddField_management(inFeatures, fieldName2, "DOUBLE")
arcpy.AddField_management(inFeatures, fieldName3, "DOUBLE")
arcpy.AddField_management(inFeatures, fieldName4, "DOUBLE")
arcpy.AddField_management(inFeatures, fieldName5, "DOUBLE")
arcpy.AddField_management(inFeatures, fieldName5, "DOUBLE")

# create unqique integer id field for OLS to use
updateFieldsList = ['OBJECTID', 'ID']
with arcpy.da.UpdateCursor(inFeatures, updateFieldsList) as cursor:
    for row in cursor:
        row[1] = row[0]
        cursor.updateRow(row)

## --------------------------------------------------------------------------------------------##

# Interpolate a series of nitrate well values
#  saved as point features onto a rectangular
#   raster using Inverse Distance Weighting (IDW).
# Requirements: Spatial Analyst Extension

# Set local variables
outRaster = os.path.join("in_memory", "idw_nitrate")
inPointFeatures = "Well_Nitrate"
zField = "nitr_con"
cellSize = 0.01
k = arcpy.GetParameter(0)
searchRadius = RadiusVariable('', 12)

# delete output if it exists
if arcpy.Exists(outRaster):
    arcpy.Delete_management(outRaster)

# Execute IDW
outIDW = Idw(inPointFeatures, zField, cellSize, k, searchRadius)

# # Save the output
# outIDW.save(outRaster)

## --------------------------------------------------------------------------------------------##

## ZonalStatisticsAsTable - Aggregate Nitrate Values to Hexbins
# Set local variables
inZoneData = "hexbin"
zoneField = "GRID_ID"
inValueRaster = outIDW
outTable = "zonalstat_nitrate"

if arcpy.Exists(outTable):
    arcpy.Delete_management(outTable)

# Execute ZonalStatisticsAsTable
outZSaT = ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster,
                                 outTable, "NODATA", "MEAN")

# transfer attributes to hexbins
fieldJoinCalc('hexbin', ['GRID_ID', 'NitrateRate'], 'zonalstat_nitrate', ['GRID_ID', 'MEAN'])
print ("The 'Nitrate Rate' field in the hexbin data has been updated")

## --------------------------------------------------------------------------------------------##

## Convert Cancer Rate feature class to raster

# Delete output raster if it exists
if arcpy.Exists("cancer_rate"):
    arcpy.Delete_management("cancer_rate")
# Set local variables
inFeature = "CancerRate_CensusTract"
outRaster = os.path.join("in_memory", "cancer_rate")
# cellSize = 0.01
field = "GEOID10"

# Execute FeatureToRaster
arcpy.FeatureToRaster_conversion(inFeature, field, outRaster)

## --------------------------------------------------------------------------------------------##

## ZonalStatisticsAsTable - Aggregate Cancer Rates to Hexbins
# Set local variables
inZoneData = "hexbin"
zoneField = "GRID_ID"
inValueRaster = "cancer_rate"
outTable = "zonalstat_tract_cancer"

if arcpy.Exists("zonalstat_tract_cancer"):
    arcpy.Delete_management("zonalstat_tract_cancer")

# Execute ZonalStatisticsAsTable
outZSaT = ZonalStatisticsAsTable(inZoneData, zoneField, inValueRaster,
                                 outTable, "NODATA", "MEAN")
# transfer attributes to hexbins
fieldJoinCalc('hexbin', ['GRID_ID', 'CancerRate'], 'zonalstat_tract_cancer', ['GRID_ID', 'MEAN'])
print ("The 'Cancer Rate' field in the hexbin data has been updated")

## --------------------------------------------------------------------------------------------##

## Perform Oridanry Least Squares Regression

olsCoefTab = 'olsCoefTab'
olsDiagTab = 'olsDiagTab'
ols = "OLS"
try:
    if arcpy.Exists(ols):
        arcpy.Delete_management(ols)
    if arcpy.Exists(olsCoefTab):
        arcpy.Delete_management(olsCoefTab)
    if arcpy.Exists(olsDiagTab):
        arcpy.Delete_management(olsDiagTab)

    arcpy.OrdinaryLeastSquares_stats("hexbin", "ID", ols,
                                     "CancerRate", "NitrateRate",
                                     olsCoefTab, olsDiagTab)
    # transfer attributes to Parcel Layer
    fieldJoinCalc('hexbin', ['ID', 'Residual'], 'OLS', ['ID', 'Residual'])
    print("The 'Residual' field in the hexbin data has been updated")

    # transfer attributes to Parcel Layer
    fieldJoinCalc('hexbin', ['ID', 'Estimated'], 'OLS', ['ID', 'Estimated'])
    print("The 'Estimated' field in the hexbin data has been updated")
except:
    # If an error occurred when running the tool, print out the error message.
    print(arcpy.GetMessages())

## --------------------------------------------------------------------------------------------##

