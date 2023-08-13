// File: portfolio/tools.js
// Version: v1.0
// Date: 2022-05-13
// Authors: Asifat Haruna Olawale
//email: olawaleasifat@gmail.com

//_____________________________________________________________________________________________

//for batch export into a drive on GEE

// Define the function for batch export
exports.batchExportToTiff = function(collection, exportOptions, suffix) {
  // Iterate over each image in the collection
  collection.aggregate_array('system:index').evaluate(function(indexes) {
      indexes.forEach(function(index) {
      //var image=collection.filterMetadata('system:index', 'equals', index).first();
      var image=collection.filterMetadata('system:index', 'equals', index).median();
      var name = index;
      var fileName = name + suffix;
      var crs = exportOptions.crs || 'EPSG:4326' //'EPSG:32632'
      var scale=exportOptions.scale|| 30
      // Set the export parameters
      var exportParams = {
        image: image.id(),
        description: fileName,
        folder: exportOptions.folder,
        scale: scale,
        crs:  crs,//exportOptions.crs,
        fileFormat: 'GeoTIFF',
        maxPixels: exportOptions.maxPixels //|| 1e13
        // Add more export options here
      };
      
      // Export the image as a GeoTIFF file
      Export.image.toDrive(exportParams);
    });
  });
}

