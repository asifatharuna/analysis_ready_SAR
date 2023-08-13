var senraw = ee.ImageCollection("COPERNICUS/S1_GRD")
var boundary = ee.FeatureCollection("users/olawaleasifat/extended_boundary_new")

// File: portfolio/sentinel 1 preprocessing.js
// Version: v1.0
// Date: 2022-05-15
// Authors: Asifat Haruna Olawale
//email: olawaleasifat@gmail.com

//Copyright © 2022 Asifat Haruna. All rights reserved. 
//This work and its accompanying resources, including but not limited to code, documentation, and any associated materials, are protected by copyright laws and international treaties.
// Unauthorized use, reproduction, distribution, or modification of this work, in whole or in part, is strictly prohibited without prior written permission from the copyright holder. 
//Any unauthorized use may result in legal action and be subject to applicable penalties and damages. For inquiries or permissions, please contact Asifat Haruna at olawaleasifat@gmail.com
            
                  //****************************
                        // PREPROCESSING
                //****************************
           
  //Preprocessing SAR (Synthetic Aperture Radar) imagery is essential to enhance data quality and extract meaningful information. 
  //It involves 1) calibration, 2) speckle filtering, 3) geometric correction, 4) radiometric normalization, 5) and terrain correction. 
  //Calibration removes system-specific effects, while speckle filtering reduces noise. Geometric correction corrects distortions caused by sensor position and platform motion. 
  //Radiometric normalization ensures consistent pixel values, and terrain correction accounts for ground elevation variations. 
  //These preprocessing steps improve the accuracy, reliability, and usability of SAR imagery for various applications, including land cover mapping and monitoring natural disasters.
  //This collection includes the S1 Ground Range Detected (GRD) scenes, processed using the Sentinel-1 Toolbox to generate a calibrated, ortho-corrected product. 

//_____________________________________________________________________________________________
  //The analysis was conducted using data from Hunsrück-Hochwald National Park in Germany. 
  //The boundary of the park is represented by a shape file named "boundary." 
  //If you wish to use the same boundary, you can access it through the following URL link:  
  //https://code.earthengine.google.com/?asset=users/olawaleasifat/extended_boundary_new
  // After importing it into your code editor, remember to rename it as "boundary."
  //If you have a different area of interest, you have the option to replace the input boundary shape file 
  //or use the drawing tools provided on GEE to delineate your desired area.
//_____________________________________________________________________________________________________


//Modify the date range to fit your study period
var s1_db= senraw.filterBounds(boundary)
        .filterDate("2022-11-28","2023-04-14")// 
        .filterMetadata("relativeOrbitNumber_start", 'equals' , 15)
        .filter(ee.Filter.eq('transmitterReceiverPolarisation', ['VV','VH'])) 
        .filter(ee.Filter.eq('resolution_meters', 10)) 
        .filter(ee.Filter.eq('instrumentMode', 'IW')) 
        .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING')) 
        .sort('system:time_start').select("VV","VH").aside(print)
    

// format the start date 
var format_date= function(image){
         var date_fmt= ee.Date(image.get("system:time_start")).format("YYYY-MM-dd")
         return image.set("date", date_fmt)
}



//Preparing the date format for facilitating spatial aggregation or  mosaic of contiguous tiles.
s1_db= s1_db.map(format_date).aside(print)

var dates= ee.List(s1_db.aggregate_array("date"))
var unique_dates= dates.distinct().aside(print)

var obs_dates_pairs=unique_dates.slice(0,-1).zip(unique_dates.slice(1)).aside(print)

var observed_daily_Images= obs_dates_pairs.map(function(d) {
    
    var start_date=ee.List(d).get(0)
    var end_date = ee.List(d).get(1)
    var filtered=  s1_db.filterDate(start_date, end_date)
    
    var aggre_t0=filtered.aggregate_array("system:time_start").get(-1)
    //var aggre_t1 = filtered.aggregate_array("system:time_end").get(-1)
    var aggre_index= filtered.aggregate_array("system:index").get(-1)
    
    
    var dict =ee.Dictionary({
                       "system:time_start":aggre_t0, "system:index":aggre_index
                      })
   
   //var monthly = filtered.median().set(_dict)
    var observed_daily = filtered.mosaic().set({
                     "system:time_start":aggre_t0, "system:index":aggre_index
                      })  
   
  return observed_daily 
   
}).aside(print)


// return back the variable name
s1_db=ee.ImageCollection(observed_daily_Images)




// Convert backscatter from dB to linear. 
// It is recommended to perform the analysis using linear units.
var db_to_lin = function(image) {
  var bandNames = image.bandNames().remove('angle');
  var lin = ee.Image.constant(10).pow(image.select(bandNames).divide(10)).rename(bandNames)
  return image.addBands(lin, null, true)
};


// Convert backscatter from linear to dB
// For converting the  backscatter to dB before exporting the analysis-ready product 
var lin_to_db = function(image) {
  var bandNames = image.bandNames().remove('angle');
  var db = ee.Image.constant(10).multiply(image.select(bandNames).log10()).rename(bandNames)
  return image.addBands(db, null, true)
};

//---------------------------------------------------------------------------//
// GAMMA MAP filter 
//---------------------------------------------------------------------------//
/** Gamma Maximum a-posterior Filter applied to one image. It is implemented as described in 
Lopes A., Nezry, E., Touzi, R., and Laur, H., 1990.  Maximum A Posteriori Speckle Filtering and First Order texture Models in SAR Images.  
International  Geoscience  and  Remote  Sensing  Symposium (IGARSS).  */

// Source:Mullissa, A.; Vollrath, A.; Odongo-Braun, C.; Slagter, B.; Balling, J.; Gou, Y.; Gorelick, N.; Reiche, J. Sentinel-1 SAR Backscatter Analysis Ready Data Preparation in Google Earth Engine. Remote Sens. 2021, 13, 1954. 
// https://doi.org/10.3390/rs13101954

var gammamap =  function(image) { 
        var KERNEL_SIZE=7
        var enl = 5;
        var bandNames = image.bandNames().remove('angle');
        //Neighbourhood stats
        var reducers = ee.Reducer.mean().combine({
                      reducer2: ee.Reducer.stdDev(),
                      sharedInputs: true
                      });
        var stats = image.select(bandNames).reduceNeighborhood({reducer: reducers,kernel: ee.Kernel.square(KERNEL_SIZE/2,'pixels'), optimization: 'window'})
        var meanBand = bandNames.map(function(bandName){return ee.String(bandName).cat('_mean')});
        var stdDevBand = bandNames.map(function(bandName){return ee.String(bandName).cat('_stdDev')});
        
        var z = stats.select(meanBand);
        var sigz = stats.select(stdDevBand);
        
        // local observed coefficient of variation
        var ci = sigz.divide(z);
        // noise coefficient of variation (or noise sigma)
        var cu = 1.0/Math.sqrt(enl);
        // threshold for the observed coefficient of variation
        var cmax = Math.sqrt(2.0) * cu
  
        cu = ee.Image.constant(cu);
        cmax = ee.Image.constant(cmax);
        var enlImg = ee.Image.constant(enl);
        var oneImg = ee.Image.constant(1);
        var twoImg = ee.Image.constant(2);
  
        var alpha = oneImg.add(cu.pow(2)).divide(ci.pow(2).subtract(cu.pow(2)));

        //Implements the Gamma MAP filter described in equation 11 in Lopez et al. 1990
        var q = image.select(bandNames).expression("z**2 * (z * alpha - enl - 1)**2 + 4 * alpha * enl * b() * z", {z: z, alpha: alpha,enl: enl})
        var rHat = z.multiply(alpha.subtract(enlImg).subtract(oneImg)).add(q.sqrt()).divide(twoImg.multiply(alpha));
  
        //if ci <= cu then its a homogenous region ->> boxcar filter
        var zHat = (z.updateMask(ci.lte(cu))).rename(bandNames)
        //if cmax > ci > cu then its a textured medium ->> apply Gamma MAP filter
        rHat = (rHat.updateMask(ci.gt(cu)).updateMask(ci.lt(cmax))).rename(bandNames)
        //if ci>=cmax then its strong signal ->> retain
        var x = image.select(bandNames).updateMask(ci.gte(cmax)).rename(bandNames)
  
        // Merge
        var output = ee.ImageCollection([zHat,rHat,x]).sum();
        return image.addBands(output, null, true);
  }   

//---------------------------------------------------------------------------//
//convert to linear scale->apply the speekle filter function
var s1_ln=s1_db.map(db_to_lin).map(gammamap).aside(print)

                    
//convert back to db for export 
s1_db=s1_ln.map(lin_to_db).aside(print)

//_______________________________________________________________________________________________________

                //********************************************************
                                  //Visualization
                //********************************************************
//_______________________________________________________________________________________________________


//A time series plot is generated for a random point within the study area, 
//which can be replaced with a point of interest chosen by the user.
var roi=ee.Geometry.Point([ 6.945508,49.782008])  

//Map.centerObject(boundary);
Map.centerObject(roi);

//setting visualization parameters
var visparam = {bands:['VV','VH','VV'],min: [-20, -25, -25],max: [0, -5, -5]}

Map.addLayer( s1_db.first(),visparam , 'sentinel_1')
Map.addLayer(boundary,{color:"red"},"boundary");


// Make another chart.
var chart = ui.Chart.image.series({
  imageCollection: ee.ImageCollection(s1_ln).select("VV"),
  
  region:roi, //again, here the ROI is a single point
  reducer: ee.Reducer.first(), 
  scale: 10
});

// Define custom options for the chart. See:
// https://developers.google.com/chart/interactive/docs/reference
var options = { 
  title: 'Gamma_ln over time', 
  hAxis: { title: 'time' },
  vAxis: { title: 'Backscatter_VV' },
  series: {
    0: { color: 'blue' }
  }
};

// Set the options of the chart and print it.
chart = chart.setOptions(options);
print(chart); // Print out in the console



// for easy comparsion of linear and db scale, make another chart for backscatter in db scale 
var chart = ui.Chart.image.series({
  imageCollection: ee.ImageCollection(s1_db).select("VV"),
  
  region:roi, //again, here the ROI is a single point
  reducer: ee.Reducer.first(), 
  scale: 10
});

// Define custom options for the chart. See:
// https://developers.google.com/chart/interactive/docs/reference
var options = { 
  title: 'Gamma_db over time', 
  hAxis: { title: 'time' },
  vAxis: { title: 'Backscatter_VV' },
  series: {
    0: { color: 'blue' }
  }
};

// Set the options of the chart and print it.
chart = chart.setOptions(options);
print(chart); // Print out in the console


//_______________________________________________________________________________________________________

                //********************************************************
                          //Export prepossed products
                //********************************************************
//_______________________________________________________________________________________________________

//import a function from the tools module and set
// export options and suffix to the outname if necessary otherwise use '' for the suffix

var utilities = require('users/olawaleasifat/portfolio:tools')

var exportOptions = {
  folder: 'GEE_EXPORT_SAR',
  scale: 10,
 crs: 'EPSG:25832',
  maxPixels: 1e13
};

var suffix = "_Gamma_db_10m"; // Add your desired suffix here

// The preprocess products can be exported in batch into a google cloud storage(drive)
//where it could be mounted and used for further analysis
//uncomment if you want to export the product
            //||||||
            //vvvvvv
//utilities.batchExportToTiff(s1_db, exportOptions,suffix)







