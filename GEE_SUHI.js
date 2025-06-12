var comune_bologna = ee.FeatureCollection("projects/hi-marikadagostini/assets/shape_BO");

var L8 = ee.ImageCollection("LANDSAT/LC08/C02/T1");

// Function to mask clouds and shadows
function maskL8sr(image) {
  var dilatedCloudBitMask = (1 << 1);
  var cirrusBitMask = (1 << 2);
  var cloudsBitMask = (1 << 3);
  var cloudShadowBitMask = (1 << 4);
  
/*  Bitmask for QA_PIXEL
    Bit 1: Dilated Cloud
      0: Cloud is not dilated or no cloud
      1: Cloud dilation
    Bit 2: Cirrus
      0: Cirrus Confidence: no confidence level set or Low Confidence
      1: High confidence cirrus
    Bit 3: Cloud
      0: Cloud confidence is not high
      1: High confidence cloud
    Bit 4: Cloud Shadow
      0: Cloud Shadow Confidence is not high
      1: High confidence cloud shadow
*/


  // filter pixel using QA band
  var qa = image.select('QA_PIXEL');
  
  var mask = qa.bitwiseAnd(dilatedCloudBitMask).eq(0)
                 .and(qa.bitwiseAnd(cirrusBitMask).eq(0))
                 .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
                 .and(qa.bitwiseAnd(cloudShadowBitMask).eq(0));
  
  // store QA band
  var q2 = image.select('QA_PIXEL').rename('QA')
  return image.updateMask(mask).addBands(q2);
}

// Function to calculate TEMP and extract QA
function calcTEMP(image) {
  
  // top of atmosphere spectral radiance (TOA)
  var ML = ee.Number(image.get('RADIANCE_MULT_BAND_10'));
  var AL = ee.Number(image.get('RADIANCE_ADD_BAND_10'));
  var L_lambda = image.select('B10').multiply(ML).add(AL);
  
  // at-satellite brightness temperature [K]
  var K1 = ee.Image.constant(image.get('K1_CONSTANT_BAND_10'));
  var K2 = ee.Image.constant(image.get('K2_CONSTANT_BAND_10'));
  var Tb = K2.divide(K1.divide(L_lambda).add(1).log());
  
  // Normalized Different Vegetation Index (NDVI)
  var ML_5 = ee.Number(image.get('RADIANCE_MULT_BAND_5'));
  var AL_5 = ee.Number(image.get('RADIANCE_ADD_BAND_5'));
  var B5_scaled = image.select('B5').multiply(ML_5).add(AL_5);
  
  var ML_4 = ee.Number(image.get('RADIANCE_MULT_BAND_4'));
  var AL_4 = ee.Number(image.get('RADIANCE_ADD_BAND_4'));
  var B4_scaled = image.select('B4').multiply(ML_4).add(AL_4);
  
  var NDVI = (B5_scaled.subtract(B4_scaled)).divide(B5_scaled.add(B4_scaled));
  
  var NDVI_soil = ee.Image.constant(0.1); // NDVI value for bare soils
  var NDVI_veg = ee.Image.constant(0.65); // NDVI value for fully vegetated soils
  
  // Fractional Vegetation Cover (FVC)
  var FVC = (NDVI.subtract(NDVI_soil)).divide(NDVI_veg.subtract(NDVI_soil));
  
  // Land Surface Emissivity (LSE)
  var epsilon_soil = ee.Image.constant(0.93); // typical soil emissivity
  var epsilon_veg = ee.Image.constant(0.98); // typical vegetation emissivity
  
  var epsilon = (epsilon_soil.multiply(ee.Image.constant(1).subtract(FVC)))
          .add(epsilon_veg.multiply(FVC));
          
  var Lambda = ee.Image.constant(10.895e-6); // wavelength of the emitted radiance [μm]
  
  var s = ee.Image.constant(1.38e-23); // Boltzmann constant [J/K]
  var h = ee.Image.constant(6.626e-34); // Planck’s constant [Js]
  var c = ee.Image.constant(2.998e8); // velocity of light [m/s]
  
  //var rho = h.multiply(c).divide(s)
  var rho = ee.Image.constant(1.4388e-2); // [mK]
  
  // Estimation of Land Surface Temperature
  var LST = Tb.divide(ee.Image.constant(1)
    .add(Lambda.multiply(Tb).divide(rho)
    .multiply(epsilon.log())));

  // from K to Celsius
  var temp = LST.subtract(273.15).rename('TEMP').toFloat(); 
  var qa3 = image.select('QA').toFloat(); 

  return temp.addBands(qa3);  // Ensure only TEMP & QA_PIXEL are exported
}

// Filter images by year, apply cloud mask and temperature calculation
var filteredCollection = L8
  .filterDate('2013-05-01', '2024-09-30')
  .filterBounds(comune_bologna)
  .map(maskL8sr)
  .map(calcTEMP);

// Convert to a list for iteration
var imageList = filteredCollection.toList(filteredCollection.size());

imageList.size().evaluate(function(count) {
  for (var i = 0; i < count; i++) {
    var img = ee.Image(imageList.get(i)); // Get image from list
    var id = img.get('system:index').getInfo(); // Get image ID
    
    Export.image.toDrive({
      image: img.clip(comune_bologna),
      description: 'Bologna_Temp_QA_' + id,
      region: comune_bologna,
      scale: 30,
      folder: 'PhD_Temp_QA_Bologna',
      fileFormat: 'GeoTIFF',
      maxPixels: 1e13
    });
  }
});




