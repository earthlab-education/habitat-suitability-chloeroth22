## Habitat suitability under climate change
### Study Overview
#### Question
In this coding challenge my goal was to create a habitat suitability model for Cylindropuntia fulgida, or Jumping Cholla, a cactus native to the Southwest United States and Northern Mexico. I aimed to answer the question of will Jumping Cholla's range be constrained or expanded by the changing temperature and precipitation caused by climate change. 

#### Species
I decided to do my study on Cylindropuntia fulgida, or Jumping Cholla. Jumping Cholla are a type of cactus native to Arizona and Northern Mexico, and are popular in the Sonoran desert. They are known as Jumping Chollas because of how easy it is for the stems to detach when brushed, which makes it seem like the stems are nearly jumping off the cactus. This ease of detachment lets the cholla attach itself to desert animals and seeds are dispersed.

<img 
  src="img/Opuntia_fulgida_1_-_Desert_Botanical_Garden.jpg" 
  alt="Jumping Cholla" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

<img 
  src="img/Opuntia_fulgida_range_map.jpg" 
  alt="Jumping Cholla Distribution Map" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

Since cholla's exist in already harsh and unforgiving environments, their population is stable. I wondered if the highest climate extremes would be able to push even a stable plant's suitability zones. Information on the optimal soil, topographic, and climate variables was difficult to ascertain so I mixed information on suitable variables for both jumping cholla specifically and cactus in general. 

##### Soil Variables
The best soil conditions for jumping cholla are ones that are fertile, gravel-rich, water-permeable soil with a pH of between 6 and 7. Generally cactus will do well in soil with a pH of between 5 and 7, as they prefer slight acidity. Since the POLARIS data set does not include percent gravel in their data I collected the percent clay and determined if the clay percentage was under 35% then it would be more optimal and more likely to be a well drained soil. I also used the mean clay % and pH level at 5-15 cm below the ground to account for cactus shallow roots.

<img 
  src="img/HSG_USDA_overlap.png" 
  alt="Soil Triangle" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

##### Topographic Variables
Jumping chollas grow from 980 to 3,280 feet in elevation and likely can withstand steep slopes, so I selected 50% slope as my maximum slope with under 30% being ideal.

##### Climate Variables
For temperature I determined the ideal range to be between 15 - 38 ℃ and the temperature tolerance to be between -5 - 43 ℃. For precipitation I was unable to find research on specifically Jumping Cholla but found that saguaro cacti survive in an ideal range 180 to 425 mm and can handle up to 500 mm of rainfall per year and down to 25 mm of rainfall annually. 

The major limiting factors for jumping cholla are cold winters, dry summers, and too much rainfall with low soil permeability. 

#### Sites
To select sites I downloaded GBIF data and compared it to the Preservation Areas in the United States to determine which areas had the highest number of Jumping Cholla recorded. 

<img 
  src="img/GBIF_oc.png" 
  alt="GBIF" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

I was able to determine that Organ Pipe Cactus National Monument and Tucson Mountain Park had the highest number of Jumping Cholla recorded.

<img 
  src="img/orpi_site.png" 
  alt="ORPI Site" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

Organ Pipe Cactus National Monument is located in the very Southern part of Arizona and shares a border with Mexico. It is 517 sq mi big, most of which is a wilderness area. It is the only place in the U.S. where senita and organ pipe cactus grow wild, and is home to many different cactus species. 

<img 
  src="img/tmp_site.png" 
  alt="TMP Site" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">

Tucson Mountain Park is a 20,000 acre park located to the West of Tucson, Arizona, and shares its Northern border with Saguaro National Park West. 

#### Time periods
I wanted to compare the extremes on both ends of the spectrum, so I chose to look at a historical time period of 1950-1979 and a late century time period of 2064-2096. 

### Climate models
Due to the uncertainty among global climate models, four different climate models were chosen to include scenarios that were warm and wet, warm and dry, cold and wet, and cold and dry. To select these I used the Climate Futures Toolbox Future Climate Scatter tool. 

<img 
  src="img/scatter.png" 
  alt="Climate Scatter" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;"> 

The climate models I decided to use were:

Cold and wet: MRI CGCM3

Hot and dry: IPSL CM5A MR

Hot and wet: CanESM2

Cold and dry: inmcm4

##### Citations
Climate Toolbox. (n.d.). Future climate scatter. Climate Toolbox. https://climatetoolbox.org/tool/Future-Climate-Scatter

Helmy, O. (2021). Carnegiea gigantea, saguaro. Fire Effects Information System. U.S. Department of Agriculture, Forest Service, Rocky Mountain Research Station. https://www.fs.usda.gov/database/feis/plants/cactus/cargig/all.pdf

Iowa State University Extension and Outreach. (n.d.). Gardening on slopes and hillsides. Iowa State University Extension and Outreach. https://yardandgarden.extension.iastate.edu/how-to/gardening-slopes-and-hillsides

PictureThis. (n.d.). How to grow and care for jumping cholla (Cylindropuntia fulgida). PictureThis. https://www.picturethisai.com/care/Cylindropuntia_fulgida.html

Wikipedia contributors. (n.d.). Cylindropuntia fulgida. Wikipedia. https://en.wikipedia.org/wiki/Cylindropuntia_fulgida

Wikipedia contributors. (n.d.). Organ Pipe Cactus National Monument. Wikipedia. https://en.wikipedia.org/wiki/Organ_Pipe_Cactus_National_Monument

### Data Access
#### Soil data
<img 
  src="img/soil.png" 
  alt="POLARIS" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">
Using the POLARIS dataset I downloaded mean clay % and soil pH at a depth of 5-15 cm for each site. 
     
#### Topographic data
<img 
  src="img/elevation.png" 
  alt="Elevation" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">
Using the SRTM dataset, downloaded through the earthaccess API, I downloaded elevation data for each site. 

<img 
  src="img/slope.png" 
  alt="Slope" 
  style="max-width:100%; height:auto; display:block; margin-left:auto; margin-right:auto;">
Next I calculated the slope of the SRTM elevation dataset using the xrspatial package. 
     
#### Climate model data
I downloaded MACAv2 data for historical data (1950-1979) and higher emissions RCP8.5 future climate data (2064-2096), using the four climate models I selected using the Climate Futures Toolbox Future Climate Scatter tool. I downloaded the max temperature data for each one of my sites, models, and time periods. 

Northwest Knowledge Network. (n.d.). REACCH climate CMIP5 MACAV2 catalog. Northwest Knowledge Network. https://www.reacchpna.org/thredds/reacch_climate_CMIP5_macav2_catalog2.html

#### Develop a Fuzzy Logic Model
I used the scikit-fuzzy library trapezoidal function to build a fuzzy model based on the known optimal values. To train the habitat suitability model I determined the optimal values, then for each **digital number** in each raster, I assigned a **continuous** value from 0 to 1 for how close that grid square is to the optimum range (1=optimal, 0=incompatible). I combined my layers by multiplying them together to get a single suitability score for each grid square on my raster. 

Insert Plots here

Looking at my plots I have a lot of concerns about the quality of my data and processing. It looks like clay and slope were very compatible with the regions, but pH, temperature, and precipitation were not very compatible, even in the historic models. In the future I would like to overlay the GBIF occurrences data to get a better sense for where the species is currently occuring so I can assess the accuracy of my model better. Also, I think I should have picked a more marginal species and more marginal sites so I could have better seen the effects of changes in temperature and precipitation. Additionally I had some issues clipping so there are spots in my plots that are marked as compatible, but they don't actually contain data. There were reprojecting issues, which may be due to order of operations issues and being unable to clip the climate data as tightly as I was able to clip the soil and elevation data. Finally, I would pick a species that had more accurate and readily available optimal conditions data available. 


