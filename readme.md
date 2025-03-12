Habitat suitability under climate change
Our changing climate is changing where key grassland species can live, and grassland management and restoration practices will need to take this into account.

In this coding challenge, you will create a habitat suitability model for a species of your choice that lives in the continental United States (CONUS). We have this limitation because the downscaled climate data we suggest, the MACAv2 dataset, is only available in the CONUS – if you find other downscaled climate data at an appropriate resolution you are welcome to choose a different study area. If you don’t have anything in mind, you can take a look at Sorghastrum nutans, a grass native to North America. In the past 50 years, its range has moved northward.

Your suitability assessment will be based on combining multiple data layers related to soil, topography, and climate. You will also need to create a modular, reproducible, workflow using functions and loops. To do this effectively, we recommend planning your code out in advance using a technique such as pseudocode outline or a flow diagram. We recommend planning each of the blocks below out into multiple steps. It is unnecessary to write a step for every line of code unles you find that useful. As a rule of thumb, aim for steps that cover the major structures of your code in 2-5 line chunks.

STEP 1: STUDY OVERVIEW
Before you begin coding, you will need to design your study.

Reflect and Respond
What question do you hope to answer about potential future changes in habitat suitability?

YOUR QUESTION HERE

Species
Try It
Select the species you want to study, and research it’s habitat parameters in scientific studies or other reliable sources. You will want to look for reviews or overviews of the data, since an individual study may not have the breadth needed for this purpose. In the US, the National Resource Conservation Service can have helpful fact sheets about different species. University Extension programs are also good resources for summaries.

Based on your research, select soil, topographic, and climate variables that you can use to determine if a particular location and time period is a suitable habitat for your species.

Reflect and Respond
Write a description of your species. What habitat is it found in? What is its geographic range? What, if any, are conservation threats to the species? What data will shed the most light on habitat suitability for this species?

YOUR SPECIES DESCRIPTION HERE

Sites
Try It
Select at least two site to study, such as two of the U.S. National Grasslands. You can download the USFS National Grassland Units and select your study sites. Generate a site map for each location.

When selecting your sites, you might want to look for places that are marginally habitable for this species, since those locations will be most likely to show changes due to climate.

Reflect and Respond
Write a site description for each of your sites, or for all of your sites as a group if you have chosen a large number of linked sites. What differences or trends do you expect to see among your sites?

YOUR SITE DESCRIPTION HERE

Time periods
In general when studying climate, we are interested in climate normals, which are typically calculated from 30 years of data so that they reflect the climate as a whole and not a single year which may be anomalous. So if you are interested in the climate around 2050, download at least data from 2035-2065.

Reflect and Respond
Select at least two 30-year time periods to compare, such as historical and 30 years into the future. These time periods should help you to answer your scientific question.

YOUR TIME PERIODS HERE

Climate models
There is a great deal of uncertainty among the many global climate models available. One way to work with the variety is by using an ensemble of models to try to capture that uncertainty. This also gives you an idea of the range of possible values you might expect! To be most efficient with your time and computing resources, you can use a subset of all the climate models available to you. However, for each scenario, you should attempt to include models that are:

Warm and wet
Warm and dry
Cold and wet
Cold and dry
for each of your sites.

To figure out which climate models to use, you will need to access summary data near your sites for each of the climate models. You can do this using the Climate Futures Toolbox Future Climate Scatter tool. There is no need to write code to select your climate models, since this choice is something that requires your judgement and only needs to be done once.

If your question requires it, you can also choose to include multiple climate variables, such as temperature and precipitation, and/or multiple emissions scenarios, such as RCP4.5 and RCP8.5.


Try It
Choose at least 4 climate models that cover the range of possible future climate variability at your sites. How did you choose?
LIST THE CLIMATE MODELS YOU SELECTED HERE AND CITE THE CLIMATE TOOLBOX

STEP 2: DATA ACCESS
Soil data
The POLARIS dataset is a convenient way to uniformly access a variety of soil parameters such as pH and percent clay in the US. It is available for a range of depths (in cm) and split into 1x1 degree tiles.

Try It
Write a function with a numpy-style docstring that will download POLARIS data for a particular location, soil parameter, and soil depth. Your function should account for the situation where your site boundary crosses over multiple tiles, and merge the necessary data together.

Then, use loops to download and organize the rasters you will need to complete this section. Include soil parameters that will help you to answer your scientific question. We recommend using a soil depth that best corresponds with the rooting depth of your species.


# Download soil data
     
Topographic data
One way to access reliable elevation data is from the SRTM dataset, available through the earthaccess API.

Try It
Write a function with a numpy-style docstring that will download SRTM elevation data for a particular location and calculate any additional topographic variables you need such as slope or aspect.

Then, use loops to download and organize the rasters you will need to complete this section. Include topographic parameters that will help you to answer your scientific question.

Warning

Be careful when computing the slope from elevation that the units of elevation match the projection units (e.g. meters and meters, not meters and degrees). You will need to project the SRTM data to complete this calculation correctly.


# Download soil data
     
Climate model data
You can use MACAv2 data for historical and future climate data. Be sure to compare at least two 30-year time periods (e.g. historical vs. 10 years in the future) for at least four of the CMIP models. Overall, you should be downloading at least 8 climate rasters for each of your sites, for a total of 16. You will need to use loops and/or functions to do this cleanly!.

Try It
Write a function with a numpy-style docstring that will download MACAv2 data for a particular climate model, emissions scenario, spatial domain, and time frame. Then, use loops to download and organize the 16+ rasters you will need to complete this section. The MACAv2 dataset is accessible from their Thredds server. Include an arrangement of sites, models, emissions scenarios, and time periods that will help you to answer your scientific question.


# Download climate data
     
Reflect and Respond
Make sure to include a description of the climate data and how you selected your models. Include a citation of the MACAv2 data

YOUR CLIMATE DATA DESCRIPTION AND CITATIONS HERE

STEP 3: HARMONIZE DATA
Try It
Make sure that the grids for all your data match each other. Check out the ds.rio.reproject_match() method from rioxarray. Make sure to use the data source that has the highest resolution as a template!

Warning

If you are reprojecting data as you need to here, the order of operations is important! Recall that reprojecting will typically tilt your data, leaving narrow sections of the data at the edge blank. However, to reproject efficiently it is best for the raster to be as small as possible before performing the operation. We recommend the following process:

1. Crop the data, leaving a buffer around the final boundary
2. Reproject to match the template grid (this will also crop any leftovers off the image)

# Download soil data
     
STEP 4: DEVELOP A FUZZY LOGIC MODEL
A fuzzy logic model is one that is built on expert knowledge rather than training data. You may wish to use the scikit-fuzzy library, which includes many utilities for building this sort of model. In particular, it contains a number of membership functions which can convert your data into values from 0 to 1 using information such as, for example, the maximum, minimum, and optimal values for soil pH.

Try It
To train a fuzzy logic habitat suitability model:

1. Research S. nutans, and find out what optimal values are for each variable you are using (e.g. soil pH, slope, and current climatological annual precipitation). 
2. For each **digital number** in each raster, assign a **continuous** value from 0 to 1 for how close that grid square is to the optimum range (1=optimal, 0=incompatible). 
3. Combine your layers by multiplying them together. This will give you a single suitability number for each square.
4. Optionally, you may apply a suitability threshold to make the most suitable areas pop on your map.
Tip

If you use mathematical operators on a raster in Python, it will automatically perform the operation for every number in the raster. This type of operation is known as a vectorized function. DO NOT DO THIS WITH A LOOP!. A vectorized function that operates on the whole array at once will be much easier and faster.


# Create fuzzy logic suitability model
     
STEP 5: PRESENT YOUR RESULTS
Try It
Generate some plots that show your key findings. Don’t forget to interpret your plots!


# Create plots
     
YOUR PLOT INTERPRETATION HERE

