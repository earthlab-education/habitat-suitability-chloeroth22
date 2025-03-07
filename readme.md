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
