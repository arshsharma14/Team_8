## October 2, 2024
### Agenda
- Decide on common variables to use between the three datasets (elevation, pH, etc.)
- Ask what TOC is in atacama soil dataset
- What is the difference between SOC and TOC?
- Potential Questions:
    - How do elevation and temperature impact the bacterial community composition and functional traits in soil vs. sand environments?
    - How do moisture content and relative humidity shape microbial functional diversity in soil vs. sand environments?
    - How does the difference in pH between soil and sand samples affect the diversity and functional capabilities of microbial communities?
    - Does the presence or absence of vegetation influence the functional diversity of microbial communities in soil vs. sand?
    - How does nitrogen and carbon composition of soil vs sand microbiome associate with the ability to grow/rotate crops?
    - How does microbial diversity vary with changes in pH across different elevations in arid (Atacama), forest, and wetland ecosystems?
### Meeting Minutes
- TOC & SOC : carbon content of soil
- Research Question: How does microbial diversity vary with changes in pH across different elevations in arid (Atacama), forest, and wetland ecosystems?
- pH bin -> acidic, neutal, basic,
- Difference between datasets overall?
- regions to compare in dataset
- lot of data wranggling before
- table -> qualities of dataset, sample each dataset, variable regions
      - combine table level
- 3 new metadata -> catagorical
- metadata seperate and combine
- wrangle data -> keep or delete? -> don't delete to keep exploring
- make sure all common headers named same -> so can combine via datasets
- carbon, nitrogen ,mositure, elevation
- each datasets -> confounding variables -> keep from one region
- maifest & metadata both will have to trimmed
- don't delete -> don't correct for any confounding variakes -> atacoma soil
- don’t remove
- wrangling -> metadata combine
- filter regions, 
- 200, 100, 100
- backup: functional analysis

data wrangling

Determine paired end reads vs single end reads for each dataset
How many samples?
What is the variable region?
	If they are all the same, we can pool them
	If they're not, we will have to make a table
Have to control for certain variables
For soil and wetlands, we have to filter the metadata and manifest
Maybe some species react differently to pH in different environments (to consider) 
