## Informal write-up of broad concepts that will be used.
1. Geometric measurements in FRW Cosmology : background expansion of a homogeneous and isotropic universe (FRW universe), and the impact on measuring quantities like brightness or angles from the earth. This is important for our studies because supernova cosmology is based on the effect of the expansion history of the universe on the measured brightness of objects.
2. Introductory Bayesian Statistics : There are two basic kinds of questions we will play with: (a) Parameter Estimation: Given a data set with measurements and uncertainties of measurements,  and assuming the data comes from a model (with unknown parameters) and noise distribution, how can one make probabilitic inferences on the unknown parameters? (b) Model Comparison (also sometimes called Model Selection) If we had the same dataset, but did not know the model, and had to choose from a few models, how can we compare these models?
3. These will then be applied to models of dark energy based on simulations of future datasets. In particular, we will be interested in how the addition of a low redshift dataset like ZTF enhances our capabilities to answer the these questions for the same higher redshift dataset.

## Rough Plan for the Summer
### Section 1 : Study of basic concepts and hands on numerical exploration of the concepts using python libraries (Weeks 1-3)
- Week 1 : Setting up, discussions and learning 1. and 2. as well as goals of for  weeks 2 and 3.
- Week 2 : Numerical exlporation of geometries and FRW universes using astropy for standard models, scipy for newer models.
- Week 3 : Numerical exploration of Statistical concepts using simpler mathematical functions
### Section 2: Use existing code to compute the cosmological constraints from future surveys
- Week 4 : Model Comparison of dark energy models (Suhail has some code setup that can be used) 
- Week 5 : Use the experience of week 2 to add in new models of dark energy
### Section 3: Analyse the impact of different datasets, interpret results and organise codebase
- Week 6, 7, 8, 9 : Setup the same studies for a few combinations of datasets, and models, and finally look into the results and interpret them. It may be a good idea to clean up and re-organize the codebase a little more. This will help in our ability to run different datasets, add probes, and models in a cleaner and more reproducible manner. We can decide on the details and extent of this exercise depending on the status of the project.
- After that, the time should be used to attend the ZTF meeting, write up reports and explore next steps.
