# CombinationTool

## Idea:
Tool to constrain the parameters of physics models using Bayesian inference by combining measurements of (different) observables. Especially suitable for EFT models. Bayesian inference is performed using the <a target="_blank" href="https://github.com/bat/BAT.jl">*Bayesian Analysis Toolkit - BAT.jl*</a>.

The likelihood used in the CombinationTool is based on the
<a target="_blank" href="https://link.springer.com/article/10.1140/epjc/s10052-016-4280-9">EFT*fitter* </a>

*  **Assumption**: Measurements of physical quantities are approximately gaussian. This allows to combine the measurements using the following likelihood:

<img src="http://latex2png.com/pngs/9fcbfbebe11bc63f39c64ac44e6bf790.png" width="450" height="55" />

with the **covariance matrix**:

 <img src="http://latex2png.com/pngs/5ceee38a8fa4eefd3b91ec36ce075746.png" width="350" height="55" />


## Input needed:
* **Model:** Predicting physical observables as a function of the model parameters
* **Measurements:** Measured values of the observables, including (various types of) uncertainties
* **Correlations:** Correlation matrices for all types of uncertainties


## Current structure of the code:
* /
  * *runCombinationTool.jl:* Main routine for performing a run of the CombinationTool
  * *inputs.jl:* Input file. Providing observables, measurements, uncertainties and correlations.

* /src
  * *CombinationTool.jl:* module definition

  * *datahandling.jl:* type definitions for observables, measurements, uncertainties and correlations. Functions to add these objects to arrays.

  * *combination.jl:* Calculation of the covariance matrix and the combination likelihood

  * *callingBAT.jl:* passing the combination likelihood to BAT.jl

* /test:
   * Unit tests for parts of the code
