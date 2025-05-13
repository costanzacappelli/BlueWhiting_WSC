**Ocean-climate conditions drive the dynamics of a large fish population and forecast reproductive success one year before spawning**

This repository contains the data and R scripts used in the analyses for the manuscript entitled "Ocean-climate conditions drive the dynamics of a large fish population and forecast reproductive success one year before spawning.". The focus of this work is to investigate how wind stress curl (WSC), a key indicator of ocean-climate variability, influences the recruitment dynamics of blue whiting (Micromesistius poutassou) in the Northeast Atlantic. The analyses aim to identify climate-recruitment relationships and assess their potential for forecasting reproductive success up to one year in advance.

### DATA:
This folder contains all input datasets and model outputs needed to perform the analyses described in the scripts. 

Loaded in script _Preparing_WSC_dataset_script2025.R_:  
1. **windstress_NCARNCEPreanalysis_20241002/**  
Folder containing NetCDF files for monthly mean wind stress components at 10m resolution:
- uwnd = Eastward wind component
- vwnd = Northward wind component  
Source: NCAR/NCEP Reanalysis. Downloaded in October 2024; data available from January 1948 to September 2024.    


Loaded in script _BW_WSC_analyses_script2025.R_:   
2. **BW_WSC_data_19812023.xlsx**  
Time series data used for regression models.  
Columns:
- _year_ (1981–2023)
- _R_: recruitment (thousands), shifted back to birthyear
- _SSB_: spawning stock biomass (tonnes)
- _logSI_: ln(R/SSB) (natural log of survival index (SI = R/SSB))
- _CANUM1_: catch numbers of age-1, birth year, shifted back to birthyear
- _logCANUM1_: log(CANUM1) (natural log of catch numbers of age-1)
- _WSClag1_: annual WSC avg. lag-1 (WSC Jan-Dec avg. lagged 1 year ahead of spawning (e.g. WSC in 1980 shifted to year 1981))
- _WSClag1_MAM_: spawning season WSC avg. lag-1 (WSC Mar-May avg. lagged 1 year ahead of spawning)

3. **SAM_BW_2024.RData**  
SAM stock assessment output for blue whiting.   
Used for forecasting future SSB under 4 recruitment prediction scenarios (GM, REA1, WSCmax, WSCmin).

4. **myforecast.R**  
Custom function to forecast future SSB values using the SAM model output under fixed recruitment scenarios.  
Used in retrospective forecasting from 2005–2023, where recruitment predictions (from REA1, GM or WSC scenarios) are inserted into the SAM framework to evaluate how SSB would evolve under those recruitment conditions.  
Inputs: a) SAM model output, b) Recruitment value to impose.  

**Preparing_WSC_dataset_script2025.R**  
This script processes uwnd and vwnd wind components to compute wind stress curl (WSC).  
Outputs:

A) Spatial explicit WSC datasets (1981–2023) over ROI (Lon: -48 to 8, Lat: 40 to 75):  
- `df_monthly_roi`: monthly mean WSC  
- `df_annual_roi`: annual mean WSC  

B) WSClag1 indices (1981–2023) over BOX (Lon: -17 to -12, Lat: 52 to 55):  
- `WSClag1_JanDec_DF`: annual (Jan–Dec) WSC avg. lag-1 (e.g., WSC Jan–Dec in 1980 → WSClag1 in 1981)  
- `WSClag1_MAM_DF`: spawning season (Mar–May) WSC avg. lag-1
 

2) **BW_WSC_analyses_script2025.R**  
This script performs the main analyses: regression modeling, predictive skill assessment, and SSB forecasts.  
Key components:  
- Fits recruitment models: RI, E, A1, RE, REA1, GM
- Compares model performance: AICc, LogLik, DeltaAICc, Evidence Ration (ER)
- Predictive skill assessment: time series cross-validation of rescruitment predictions
- One-step ahead prediction errors (MSE, SD, mean error; diagnostic checks)
- Forecasts SSB using SAM under 4 rescruiment scenarios (GM, REA1, WSCmax, WSCmin)
- Fit CANUM1 models (A1, E, EA1) and repeat (model performance, predictive skill assessment, diagnostic on prediction errors)
- Evaluates WSC averaging timing: annual avg. (Jan–Dec) vs spawning season avg. (Mar–May) (lag-0 and lag-1)
