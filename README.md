# Data and code for *Designing epigenetic clocks for wildlife research*

**Abstract**: The potential applications of epigenetic clocks are expanding in wildlife conservation and management. The pace at which they are being adopted highlights the need for field-specific design best practices. Epigenetic clocks were originally developed for human studies, presenting challenges for their adoption in wildlife research. Most notably, the estimated ages of sampled wildlife can be unreliable, and sampling restrictions limit the number and variety of available samples, which can reduce the accuracy of epigenetic clocks for wildlife. In this article, we present a detailed workflow for designing, validating, and applying wildlife epigenetic clocks in a way that accounts for sampling constraints. We provide recommendations for two main applications of wildlife epigenetic clocks: estimating unknown ages and assessing cumulative biological aging. Our simulations and analyses applied to an extensive polar bear dataset from across the Canadian Arctic demonstrate that accurate epigenetic clocks can be constructed and validated using wildlife samples despite their limitations. With our workflow and examples, we hope to make the use of epigenetic clocks more accessible and widespread in wildlife conservation and management.

The code is separated into nine R scripts that recreate the analyses and figures presented in the paper. There are eight functions in the /functions folder required to run these analyses. Scripts 1-2 create simulated data and run the simulations in boxes 2-4. Scripts 3-4 process real DNA methylation data from 10 polar bear subpopulations in the Canadian Arctic into the format required for creating epigenetic clocks (see workflow steps A and B in the manuscript). Scripts 5-9 demonstrate the last two parts of the workflow (C and D) covered in the paper using the polar bear data. An overview of the data are in Box 1. Due to the size of the polar bear dataset, it is not available in this repository, but will be made available on Dryad.

All code was run in R v4.3.1 and runs on macOS Monterey v12.7.3

**CRAN Package versions**:

* tidyverse v2.0.0
* cowplot v1.1.3
* glmnet v4.1-8
* limma v1.1.3
* BiocManager v1.30.22

**Developmental Package versions**:

* packages/HorvathMammalMethylChip40anno.test.unknown_0.2.2.tar.gz
* packages/HorvathMammalMethylChip40manifest_0.2.2.tar.gz

Note that both package versions were downloaded and are available from https://github.com/shorvath/MammalianMethylationConsortium

**/functions**:

* fitClock.R: Fit clock using training and testing separate training and testing data, either returning only accuracy metrics or plot as well
* fitLimma.R: EWAS for feature selection
* fitLOOClock.R: Compare leave-one-out cross validation and hold-out clocks
* groupLimma.R: Fit EWAS for primary variable controlling for other variables
* makePpnPlot.R: Helps visualize change in accuracy metrics with different training and testing sets
* sampleGrps2.R: Helps sample observations by group
* simpCpGs.R: simulate methylation data
* testAccuracy.R: Test accuracy depending on bias
* normalizeBetas.R: Normalizes raw methylation data from idat files

**/R**: 

* 01 -- Create Fig.1 plots to demonstrate clock accuracy and the difference between a clock with high R2 and high MAE vs. high R2 and low MAE
* 02 -- Run simulations and make plots for boxes: class/age biases and feature selection
* 03 -- Normalize raw DNA methylation data in idat files
* 04 -- Combine sample information with normalized DNA methylation betas
* 05 -- Effects of sampling bias on clock performance. Shows how age bias, sex bias, tissue bias, genetic variation, and aging error can affect clock performance using the polar bear data
* 06 -- EWAS to find CpGs correlated with sample characteristics. Remove any probes correlated (p < 0.05) with tissue, sex, and population andkeep top ~ 10% of probes correlated with age
* 07 -- Effects of feature selection on clock accuracy. Cut CpG sites based on alignment to the polar bear genome and EWAS and then test implications for polar bear clock accuracy
* 08 -- Plot approaches for feature selection. Shows how EWAS and aligning the genome can impact clock accuracy
* 09 -- Clock validation procedures. Simulate the process of sampling a population and validating a clock using leave-one-out cross-validation versus a true hold-out set

**/input**:

* full_sibs.rds: Character vector of individual IDs with siblings in dataset. These individuals are normally excluded from clock design because of potential bias due to genetic relatedness

* low_qc_positions.rds: Data frame with low-quality detection p-values (high values indicate quality issues) for screening samples
  * position: chip position
  * chip.ID.loc: chip IDs corresponding to name of position
  * detection_p: detection p-value; > 0.05 is high
  * batch number: Ranges from 1-8

* norm_betas.rds: Data frame with polar bear DNA methylation data (**NOTE:** this file is > 100 MB and therefore not available on GitHub. It will eventually be made available on Dryad, but can also be recreated locally by following the normalization and data cleaning steps in scripts 03 and 04.)
  * sampleId: Unique sample name
  * chip.ID.loc: chip IDs corresponding to name of position
  * id: ID for individual bear (might be multiple sampleIds/bear)
  * Spec: Tissue type (blood, skin, muscle)
  * YMD: Year-month-day of sample
  * sex: Sex of individual
  * age: Age of individual at time of sample
  * Born: Birth year of individual
  * Population: Code of population to which individual belongs (correspond to populations in Box 1)
  * cgXXXXXXXX: Normalized betas (methylation) at CpG sites

* pb_alignment.rds: Data frame with sites aligning to the polar bear genome

**/sample_sheets**:

* batchXX_samples.rds: Data frame with sample information corresponding to each sample in sample sheets.
  * sampleId: Unique sample name
  * id: ID for individual bear (might be multiple sampleIds/bear)
  * Spec: Tissue type (blood, skin, muscle)
  * YMD: Year-month-day of sample
  * sex: Sex of individual
  * age: Age of individual at time of sample
  * Born: Birth year of individual

* PB_arrayXX_sample_sheetXX.rds: Data frame with information about the array location corresponding to each sample in *batchXX_samples.rds*.
  * Sample_Name: Unique sample name (corresponds to sampleId in *batchXX_samples.rds*)
  * Sample_Well: Location on 96-well plate
  * Sample_Plate: Identifies the plate for batches with multiple plates
  * chip.No: Unique identifier for the chip (four per array)
  * chip.ID: Numeric identifier for the chip (four per array)
  * stripe: Row and column position of the sample on the chip (RXX = rows 1-6, CXX = columns 1-2)
  * row: Row identifier for 96-well plate (8 rows, a-h)
  * column: Column identifier for 96-well plate (12 rows, numbered)

* updated_sample_sheet_PB_arrayXX.rds: Data frame identical to corresponding original sample sheets (*PB_arrayXX_sample_sheetXX.rds*), with the addition of 'Basename' column pointing to the folder location of the idat file.

**/iscans**: Contains raw DNA methylation image files organized in batches

**/packages**: Contains developmental versions of packages needed for normalization

**/temp**: Resampled biased polar bear clocks from 03-TestClockBiases.R, for reproducing figures

* age_error_model.rds: Age estimation error
* ageGrp_bias_model.rds: Age bias (old predict young)
* Population_bias_model.rds: Population/genetic bias
* sex_bias_model.rds: Sex bias
* Tissue_bias_model.rds: Tissue bias (blood/skin/muscle)