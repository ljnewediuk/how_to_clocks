# Data and code for *Designing epigenetic clocks for wildlife research*

**Abstract**: The potential applications of epigenetic clocks are expanding in wildlife conservation and management. The pace at which they are being adopted highlights the need for field-specific design best practices. Epigenetic clocks were originally developed for human studies, presenting challenges for their adoption in wildlife research. Most notably, the estimated ages of sampled wildlife can be unreliable, and sampling restrictions limit the number and variety of available samples, which can reduce the accuracy of epigenetic clocks for wildlife. In this article, we present a detailed workflow for designing, validating, and applying wildlife epigenetic clocks in a way that accounts for sampling constraints. We provide recommendations for two main applications of wildlife epigenetic clocks: estimating unknown ages and assessing cumulative biological aging. Our simulations and analyses applied to an extensive polar bear dataset from across the Canadian Arctic demonstrate that accurate epigenetic clocks can be constructed and validated using wildlife samples despite their limitations. With our workflow and examples, we hope to make the use of epigenetic clocks more accessible and widespread in wildlife conservation and management.

All code was run in R v4.3.1 and runs on macOS Monterey v12.7.3

**Package versions**:

* tidyverse v2.0.0
* cowplot v1.1.3
* glmnet v4.1-8
* limma v1.1.3

**/functions**:

* fitClock.R: Fit clock using training and testing separate training and testing data, either returning only accuracy metrics or plot as well
* fitLimma.R: EWAS for feature selection
* fitLOOClock.R: Compare leave-one-out cross validation and hold-out clocks
* groupLimma.R: Fit EWAS for primary variable controlling for other variables
* makePpnPlot.R: Helps visualize change in accuracy metrics with different training and testing sets
* sampleGrps2.R: Helps sample observations by group
* simpCpGs.R: simulate methylation data
* testAccuracy.R: Test accuracy depending on bias

**/R**: 

* 01 -- Create Fig.1 plots to demonstrate clock accuracy and the difference between a clock with high R2 and high MAE vs. high R2 and low MAE
* 02 -- Run simulations and make plots for boxes: class/age biases and feature selection
* 03 -- Effects of sampling bias on clock performance. Shows how age bias, sex bias, tissue bias, genetic variation, and aging error can affect clock performance using the polar bear data
* 04 -- EWAS to find CpGs correlated with sample characteristics. Remove any probes correlated (p < 0.05) with tissue, sex, and population andkeep top ~ 10% of probes correlated with age
* 05 -- Effects of feature selection on clock accuracy. Cut CpG sites based on alignment to the polar bear genome and EWAS and then test implications for polar bear clock accuracy
* 06 -- Plot approaches for feature selection. Shows how EWAS and aligning the genome can impact clock accuracy
* 07 -- Clock validation procedures. Simulate the process of sampling a population and validating a clock using leave-one-out cross-validation versus a true hold-out set

**/input**:

* full_sibs.rds: Character vector of individual IDs with siblings in dataset. These individuals are normally excluded from clock design because of potential bias due to genetic relatedness

* low_qc_positions.rds: Data frame with low-quality detection p-values (high values indicate quality issues) for screening samples
  * position: chip position
  * chip.ID.loc: chip IDs corresponding to name of position
  * detection_p: detection p-value; > 0.05 is high
  * batch number: Ranges from 1-8

* norm_betas: Data frame with polar bear DNA methylation data (**NOTE:** this file is > 100 MB and therefore not available on GitHub. It will be made available on Dryad.)
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

* pb_alignment: Data frame with sites aligning to the polar bear genome

**/temp**: Resampled biased polar bear clocks from 03-TestClockBiases.R, for reproducing figures

* age_error_model.rds: Age estimation error
* ageGrp_bias_model.rds: Age bias (old predict young)
* Population_bias_model.rds: Population/genetic bias
* sex_bias_model.rds: Sex bias
* Tissue_bias_model.rds: Tissue bias (blood/skin/muscle)