# how_to_clocks

## Outline:

* Normal pre-processing workflow will come from previous projects, so we will start here with the normalized betas from all populations

### Box 1

Shows how age bias, sex bias, tissue bias, genetic variation, and aging error can affect clock performance

* Script 1: WH clock with no bias for comparison (to clocks from scripts 1-5)
* Script 2: Build clocks using only young or only old individuals and use them to predict the other age class
* Script 3: Build clocks using only females or only males and use them to predict the other sex
* Script 5: Introduce error in ages, then build a clock with the inaccurately-aged individuals
* Script 6: Build clock using only Beaufort bears and use it to predict ages of bears in the WH clock and vice-versa
* Script 4: Build clocks using only skin, blood, or muscle from all populations and use them to predict the other tissues

### Box 2

Shows how feature selection methods might improve accuracy of clock

* Script 1: Clock built with samples from all populations and all tissues with EWAS and genome alignment for comparison
* Script 2: Clock with only EWAS (to test how genome alignment affects performance)
* Script 3: Clock with only genome alignment (to test how EWAS affects performance)

### Box 3

Shows how different validation methods perform for validation

* Script 1: Clock built with all samples from all populations with small hold-out set for comparison
* Script 2: Validate clock using cross-validation, then use it to predict the hold-out, and compare accuracy. This will test how much accuracy is inflated when we only compare within the same dataset and don't test on new data.
* Script 3: Build clock using stability selection, then use it to predict the hold-out. This will test how much stability selection can improve our clocks.
