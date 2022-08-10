# [FeatureCloud sPLINK](https://featurecloud.ai/app/splink)
## sPLINK workflow in FeatureCloud platform: Privacy-aware Tool for Genome-wide Association Studies__

**sPLINK** is a tool for federated, privacy-aware genome-wide association studies (GWAS). **sPLINK** is based on the [HyFed framework](https://github.com/tum-aimed/hyfed) and supports three association tests: chi-square, linear regression, and logistic regression. The sPLINK paper was published in Genome Biology: https://doi.org/10.1186/s13059-021-02562-1.


![state diagram](/data/images/state_diagram.png)


## Config file

```python
splink:
  local_dataset:
    data: hapmap1_100_1.bed
    phenotype: hapmap1_100_1_quant.pheno
    covariate: hapmap1_100_1.cov
  phenotype_name: pheno
  covariate_names: Age,Smoking
  logic:
    mode: file
    dir: .
  axis: 0
  use_smpc: false
  algorithm: linear-regression # logistic-regression # linear-regression # chi-square
  covariates: null
  chunk_size: 10 #  (x1000)
  max_iterations: 20
  cpu_cores: 1
  result:
    manhattan: manhattan-plot
    output: linear-regression-result.csv # logistic-regression-result.csv or linear-regression-result.csv or chi-square-result.csv

```


## Experiments
Different experiments can be done using FC sPLINK app based on the association tests. For each test, one should use
other parameters in the config file. All the experiments are covered in the sample data directory.

### Chi-Square

```python
splink:
  local_dataset:
    data: hapmap1_100_1.bed
    phenotype: hapmap1_100_1.pheno # should be binary
    covariate: hapmap1_100_1.cov # Irrelevant
  phenotype_name: pheno
  covariate_names: Age,Smoking # Irrelevant
  logic:
    mode: file
    dir: .
  axis: 0 # Irrelevant
  use_smpc: false
  algorithm: chi-square
  covariates: null # Irrelevant
  chunk_size: 10 
  max_iterations: 20
  cpu_cores: 4
  result:
    manhattan: manhattan-plot
    output: chi-square-result.csv # logistic-regression-result.csv or linear-regression-result.csv or chi-square-result.csv

```

### Linear regression
```python
splink:
  local_dataset:
    data: hapmap1_100_1.bed
    phenotype: hapmap1_100_1_quant.pheno
    covariate: hapmap1_100_1.cov
  phenotype_name: pheno
  covariate_names: Age,Smoking
  logic:
    mode: file
    dir: .
  axis: 0
  use_smpc: false
  algorithm: linear-regression
  covariates: null
  chunk_size: 10 #  (x1000)
  max_iterations: 20
  cpu_cores: 4
  result:
    manhattan: manhattan-plot
    output: linear-regression-result.csv # logistic-regression-result.csv or linear-regression-result.csv or chi-square-result.csv

```

### Logistic Regression

```python
splink:
  local_dataset:
    data: hapmap1_100_1.bed
    phenotype: hapmap1_100_1.pheno # should be binary
    covariate: hapmap1_100_1.cov
  phenotype_name: pheno
  covariate_names: Age,Smoking
  logic:
    mode: file
    dir: .
  axis: 0
  use_smpc: false
  algorithm: logistic-regression
  covariates: null
  chunk_size: 10 #  (x1000)
  max_iterations: 20 # important
  cpu_cores: 1
  result:
    manhattan: manhattan-plot
    output: logistic-regression-result.csv # logistic-regression-result.csv or linear-regression-result.csv or chi-square-result.csv

```


### Run sPLINK

#### Prerequisite

To run the sPLINK application you should install Docker and featurecloud pip package:

```shell
pip install featurecloud
```

Then either download the sPLINK image from featurecloud docker repository:

```shell
featurecloud app download featurecloud.ai/splink
```

Or build the app locally:

```shell
featurecloud app build featurecloud.ai/splink
```

You can use provided example data or your own data. And provide the desired settings in the `config.yml` file.

#### Running app

You can run sPLINK as a standalone app in the [FeatureCloud test-bed](https://featurecloud.ai/development/test) or [FeatureCloud Workflow](https://featurecloud.ai/projects). You can also run the app using CLI:

```shell
featurecloud test start --app-image featurecloud.ai/splink --client-dirs './sample_data/lg/c1,./sample_data/lg/c2' --generic-dir './sample_data/lg/generic'
```
For running sPLINK with Linear regression and Chi-square tests, use `lr` and `chi-square` subdirectories, correspondingly, 

```angular2html
@article{nasirigerdeh2022splink,
title={sPLINK: a hybrid federated tool as a robust alternative to meta-analysis in genome-wide association studies},
author={Nasirigerdeh, Reza and Torkzadehmahani, Reihaneh and Matschinske, Julian and Frisch, Tobias and List, Markus and Sp{\"a}th, Julian and Weiss, Stefan and V{\"o}lker, Uwe and Pitk{\"a}nen, Esa and Heider, Dominik and others},
journal={Genome Biology},
volume={23},
number={1},
pages={1--24},
year={2022},
publisher={BioMed Central}
}


@misc{nasirigerdeh2021hyfed,
    title={HyFed: A Hybrid Federated Framework for Privacy-preserving Machine Learning},
    author={Reza Nasirigerdeh and Reihaneh Torkzadehmahani and Julian Matschinske and Jan Baumbach and Daniel Rueckert and Georgios Kaissis},
    year={2021},
    eprint={2105.10545},
    archivePrefix={arXiv},
    primaryClass={cs.LG}
}
```
