# Sleep Releases Action-Mode Network Revealed by Human Connectome



This repository contains the analysis code, atlas resources, example workflow, and released results used in this study. It includes the core procedures required to preprocess the data, construct individualized and group-level whole-brain functional atlases across sleep stages, quantify sleep-related topological reorganization, perform statistical modeling and classification, and reproduce the main analyses reported in the manuscript.

<p align="center">
  <img src="result/Figure/fig1.png" alt="Figure 1" width="1000">
</p>

---

## Overview

Sleep is an essential biological process that restores the brain. Although numerous advances have been made in uncovering the circuits regulating the sleep-wake cycle, the sleep-related dynamic reconfiguration of brain functional connectome topography remains unknown. Here, we collected simultaneous electroencephalogram and functional magnetic resonance imaging data, along with multiple behavioral assessments, from 130 healthy adults during nocturnal sleep. This rich dataset enabled us to construct a large-scale, sleep-state-dependent suite of system-level brain atlases and to establish their associations with sleep quality and pressure. We found that the action-mode network acts as a highly dynamic system that contracts as sleep deepens, most prominent in subcortical structures. Similar contractions were observed in humans and macaques under anesthesia, reflecting cross-species homologues of network reorganization during low-arousal states. Furthermore, the topological reconfiguration of the action-mode network is a stable biological marker that distinguishes sleep stages, reflects sleep quality, and predicts sleep pressure. Together, our work identifies the action-mode network as a dynamic core that reconfigures across the sleep-wake cycle, potentially serving as a neuro-target for intervening sleep disorders and tracking brain restoration.

This repository contains:

- code for data processing and analysis
- resources required for atlas construction
- example scripts demonstrating the full workflow
- released functional brain atlases for different sleep stages
- figure outputs used in the manuscript

---

## Table of Contents

- [Sleep Releases Action-Mode Network Revealed by Human Connectome](#sleep-releases-action-mode-network-revealed-by-human-connectome)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [Repository Structure](#repository-structure)
  - [Directory descriptions](#directory-descriptions)
  - [System Requirements](#system-requirements)
    - [Recommended operating system](#recommended-operating-system)
    - [Recommended hardware](#recommended-hardware)
  - [Software Dependencies](#software-dependencies)
    - [Core software](#core-software)
    - [Main neuroimaging and analysis dependencies](#main-neuroimaging-and-analysis-dependencies)
    - [Additional tools used in figure generation and post-processing](#additional-tools-used-in-figure-generation-and-post-processing)
  - [Official links](#official-links)
  - [Installation and Setup](#installation-and-setup)
    - [2. Install required external software](#2-install-required-external-software)
    - [3. Confirm required commands are available](#3-confirm-required-commands-are-available)
    - [4. Configure local paths](#4-configure-local-paths)
    - [5. Example Workflow](#5-example-workflow)
  - [Resources and Released Results](#resources-and-released-results)
  - [Notes](#notes)
  - [Citation](#citation)
  - [Contact](#contact)

---

## Repository Structure

```text
Sleep/
├── code/                              # Core functions used by scripts
│   └── postprocess/                   # fMRI post-processing and censoring pipeline
├── example/                           # Example workflow scripts
├── resource/                          # Resources required by the pipeline
│   └── Template/                      # Template atlases used for atlas construction
├── result/                            # Released results
│   └── Figure/                        # Figures used in the manuscript
├── scripts/                           # Main step-by-step analysis scripts
├── README.md                          # Repository documentation
├── sleep_info.mat                     # Paths and session identifiers for all raw sleep data
├── sleep_subject_list.txt             # Sleep participant list
├── subject_contral_select_7.mat       # Control-group participant list
└── subject_info.mat                   # Physiology, staging, and file-path metadata used in final analyses
```
---
## Directory descriptions

**code/**  
Stores reusable functions for processing and analysis.

**code/postprocess/**  
Implements the post-processing procedure applied after CIFTI conversion. The additional post-processing steps used in this study were:  
1. nuisance regression using global signal, averaged ventricular signal, averaged white matter signal, six motion parameters, and their temporal derivatives (18 regressors in total);  
2. removal of epochs with mean FD > 0.25 mm and frames with FD > 0.25 mm or DVARS over 50 > 75%;  
3. CIFTI smoothing using 6 mm FWHM surface smoothing for cortex and 8 mm FWHM parcel-constrained smoothing for subcortical volumes.

**example/**  
Provides an end-to-end example of the project workflow. The main entry file is:

```bash
example/example.sh
```

resource/Template/  
Contains the template atlases used to guide atlas construction.

result/Figure/  
Contains manuscript figures and visualization outputs.

scripts/  
Contains the main stepwise analysis scripts. These scripts perform the project logic and call functions from `code`.

---
## System Requirements

This repository is designed primarily for **Linux-based environments** and HCP/CIFTI-style neuroimaging workflows.

### Recommended operating system

- 64-bit Linux

### Recommended hardware

Runtime depends on dataset size and the number of epochs/sessions. For practical use, we recommend:

- multi-core CPU
- at least **32 GB RAM**
- sufficient disk space for raw, intermediate, and output neuroimaging files
- Docker support for preprocessing

---
## Software Dependencies

The code for analyzing the data and the functional brain atlases for each sleep stage are available in this repository. The study also used the following external tools and software.

### Core software

- **MATLAB R2022a**
- **R 4.4.2**
- **Ime4 1.1-37**
- **FSL 6.0.5**
- **Connectome Workbench 1.5.0**
- **Docker**

### Main neuroimaging and analysis dependencies

- **ciftify 2.3.3**  
  Used for structural and functional MRI preprocessing in CIFTI format.
- **Homologous-Functional-Regions**  
  Used for individualized cortical parcellation.
- **Neurosynth**  
  Used for meta-analysis of encroached regions.
- **SUIT toolbox**  
  Used for cerebellar surface visualization.
- **plot_fig_subcortex**  
  Used for 3-D scatter visualization of subcortical volumes.

### Additional tools used in figure generation and post-processing

- **Gephi**  
  Used for spring-embedded network visualization.
- **GraphPad Prism**  
  Used for statistical plotting.
- **Microsoft PowerPoint**  
  Used for manual final panel layout and figure composition.


---

## Official links

- [ciftify](https://ciftify.com/)  
  https://ciftify.com/edickie/ciftify

- [Homologous-Functional-Regions](https://github.com/MeilingAva/Homologous-Functional-Regions)

- [Neurosynth](https://github.com/neurosynth/neurosynth)

- [SUIT toolbox](https://www.diedrichsenlab.org/imaging/suit.htm)

- [plot_fig_subcortex](https://github.com/wd-veloce/plot_fig_subcortex)

- [Connectome Workbench](https://www.humanconnectome.org/software/connectome-workbench)

- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)

- [MATLAB](https://www2.mathworks.cn/products/matlab.html)

- [R](https://www.r-project.org)

- [Ime4](https://github.com/Ime4/Ime4)


---

## Installation and Setup

1. Clone the repository
```Bash
git clone https://github.com/BIT-YangLab/Sleep.git
cd Sleep
```


### 2. Install required external software

Please install the software listed in Software Dependencies, especially:

- Docker
- ciftify
- Connectome Workbench
- FSL
- MATLAB
- R and `lme4`

### 3. Confirm required commands are available

Depending on the modules you run, the following commands should be available in your environment:
```text
docker
wb_command
matlab
Rscript
fsl
```
### 4. Configure local paths

Before running the workflow on your own data, update the path definitions in scripts to match your local environment, including:
- raw data directory
- output directory
- temporary directory
- FreeSurfer license path
- software installation paths

**Important**  
Several scripts assume project-specific local paths and directory structures. Please review and modify path variables before execution.



### 5. Example Workflow

A complete example is provided in:
```Bash
example/example.sh
```

---

## Resources and Released Results

**resource/Template/**  
Contains template atlases required for atlas construction.

**result/**  
Contains released atlas results and other outputs generated in this study.

**result/Figure/**  
Contains figure files for the manuscript.

The repository includes the functional brain atlases for each sleep stage, including group-average and individualized atlas outputs used to support the study findings.

---

## Notes

- The example workflow shows the overall logic of the project, but local path configuration is required before execution.
- Some scripts are designed for project-specific data organization and may require adaptation for external datasets.
- For full logic and script-specific details, please inspect comments inside each script and the `INSTRUCTION` files inside directories.


---
## Citation

If you use this repository, please cite the corresponding manuscript:
```text
Yang, G.#, Zou, Q.H.#, Li, J.# et al. Sleep releases action-mode network revealed by human connectome.
```

Please also cite the external software tools used in your workflow where appropriate.

---

## Contact

For other scripts and data descriptions of Sleep repository, please see [INSTRUCTION](#) inside directory.

If you have any issues, please email **Guoyuan Yang** ([yanggy@bit.edu.cn](mailto:yanggy@bit.edu.cn)), we are willing to help you to solve your problem.

**Happy researching!**