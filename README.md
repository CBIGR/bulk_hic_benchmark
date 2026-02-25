# EVALUATION OF DEEP LEARNING TOOLS FOR CHROMATIN CONTACT PREDICTION

## Project Overview

This project evaluates the performance of five deep learning models for chromatin contact prediction:

- C.Origami  
- Chromafold  
- HiCdiffusion  
- GRACHIP  
- Epiphany  

Evaluation metrics include accuracy, image quality, and loop detection downstream analysis using true Hi-C maps.

## Methodology

For each model:
1. Preprocess data following the original paper’s method.
2. Train and perform inference using the model’s code.
3. Compare predicted contact maps with true Hi-C maps across:
   - Accuracy performance  
   - Image quality  
   - Loop detection (using 4 loop callers: HICCUPS, SIP, LASCA, Mustache)
  
![process drawio](https://github.ugent.be/vermeirssenlab/MasterThesis_Thu/assets/14965/08abcf07-9b45-437b-838c-e0fb8038e3e8)


## Repository Structure

The repository is organized as follows:
```
├── Epiphany/
├── C.Origami/
├── Chromafold/
├── HiCdiffusion/
├── GRACHIP/
├── Loop detection/
  ├── mustache/
  ├── SIP/
  ├── LASCA/
  └── HICCUPS/
└── comparison analysis/
  ├── loop analysis/
  └── Performance analysis/
```
- Model folders: Each contains code for data preprocessing, training, and inference specific to the model.

- Loop detection folder: Contains code to run four loop callers (HICCUPS, SIP, LASCA, Mustache) for downstream loop detection analysis.

- Comparison analysis folder: Contains Jupyter notebooks and scripts to perform and visualize performance comparisons and loop detection analysis results.

## Usage
Each model folder includes detailed instructions on how to preprocess data, train the model, and generate predictions. Loop detection folders provide scripts to run loop callers on predicted contact maps. Finally, the comparison analysis notebooks demonstrate how to reproduce the evaluation results.

## Reference
Source codes:
- C.Origami: [https://github.com/tanjimin/C.Origami.git](https://github.com/tanjimin/C.Origami.git)
- Epiphany: [https://github.com/arnavmdas/epiphany.git](https://github.com/arnavmdas/epiphany.git)
- HiCDiffusion: [https://github.com/SFGLab/HiCDiffusion.git](https://github.com/SFGLab/HiCDiffusion.git)
- ChromaFold: [https://github.com/viannegao/ChromaFold.git](https://github.com/viannegao/ChromaFold.git)
- GRACHIP: [https://github.com/Ruoyun-W/GRACHIP.git](https://github.com/Ruoyun-W/GRACHIP.git)
- SIP: [https://github.com/PouletAxel/SIP/wiki/SIP-Quick-Start](https://github.com/PouletAxel/SIP/wiki/SIP-Quick-Start) 
- HICCUPS: [https://github.com/aidenlab/juicer/wiki/HiCCUPS](https://github.com/aidenlab/juicer/wiki/HiCCUPS)
- LASCA: [https://github.com/ArtemLuzhin/LASCA_pipeline.git](https://github.com/ArtemLuzhin/LASCA_pipeline.git)  
- Mustache: [https://github.com/ay-lab/mustache.git](https://github.com/ay-lab/mustache.git)
