
------------------------------------------------------------------------------------------------------
# SCARP method for scATAC-seq data analysis
Incorporating network diffusion and peak location information for better single-cell ATAC-seq data analysis.

The SCARP method project is deposited [here](https://github.com/Wu-Lab/SCARP).

## Authors
- yujiating@amss.ac.cn
- lywu@amss.ac.cn


## Reproduce results
For reproducibility, we provide all the necessary scripts and data.

SCARP has a total of 5 cases, as indicated in the directory.
```
.
├── Exp1_Benchmark                  
├── Exp2_Robustness                   
├── Exp3_SNARE_seq            
├── Exp4_SOX10_Knockdown         
└── Exp5_10X_Multiome                   
```

To reproduce the code for the case you desire, simply enter the corresponding folder and run the code step by step. 

For example, if you want to reproduce the results of SCARP on benchmark scATAC-seq data, enter the folder labeled `Exp1_Benchmark`, and run the code in the following sequence:
- S01_Data_Preprocessing.ipynb
- S02_Run_SCARP.ipynb
