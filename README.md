# 3D_Rift_NCI

        
## Objective

      - Run the model script 3D_Rift_NCI.py at different resolutions and under different compute conditions.
     
## Requirements

      - Underworld2 >= v2.5.1b https://github.com/underworldcode/underworld2
          Refer to docs/install_guides/nci_raijin.sh for build instructions
      - UWGeodynamics >= 0.5.0 https://github.com/rbeucher/UWGeodynamics
          As of writing this one must use the repository to download the latest version and locally use pip to install UWGeodyanmics, command:
          git clone https://github.com/rbeucher/UWGeodynamics.git ; pip install -e UWGeodynamics --user
          
## How to run

      - A pbs_script compatible with raijin is given.
      - The user must modify the PYTHONPATH env. variable in the pbs_script to include the UWGeodynamics directory
      Please modify the following parts of the pbs_script when benchmarking the model.
        - Different model resolutions by modifying the UW_RESFACTOR environment variable. UW_RESFACTOR is used as a factor for the model resolution. 1 -> (128,64,64), 2-> (256, 128, 128).
        - Different compute conditions by modifying the modules and pbs conditions.
