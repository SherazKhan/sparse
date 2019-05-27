## Source Space and Forward_Solutions
1. Define cortical surface source spaces
a. Use standard steps (incl. patch decompositions) outlined in MNE manual 
b. Same as used by SPIGH code 
2. Define subcortical volume source spaces: 
a. In Shell
i. Reconstruct MRIs via Freesurfer Stream
ii. Setup a standard subject folder based on MNE/Freesurfer manual instructions
iii. This will have all info to extract volumes of major subcortical structures
b. In MATLAB: Create .mgz, .txt and .mat files with masks, spatial coordinates and geometric information for each subcortical region’s source spaces. Run exec_srcfwd_final. This calls
i. read_save_sfwd to perform dimensionality reduction, split into patches/subvolumes
ii. check_subdivs to create mgz files for each subdivision in mri folder
iii. srcspace_figure to generate mgz file illustrating subdivided source space figure
3. This codestream will lead you to run shell script functions as needed too:
a. Generate cortical forward solutions using standard steps from MNE manual
b. Run scort_srcfwd_exec for subcortical volumes, hippocampus and brainstem
4. Scripts to generate figures in the paper are under \FigureS1
5. In Freeview: view .mgz, label and surface files together to check if all patch/subvolume/surface source spaces are done properly, so as to be sure what regions forward solutions correspond to
6. In Freeview, visualize subcortical source space masks to ensure all done correctly

# Data_Process: Remove artifacts, Filter Data; Evoked Averages; Compute Covariances
## Inverse_Solutions
MATLAB code to read data, process source/noise covariances, define parameter choices, obtain MNE estimates, obtain SPIGH estimates for given hierarchy, transition across hierarchies
1. Scripts to run cortical inverse solution (within and across all cortical hierarchies)
2. Scripts to visualize cortical inverse solution (time series, spatiotemporal movies, performance metrics)
3. Scripts to run subcortical inverse solution (with reduced cortical source space and all deep sources)
4. Scripts to visualize end joint inverse solution (time series, spatiotemporal movies, performance metrics)
Note: running this code will display instructions on shell commands to run as and when required. 
Core_Functions: General Sub-Functions Called for Angles, Field Maps, Algorithm, Visualization, ETC
