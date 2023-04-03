# WFSim Code for wind farm simulation#

- Scripts for MATLAB simulation
	- WFSim_demo: One simulation run
    - RunWFSim_Koopman_Rsweep: Multiple runs for sweep over different Koopman 
      sys ID models and different R weighting values

- bin (TU Delft + TUHH): 
  Folders with controller code	
- controlDefinitions (TU Delft + TUHH): 
  m and mat files for control output/plant input signals
  thrust coefficient (Ct) and yaw (phi)
- developmentTools (TU Delft): 
  find optimal set point for yaw for different 
  farm layout configurations and wind speeds (currently data missing) 
- layoutDefinitions (TU Delft):
  wind farm layouts
- libraries (TU Delft + TUHH): 
  auxiliary functions for plotting, saving data as mat files and analysis
- OldCode (TUHH): 
  Some test code (e.g. compared speed of qLPV-MPC and MPC solved 
  with optimize from YALMIP)
- Yalmip-master (added by TUHH)
  
- Write a short file about the content of this folder