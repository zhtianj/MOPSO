##MOPSO

##present 
Optimization of backwashing management using a multi-objective particle swarm algorithm to simulate groundwater recharge to prevent clogging.

##code 
Feclogg.m
PSOFETSS.m
dfe3.mph
fe1.mph
fe2.mph

##Software installation and configuration
COMSOL Multiphysics（Recommended version 6.0 or higher）;
MATLAB（Recommended version R2022b or higher）;
COMSOL Multiphysics with MATLAB（LiveLink for MATLAB）;
IPhreeqcCOM.

##User Guide
1. Make sure that COMSOL Multiphysics and MATLAB are properly installed and that LiveLink for MATLAB is enabled.
2. After downloading IPhreeqcCOM from the official PHREEQC website, change the location of IPhreeqcCOM in the Feclogg file(ex:iphreeqc.LoadDatabase('E:\IPHREEQCCOM\database\phreeqc.dat').
3. The .mph file is downloaded and saved in the default path of COMSOL Multiphysics with MATLAB(ex:C:\Users\Administrator\.comsol\v62\llmatlab).
4. In MATLAB, open the main script file PSOFETSS.m.
5. Subsequent comments will modify the blocking or PSO parameters as appropriate
The script will interact with COMSOL Multiphysics to create, simulate and analyze the model. The results will be displayed in MATLAB commands and graphs.

##Testing
After running PSOFETSS, the Pareto frontier plot and other data are obtained.
