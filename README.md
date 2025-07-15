# contam-x-jr-ml

This is a MATLAB implementation of *contam-x-jr*.  

This version (1.0.0) is a simplified version of *AIRNET*.  
It only implements orifice and leakage area elements which limits calculations to naturally driven airflows.  
This version does not perform contaminant calculations.  

The following verification test cases are based on those provided in Appendix B of NISTIR-89-4072.  
- *airnet_pl1.m*: Powerlaw Element Test #1  
- *airnet_pl2.m*: Powerlaw Element Test #2 
- *airnet_pl3.m*: Powerlaw Element Test #3  
- *airnet_stack1.m*: Stack Effect Test #1  

These tests can be run using the MATLAB Test Browser.  
