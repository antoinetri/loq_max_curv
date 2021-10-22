# loq_max_curv
A small python code to calculate the lower limit of quantification (LOQ) of a pXRF calibration. To the best of our knowledge, LOQ are generally calculated using a statistical approach (based on statistical LOD times 10/3). Our appraoch is to determine it using an emprical approach. However, this implies to choose arbitrarily the LOQ. Here, we propose an original approach to do it semi-automatically by finding the maximum of curvature point on a power line trend fitted on a PoD vs ref conc plot. Feel free to contact me if you want more info about the applciability to your case study.

This is how data should look like before importing it:
![image](https://user-images.githubusercontent.com/14851413/138432440-fea9e5ea-078c-4489-b4a1-a2102c6b95d0.png)

Example for Al2O3
![image](https://user-images.githubusercontent.com/14851413/138432651-e42c2888-effd-45a9-b935-54a6db4d0a53.png)

Example for Zn: 
![image](https://user-images.githubusercontent.com/14851413/138432774-bcfc12cd-318d-46bf-8c46-a1e5c02e055d.png)

Example for the Y
![image](https://user-images.githubusercontent.com/14851413/138432898-fd1cf351-59e5-4c97-9e0e-d14dc666c870.png)

The red dots are the data, the blue line is the power line fitted to the data, the dashed black line points to the concentration at the maximum curvature point (mcp), the green dashed line points to the concentration where the percentage of difference (PoD) value equals 25%. For more information, please check out our review paper on pXRF calibrations: workflow and goog pratices to optimize the analysis of geological samples. This is prep. at the moment, link will follow asap.

