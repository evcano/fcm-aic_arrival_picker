DESCRIPTION
-----------

This repository contains the code and data used in E.V. Cano, J. Akram, D.B. Peter; Automatic seismic phase picking based on unsupervised machine learning classification and content information analysis; submitted for peer-review in GEOPHYSICS (November 2020).

This repository includes four 3C downhole microseismic datasets (three synthetic and a few real events) and three different picking algorithms. Additionally, codes for hypocenter location using damped least-squares are included.


HOW TO USE
----------
- To use a picking algorithm, execute the script 'RUN_ME.m' inside the desired algorithm's directory. Set the 'project_name' variable as 'synthetic_1', 'synthetic_2', 'synthetic_3', or 'real' depending on the dataset to be processed.

- To locate the hypocenters of the synthetic datasets events, use the script 'locate_hypocenters.m' inside the 'hypocenter_location' directory. The hypocenter location is carried using the picks computed by the FCM-AIC picking algorithm.
