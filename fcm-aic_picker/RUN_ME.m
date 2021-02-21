%     Author: Eduardo Valero Cano
%     ---------------------------
%     Supplementary material for the manuscript "Automatic seismic phase
%     picking based on unsupervised machine learning classification and
%     content information analysis" submitted for peer-review in GEOPHYSICS.
%     November 2020.

%     RUN_ME.m: creates and process one project. 

project_name = 'synthetic_1';

%% DO NOT EDIT BELOW THIS LINE --------------------------------------------
manage_project(project_name, 'create');
process_project(project_name, 'all');
process_project(project_name, 'single phase');