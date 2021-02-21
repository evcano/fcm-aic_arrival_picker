FCM-AIC WAVE ARRIVAL PICKER
---------------------------
Copyright (C) November 2020  Eduardo Valero Cano,
King Abdullah University of Science and Technology (KAUST).

See LICENSE file.
============================================================================


SYSTEM REQUIREMENTS
============================================================================
- MATLAB R2016b or new versions.
- MATLAB Signal Processing toolbox.
- MATLAb Computer Vision toolbox (for the ransac function).


TOOL DESCRIPTION
============================================================================
FMC-AIC WAVE ARRIVAL PICKER is a MATLAB tool for P- and S-wave arrival pick-
ing on three-component (3C) downhole microseismic data. It uses the fuzzy c-
means (FCM) clustering algorithm to determine P- and S-wave arrival windows
and the Akaike Information Criterion (AIC) to pick the arrival times.

For a description of the theory behind this tool read: E.V. Cano, J. Akram,
D.B. Peter. Automatic seismic phase picking based on unsupervised machine 
learning classification and content information analysis. Submitted for peer
-review in GEOPHYSICS.


HOW TO USE THIS TOOL
============================================================================
CREATING A PROJECT
------------------
The first step to pick the arrivals of a 3C downhole microseismic dataset is
to create a project by typing "manage_project('project_name','create')" on 
the MATLAB command window. To create a project, the following is required:

    - Waveforms in SAC format or MAT files (one MAT file per waveform with 
      one three-row array [first row: time axis; second row: amplitude valu-
      es; third row: header information]). The data needs to be organized as
      follows:

        data_dir/
            event_dir_1/
                    event_name.receiver_1.channel.FORMAT
                    ...
                    event_name.receiver_N.channel.FORMAT
            ...
            event_dir_N/
                    event_name.receiver_1.channel.FORMAT
                    ...
                    event_name.receiver_N.channel.FORMAT

    - Project parameters file. The file name must be "project_name.par".
      For a description of the parameters see the file "real.par".

A project file "project_name.project.mat" will be generated after the proje-
ct creation. The project file contains the following MATLAB variables:

    - "par". Struct containing the project parameters.

    - "events_info". Table with information about the microseismic events.
      The table fields for each event are:

           - "id". Number that identifies the event.

           - "name". Name of the event. It is equal to the name of the event
             directory "event_dir".

           - "dir". Path to "event_dir".

           - "status". Indicates if the picking of the event arrivals has
             been done ('D'), is pending ('P'), or the picks are unident-
             ified ('U'). 

Aditionally, a directory "project_name_results" will be created. All the pr-
oject results are stored here.


PROCESSING A PROJECT
--------------------
To run the arrival picker type "process_project('project_name,'option')" on
the MATLAB command window. The available options are:

    - "pending". Picks the events with status 'P'.

    - "single phase". Picks the events with status 'U'.

    - "all". Picks all the events disregarding the status and overwriting 
      previous results.

    - "id". Runs the arrival picker on the event with id "id" overwriting
      the previous result.

For each event, a results file "id_event_name_results.mat" will be stored on
the "project_name_results" directory. The results file contains a struct with
the following fields:

    - "(p/s/u)_window". Two element array with the start- and end-
      sample of the P-, S-, or unidentified-wave window.

    - "(p/s/u)_pick". P-, S-, or unidentified-wave arrival pick s-
      ample.

    - "(p/s)_moveout". Estimated P- or S-wave arrival moveout.

    - "(p/s)_moveout_ni". Number of inliers of the P- or S-wave
      arrival moveout.

    - "(p/s)_moveout_coef". Coefficients of the quadratic polyn-
      omial describing the P- or S-wave arrival moveout.

    - "pick_flag". (1) Unidentified pick relabeled as P pick.
      (2) Unidentified pick relabeled as S pick. (3) P pick rel-
      abeled as S pick.

    - "baz". Estimated event back-azimuth on each receiver.

If set on the parameters file, all figures will also be stored on the 
"project_name_results" directory.


CHANGING A PROJECT PARAMETERS
-----------------------------
To change the parameters of a project, edit the project parameters file and
write "manage_project('project_name','update') on the MATLAB command window.
This will update the project file.


SCRIPTS OVERVIEW
============================================================================
PSEUDO CALL TREE
----------------
- manage_project.m
- process_project.m
    - identify_single_phase_event.m
    - wave_arrival_picker.m
        - compute_trace_features.m
        - fuzzy_c_means.m
        - get_arr_times.m
        - get_arr_windows.m
        - get_candidate_arr_windows.m
    - verify_phase_labeling.m
- minmaxn.m
- plot_waveforms.m
- read_waveforms.m


DESCRIPTION
-----------
- compute_trace_features.m. computes the mean amplitude absolute value, peak
  power spectral density, and short- long-term average ratio of a one trace.

- fuzzy_c_means.m. Fuzzy c-means clustering algorithm.

- get_arr_times.m. Picks the P- and S-wave arrival times.

- get_arr_windows.m. Determines P- and S-wave arrival windows from the cand-
  idate arrival windows.

- get_candidate_arr_windows.m. Determines windows of possible (candidate) w-
  ave arrivals.

- identify_single_phase_event.m. Labels the picks of single-phase events as
  P or S picks.

- manage_project.m.  Creates a project or updates its parameters.

- minmax.m. Normalizes a vector usin min-max normalization.

- plot_waveforms.m. Plots one or more 3C records.

- process_project.m. Runs the wave arrival picker on a project.

- read_waveforms.m. Reads one or more 3C records on SAC format to a MATLAB
  array of dimensions n_traces x n_samples. The order of the components is
  "east,north,vertical" from the top to the bottom row.

- verify_phase_labeling.m. Labels unknown picks as P or S picks and corrects
  S picks labeled as P picks.

- wave_arrival_picker. Picks the P- and S-wave arrival times on one 3C reco-
  rd.
