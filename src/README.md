# Code

This directory contains the code required to replicate the computational experiments.

The RO folder contains the code for the proposed robust optimization method.

The RS folder contains the code for the proposed robust satisficing method.

To run the experiments:

    1. In the ConGenerAlg file, set data_path to the directory containing the input data and result_path to the directory where the computational results should be saved. For the RS model, set obj_path to the path of the file containing the prescribed objective values.
    
    2. Run the main function in the ConGenerAlg file.