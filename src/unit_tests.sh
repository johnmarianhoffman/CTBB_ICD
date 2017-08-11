#!/bin/bash

# Baseline Reconstructions (tests downsampled matrices, no WFBP initialization, 10 iterations)

/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_out_1p0.yaml # Out, 1.0mm
/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_in_1p0.yaml # In, 1.0mm

/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_out_1p5.yaml # Out, 1.5mm
/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_in_1p5.yaml # In, 1.5mm

/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_out_3p0.yaml # Out, 1.5mm
/home/john/Code/icd/icd /home/john/Code/icd/tests/config_acr_in_3p0.yaml # In, 1.5mm
