#pragma once

// Rotate slices from fixed cartesian grid into rotating coordinate frame
// e.g. Initial WFBP slices into ICD rotating frame
void rotate_slices_fixed2rotating(const struct recon_params * rp,struct ct_data * data);

// Return rotating slices to fixed cartesian grid
// e.g. Reconstructed ICD slices back for final save to disk
void rotate_slices_rotating2fixed(const struct recon_params * rp,struct ct_data * data);
