#pragma once

struct recon_params configure_recon_params(std::string filename);
void load_raw(struct recon_params *  rp,struct ct_data * data);
void clean_up(struct recon_params * rp, struct ct_data * data);
