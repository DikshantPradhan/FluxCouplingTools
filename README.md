# FluxCouplingTools

This repository contains the CachedFCF and FALCON method as well as other tools used in the paper "Efficient enzyme coupling algorithms identify functional pathways in genome-scale metabolic models". This package also requires the following packages:

sybil,
parallel,
gurobi,
dplyr,
pryr,
data.tree.

A list of coupled reaction sets can be calculated through the following example code:

grb_model <- convert_sybil_to_gurobi(sybil_model)
coupling_mtx <- flux_coupling_raptor(grb_model)$coupled
sets <- get_list_of_sets_from_mtx(coupling_mtx)

A list of coupled gene and reaction sets can be calculated through the following example code:

falcon_model <- GRB_generate_falcon_model(sybil_model)
coupling_mtx <- flux_coupling_raptor(falcon_model)$coupled
sets <- get_list_of_sets_from_mtx(coupling_mtx)

Scripts to reproduce the figures in the paper can be found in the Pradhan2019Analysis repository
