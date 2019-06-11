#!/bin/bash

# run Zeo++ with modified extended output, printing information about PLD and cell parameters
network  -xyz ./coree_BAZJET_clean.xyz -cellparams ./coree_BAZJET_clean.pbc -r ./UFF.rad -resex ./coree_BAZJET_clean.resex ./coree_BAZJET_clean.cif > 0_verbose.txt

# parse output for use in geometric overlap method
echo coree_BAZJET_clean $(head -n 3 coree_BAZJET_clean.resex | tail -n 1) $(head -n 1 coree_BAZJET_clean.pbc) > 0_parsed_output.txt
