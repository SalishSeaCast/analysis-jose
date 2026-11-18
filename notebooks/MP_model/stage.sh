#!/bin/bash
# usage: source stage.sh 01jan24 [02jan24 ...]
SRC=/home/jvalenti/scratch/Parcels/forcing
DST=$SLURM_TMPDIR/forcing
for d in "$@"; do
  mkdir -p $DST/SalishSeaCast/$d $DST/SHEM/$d $DST/WW3/$d
  rsync -rt --include='SalishSea_1h_????????_????????_grid_[TUVW].nc' --exclude='*' \
        $SRC/SalishSeaCast/$d/ $DST/SalishSeaCast/$d/
  rsync -rt --include='SHEM_1h_*_biol_T.nc' --include='SHEM_1d_*_prod_T.nc' \
        --include='SHEM_1d_*_graz_T.nc' --exclude='*' $SRC/SHEM/$d/ $DST/SHEM/$d/
  rsync -rt $SRC/WW3/$d/ $DST/WW3/$d/
done
du -sh $DST/*/*
export FORCING_DIR=$DST