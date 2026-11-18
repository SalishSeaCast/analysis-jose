
#!/bin/bash
# Usage:
#   source stage.sh 01jan24
#   source stage.sh 01jan24 02jan24 03jan24

SRC=/home/jvalenti/scratch/Parcels/forcing
DST=$SLURM_TMPDIR/forcing

# Check prerequisites
if [[ -z "$SLURM_TMPDIR" ]]; then
    echo "Error: SLURM_TMPDIR is not set."
    return 1
fi

if ! command -v nccopy >/dev/null 2>&1; then
    echo "Error: nccopy not found. Load the NetCDF module first."
    return 1
fi

if [[ $# -eq 0 ]]; then
    echo "Error: No dates provided."
    return 1
fi

echo "SRC: $SRC"
echo "DST: $DST"
echo "Arguments: $@"

# Stage requested days
for d in "$@"; do
    mkdir -p "$DST/SalishSeaCast/$d" \
             "$DST/SHEM/$d" \
             "$DST/WW3/$d"

    rsync -rt \
        --include='SalishSea_1h_????????_????????_grid_[TUVW].nc' \
        --exclude='*' \
        "$SRC/SalishSeaCast/$d/" "$DST/SalishSeaCast/$d/" || return 1

    rsync -rt \
        --include='SHEM_1h_*_biol_T.nc' \
        --include='SHEM_1d_*_prod_T.nc' \
        --include='SHEM_1d_*_graz_T.nc' \
        --exclude='*' \
        "$SRC/SHEM/$d/" "$DST/SHEM/$d/" || return 1

    echo "Staging WW3 for $d"
    echo "Source: $SRC/WW3/$d/"
    echo "Destination: $DST/WW3/$d/"

    mkdir -p "$DST/WW3/$d"

    rsync -av "$SRC/WW3/$d/" "$DST/WW3/$d/" || {
        echo "WW3 rsync failed for $d"
        return 1
    }
    #rsync -rt "$SRC/WW3/$d/" "$DST/WW3/$d/" || return 1
done

# Uncompress all staged NetCDF files
echo "Uncompressing staged NetCDF files..."

while IFS= read -r -d '' file; do
    tmp="${file%.nc}_uncompressed.nc"

    echo "Uncompressing: $file"

    if nccopy -d 0 "$file" "$tmp"; then
        if mv -f "$tmp" "$file"; then
            echo "Done: $file"
        else
            echo "Error replacing: $file"
            rm -f "$tmp"
            return 1
        fi
    else
        echo "Error uncompressing: $file"
        rm -f "$tmp"
        return 1
    fi

done < <(find "$DST" -type f -name '*.nc' -print0)

# Report staged data size
du -sh "$DST"/*/*

# Set forcing directory only after successful completion
export FORCING_DIR="$DST"

echo "FORCING_DIR set to $FORCING_DIR"
