#To generate the layour with 10 block rows extra for decoy-based precision ()

source /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/.venv/bin/activate
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/HSDB/bar1_zipcode
python generate_block_barcodes.py --config block_zipcode_config_253x363_extra.json --out-dir 235x363_extra
cd ../layout_generator
python dual_barcode_hierarchical.py \
    --array-rows 295 --array-cols 363 \
    --bar1 ../bar1_zipcode/235x363_extra/bar1_barcodes_28nt.txt \
    --bar2 ../bar2_library/spot_barcodes_8nt.txt \
    --out-dir ../outputs/295x363_extra \
    --seed 42

# For zipcodes
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Zipcodes
python generate_layout_235x363_zipcodes.py config_235x363_nosep_extra.json output_zipcodes_235x363_extra

#Another idea is not to apply zipcodes strategy to blocks but instead use a simple spatial gradient strategy to generate the barcodes for the blocks.
