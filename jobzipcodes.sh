#!/bin/bash
#SBATCH -J generate_layout_235x363_zipcodes_extra
#SBATCH --mem=64G
#SBATCH --time=04:00:00

source /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/.venv/bin/activate
cd /user/gent/446/vsc44685/ScratchVO_dir/layout_opt/Zipcodes
python generate_layout_235x363_zipcodes.py config_235x363_nosep_extra.json output_zipcodes_235x363_extra