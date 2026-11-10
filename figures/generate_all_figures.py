#!/usr/bin/env python3
"""
Master script to generate all method comparison figures
"""

import subprocess
import sys
from pathlib import Path

def check_matplotlib():
    """Check if matplotlib is installed"""
    try:
        import matplotlib
        print(f"✅ matplotlib {matplotlib.__version__} is installed")
        return True
    except ImportError:
        print("❌ matplotlib is not installed")
        return False

def install_matplotlib():
    """Install matplotlib if needed"""
    print("Installing matplotlib...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])
    print("✅ matplotlib installed")

def run_script(script_path):
    """Run a figure generation script"""
    print(f"\n{'='*60}")
    print(f"Running: {script_path.name}")
    print('='*60)
    try:
        subprocess.check_call([sys.executable, str(script_path)])
        print(f"✅ {script_path.name} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {script_path.name} failed with error code {e.returncode}")
        return False

def main():
    """Main execution"""
    print("="*60)
    print("METHOD COMPARISON FIGURES GENERATOR")
    print("="*60)
    
    # Check/install dependencies
    if not check_matplotlib():
        install_matplotlib()
    
    # List of figure scripts in order
    scripts = [
        'figure_1a_block_structure.py',
        'figure_1b_spatial_similarity.py',
        'figure_2a_quadtree_subdivision.py',
        'figure_2b_rll_permutation.py',
        'figure_3a_gray_code_principle.py',
        'figure_3b_spatial_comparison.py',
    ]
    
    current_dir = Path(__file__).parent
    success_count = 0
    failed_scripts = []
    
    for script_name in scripts:
        script_path = current_dir / script_name
        if script_path.exists():
            if run_script(script_path):
                success_count += 1
            else:
                failed_scripts.append(script_name)
        else:
            print(f"⚠️  Script not found: {script_name}")
            failed_scripts.append(script_name)
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print(f"✅ Successfully generated: {success_count}/{len(scripts)} figures")
    
    if failed_scripts:
        print(f"❌ Failed: {', '.join(failed_scripts)}")
    else:
        print("🎉 All figures generated successfully!")
    
    print("\nOutput files saved in current directory:")
    print("  - *.png (high resolution, 300 dpi)")
    print("  - *.pdf (vector format)")

if __name__ == "__main__":
    main()
