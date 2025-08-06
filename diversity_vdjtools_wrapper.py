#!/users/miniconda3/envs/bin/python

"""
diversity_vdjtools_wrapper.py - vdjtools diversity calculate
usage: python diversity_vdjtools_wrapper.py -m metadata.txt -o output_dir
"""


import argparse
import subprocess
import os
import sys
from pathlib import Path



def validate_jar():
    """Verify if the JAR file exists"""
    script_dir = Path(__file__).parent.absolute()
    jar_path = script_dir / "tools" / "vdjtools-1.2.1.jar"
    
    if not jar_path.exists():
        print(f"\n[Error] Jar file not found: {jar_path}")
        print("Please perform the following actions:")
        print(f"1. Create the tools directory:: mkdir -p {script_dir/'tools'}")
        print(f"2. Put vdjtools-1.2.1. jar into {jar_path}")
        print("Or download from the following site")
        print("https://github.com/mikessh/vdjtools/releases/download/1.2.1/vdjtools-1.2.1.jar")
        sys.exit(1)
    return jar_path

def validate_resamples(value):
    """Validate resamples parameter"""
    try:
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError("Resamples must be a positive integer")
        return ivalue
    except ValueError:
        raise argparse.ArgumentTypeError("Resamples must be an integer")

def main():
    parser = argparse.ArgumentParser(description='Calculate TCR diversity metrics using vdjtools',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  # Default resamples (1,000,000)
  python %(prog)s -m metadata.txt -o results
  
  # Custom resamples (10,000,000)
  python %(prog)s -m metadata.txt -o results -x 10000000
""")
    # Required arguments
    parser.add_argument("-m", "--metadata", required=True, 
                        help="Input metadata file (columns: sample_id, file_path).")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory (default: current directory).")
    # Optional resamples argument
    parser.add_argument("-x", "--resamples", type=validate_resamples, default=1000000,
                      help="Number of resample reads for diversity estimation (default: 1000000)")
    args = parser.parse_args()

    # verify path
    jar_path = validate_jar()
    os.makedirs(args.outdir, exist_ok=True)
    # verify input
    if not os.path.exists(args.metadata):
        sys.exit(f"File not exists {args.metadata}")

    os.makedirs(args.outdir, exist_ok=True)

 
    

    # shell process
    cmd = [
        "java", "-jar", str(jar_path),
        "CalcDiversityStats",
        "-x", str(args.resamples),
        "-m", os.path.abspath(args.metadata),
        "."
    ]
    
    log_path = Path(args.outdir) / "vdjtools_calstat.log"
    result_file = Path(args.outdir) / "diversity.strict.resampled.txt"
    try:
        print(f"vdjtools analysis, log output to: {log_path}")
        
        with open(log_path, "w") as log:
            subprocess.run(cmd, cwd=args.outdir, check=True, 
                         stdout=log, stderr=subprocess.STDOUT)
        
        if result_file.exists():
            print("Successfully generated result file.")
            print(f"result file path: {result_file}")
        else:
            raise FileNotFoundError(f"No result file generated, please check the log: {log_path}")
            
    except subprocess.CalledProcessError as e:
        print(f"[Error] vdjtools execution failed, return code: {e.returncode}")
        print(f"please check the log: {log_path}")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {str(e)}")
        sys.exit(1)
            
    

if __name__ == "__main__":
    main()
    
