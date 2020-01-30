#!/bin/env python3
import argparse
import subprocess
import sys
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", default='dengue.fa', help="specify fasta file to test (must be located in tests/test)")
    parser.add_argument("--big", action='store_true', help="use 64-bit")
    parser.add_argument("--sa", action='store_true', help="test full SA output")
    parser.add_argument("--rssa", action='store_true', help="test full RLSA output")
    args = parser.parse_args()

    if args.sa and args.rssa:
        print("cannot activate both -sa and -rssa at the same time!")
        exit(1)

    pfbwt = "./pfbwt-f"
    fasta_loc = os.path.join("./tests/test/", args.fasta)
    options = []
    if args.sa:
        options.append("-s")
    elif args.rssa:
        options.append("-r")
    if args.big:
        options.append("-m")

    pfbwt_cmd = [pfbwt, fasta_loc]

    # run the program and return error if it fails
    try:
        subprocess.run(pfbwt_cmd, check=True, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("ERROR: pfbwt exited with non-zero status")
        exit(1)

    # compare output if program run was succcessful
    truth_loc = os.path.join("./tests/truth/", args.fasta)

    def run_diff(ext):
        diff = ["diff", "-q"]
        try:
            subprocess.run(diff + [truth_loc + ext, fasta_loc + ext], check=True)
        except subprocess.CalledProcessError:
            print("ERROR: " + ext + " files differ!") 
            return False
        return True

    exts = ['.bwt']
    if args.sa:
        exts += ['.sa']
    elif args.rssa:
        exts += ['.ssa', '.esa']

    for ext in exts:
        if run_diff(ext):
            print("SUCCESS: " + ext)
        else:
            print("FAIL")
            exit(1)
