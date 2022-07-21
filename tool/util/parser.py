import argparse

def get_arguments() :
    parser = argparse.ArgumentParser(description = "")
    parser.add_argument('-i', '--input', type = str, nargs = '*', help = "Input file(s)")
    parser.add_argument('-o', '--output', type = str, nargs = '*', help = "Output file(s)")
    args = parser.parse_args()
    return args