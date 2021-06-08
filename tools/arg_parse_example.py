import argparse


arg_parser = argparse.ArgumentParser(description=" describe argparse")

arg_parser.add_argument("bedfile", help='bed file w/ full path')

args = arg_parser.parse_args()

F = args.bedfile
