from argparse import ArgumentParser
from pyfaidx import Faidx


def main():
    parser = ArgumentParser()
    parser.add_argument("fasta_file")
    args = parser.parse_args()
    fasta_index = Faidx(args.fasta_file)
    fasta_index.close() 


if __name__ == "__main__":
    main()
