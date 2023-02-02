import argparse
import os

def process_fasta(file : str):
    lines = []
    with open(file) as fin:
        seq = ''
        for line in fin.readlines():
            if line.startswith('>'):
                if not seq == '':
                    lines.append(f"{seq}\n")
                lines.append(line)
                seq = ''
            else:
                seq+= line.strip()
    
    with open(file, 'w') as fout:
        fout.writelines(lines)

process_fasta("/drive/test/Altcar1.fasta")
'''    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", dest="input", help="Input directory containing .fasta files needing to be processed")

    args = parser.parse_args()

    if args.input == None:
        print("Input directory not specified...")
        exit()
    input_dir = args.input

    assert os.path.exists(input_dir), "input_dir could not be found"

    for file in os.listdir(input_dir):
        process_fasta(f"{input_dir}{file}")
'''