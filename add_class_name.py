"""Creates a csv file wtih species names and class names
Input is a file with fasta sequences, where the names of the species are the IDs in JGI
It looks up the species in JGI taxonomy server

Data file can be downloaded from JGI:
https://mycocosm.jgi.doe.gov/ext-api/mycocosm/catalog/download-group?flt=&seq=all&grp=fungi&srt=released&ord=asc

Note: Some species are incorrect due to synonyms in names, or new names
Names that can't be found can often be double-checked in the NCBI taxonomy server
"""
import argparse
import re
import requests

# regex parser for data file
PARSER = re.compile(r'"{3}(\w* \w*)')


def find_class_names(file: str, output: str, data_file: str):
    """Finds class names for species using JGI ids found in fasta file

    Args:
        file (str): Input .fasta file
        output (str): Output file containing all species with their respective classes in csv format
        data_file (str): Data file containing IDs and species names from JGI
    """
    tax = {}
    long_name = {}
    names = ""
    with open(file, encoding="utf-8") as file_cont:
        for line in file_cont.readlines():
            if line.startswith(">"):
                name = line.split(">")[1].strip()
                tax[name] = "not_found"
                with open(data_file, errors="ignore", encoding="utf-8") as ref:
                    for refline in ref.readlines():
                        id_name = refline.split(",")[2]
                        l_name = refline.split(",")[1]

                        if name == id_name:
                            l_name = PARSER.match(l_name).group(1)
                            l_name = "_".join(l_name.split(" "))
                            long_name[name] = l_name
                            names += f"{l_name},"

    url = "https://taxonomy.jgi-psf.org/POST"
    res = requests.post(url, data="name/" + names, timeout=120)

    res_json = res.json()

    with open(output, "w", encoding="utf-8") as output_file:
        for entries in long_name.items():
            key = entries[0]
            try:
                class_name = res_json[long_name[key]]["class"]["name"]
                output_file.write(f"{key},{class_name}\n")
            except KeyError:
                print(f"{key}|{long_name[key]} not found in taxonomy...")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest="file", help="Fasta file")
    parser.add_argument("-o", "--out", dest="output", help="Output destination")
    parser.add_argument("-d", "--data", dest="data", help="Data file from JGI")

    args = parser.parse_args(
        [
            "-f",
            "output_fastas_per_gene_yanfang/concat.fasta",
            "-o",
            "classes_yanfang.txt",
            "-d",
            "fungi.txt",
        ]
    )

    if args.file is None:
        print("Input file not specified")
        exit()
    if args.output is None:
        print("Output file not specified")
        exit()
    if args.data is None:
        print("Data file not specified")

    __file = args.file
    __output = args.output
    __data_file = args.data
    find_class_names(__file, __output, __data_file)
