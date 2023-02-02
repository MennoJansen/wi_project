import requests
import argparse
import re

fungi = "fungi.txt"
p = re.compile('"{3}(\w* \w*)')

def find_class_names(file : str, output : str):
    tax = {}
    long_name = {}
    names = ""
    with open(file) as f:
        for line in f.readlines():
            if line.startswith(">"):
                    name = line.split('>')[1].strip()
                    tax[name] = "not_found"
                    with open(fungi, errors='ignore') as ref:
                         for refline in ref.readlines():
                            id_name = refline.split(',')[2]
                            l_name = refline.split(',')[1]

                            if name == id_name:
                                l_name = p.match(l_name).group(1)
                                l_name = "_".join(l_name.split(" "))
                                long_name[name] = l_name
                                names += f"{l_name},"
    
    url = "https://taxonomy.jgi-psf.org/POST"
    res = requests.post(url, data="name/" + names)

    res_json = res.json()

    with open(output, 'w') as o:         
        for key in long_name:
            try:
                class_name = res_json[long_name[key]]['class']['name']
                o.write(f"{key},{class_name}\n")
            except:
                print(f"{key}|{long_name[key]} not found in taxonomy...")         
    print("test")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file", dest="file", help="Fasta file")
    parser.add_argument("-o", "--out", dest="output", help="Output destination")

    args = parser.parse_args()
    file = "output_fastas_per_gene_v3/cat_trimal.fasta"
    output = "classes.txt"
    find_class_names(file, output)

