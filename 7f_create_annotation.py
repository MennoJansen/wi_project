file = "classes.txt"

colors = {
    "Eurotiomycetes": "#ce7e00",
    "Dothideomycetes": "#8fce00",
    "Lecanoromycetes": "#2986cc",
    "Leotiomycetes": "#6a329f",
    "Sordariomycetes": "#6dbbcc",
    "Pezizomycetes": "#f44336",
    "Xylonomycetes": "#990000",
    "Orbiliomycetes": "#744700",
    "Peziz_Incertae": "#000000",
    "Taph_Incertae": "#7bffbd",
    "Saccharomycetes": "#ecff0e",
    "Neolectomycetes": "#cc00ff",
    "Pneumocystidomycetes": "#8bc34a",
    "Schizosaccharomycetes": "#fad0c3",
    "Taphrinomycetes": "#d4cfb",
    "Coniocybomycetes": "#bedadc",
}

with open("annotation.txt", "w") as a:
    a.write("DATASET_COLORSTRIP\n")
    a.write("SEPARATOR COMMA\n")
    a.write("DATASET_LABEL,classes\n")
    a.write("COLOR,#ff0000\n")
    a.write("COLOR_BRANCHES,1\n\n")

    a.write("DATA\n\n")
    with open("classes.txt") as c:
        for line in c.readlines():
            split_line = line.split(",")
            name = split_line[0]
            color = colors[split_line[1].strip()]
            a.write(f"{name},{color}\n")
