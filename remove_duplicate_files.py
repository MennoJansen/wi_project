import shutil
import os

dir = "/drive/gff/"

for file in os.listdir(dir):
    if file.endswith(".gff3"):
        try:
            os.remove(f"{dir}{file[:-1]}")
            print("Deleted file")
        except:
            print("File could not be found")
        