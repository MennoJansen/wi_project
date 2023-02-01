import sys
import os
import subprocess

cdir = os.getcwd()
pdir = cdir + "/proteins/"
hmm_file = cdir + "/fungal_models.hmm"
odir = cdir + "/hmm_output2/"
counter = 1
total = 982

#check if file is ok
def check_file(file):
    with open(file, 'rb') as f:
        try:
            f.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
        except OSError:
            f.seek(0)
        last_line = f.readline().decode()
        if last_line.startswith("# [ok]"):
            return True
        return False

failed_hmms = []

for root, dirs, files in os.walk(pdir):
    for file_name in files:
        name = file_name.split(".")[0]
        print("hmmsearch for " + file_name + " (" + str(counter) + "/" + str(total) + ")")
        counter += 1
        #print(["hmmsearch", "-tblout", "hmm_out_" + name + ".txt", hmm_file, file_name])
        if not os.path.exists(odir + "hmm_out_" + name + ".txt") or not check_file(odir + "hmm_out_" + name + ".txt"):
            subprocess.run(["hmmsearch", "--cpu", str(12), "-E", str(0.0000000001), "--tblout", odir + "hmm_out_" + name + ".txt", hmm_file, pdir + file_name], stdout=subprocess.DEVNULL)
        if not os.path.exists(odir + "hmm_out_" + name + ".txt"):
            failed_hmms.append(name)

print("Hmmsearch failed for these entries:")
print(failed_hmms)


