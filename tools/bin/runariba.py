import os
import subprocess
import argparse

# ---------------------------
# Sample lists
# ---------------------------
samplesA = ["2","3","4","5","6","8"]
samplesB = ["9","11","12","13","14","18","19","22"]
samplesC = ["24","25","26","28","29","34","35","36"]
samplesD = ["39","40","43","44","45","47","48","50"]
samplesE = ["51","52","56","57","59","61","62","64"]
samplesF = ["66","67","69","70","72","76","80","82"]
samplesG = ["83","84","86","87","88","89","91","93"]
samplesH = ["94","96","97","98","101","102","103","105"]
samplesI = ["109","111","112","113","120"]
samplesJ = ["17","20","23","27","32","37","38","41","46"]
samplesK = ["49","53","54","58","75","95","100","115","116"]
samplesL = ["21","63","77","78","114"]
samplesM = ["121","122","123","124","55B","60B","65B","79B","85B","110B"]
samplesN = ["125"]

allSamples = (
    samplesA + samplesB + samplesC + samplesD +
    samplesE + samplesF + samplesG + samplesH +
    samplesI + samplesJ + samplesK + samplesL +
    samplesM + samplesN
)

# ---------------------------
# Argument parsing
# ---------------------------
parser = argparse.ArgumentParser(description="Run ARIBA on all samples")

parser.add_argument(
    "--reads_dir",
    default="/root/dire/data/Illumina/trimmed",
    help="Directory with trimmed reads"
)

parser.add_argument(
    "--db",
    default="/root/resfinder_ariba_db",
    help="ARIBA database directory"
)

parser.add_argument(
    "--outdir",
    default="/root/ariba",
    help="Output directory"
)

parser.add_argument(
    "--force",
    action="store_true",
    help="Run ARIBA even if output already exists"
)

args = parser.parse_args()

reads_dir = args.reads_dir
db = args.db
outdir = args.outdir
force = args.force

os.makedirs(outdir, exist_ok=True)

# ---------------------------
# Run ARIBA
# ---------------------------
for sample in allSamples:
    r1 = os.path.join(reads_dir, f"sample{sample}_1_val_1.fq.gz")
    r2 = os.path.join(reads_dir, f"sample{sample}_2_val_2.fq.gz")

    sample_out = os.path.join(outdir, f"sample{sample}")
    report_file = os.path.join(sample_out, "report.tsv")

    if not os.path.exists(r1) or not os.path.exists(r2):
        print(f"Skipping sample{sample}: reads not found")
        continue

    if os.path.exists(report_file) and not force:
        print(f"Skipping sample{sample}: output already exists ({report_file})")
        continue


    cmd = [
        "ariba",
        "run",
        db,
        r1,
        r2,
        sample_out
    ]

    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)

print("All samples processed.")