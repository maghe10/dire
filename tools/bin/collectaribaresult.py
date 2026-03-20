import os
import shutil

src_root = "/root/ariba"
dst_root = "/root/dire/data/Illumina/ariba"

# skapa destination om den inte finns
os.makedirs(dst_root, exist_ok=True)

for sample in os.listdir(src_root):

    src_dir = os.path.join(src_root, sample)

    if not os.path.isdir(src_dir):
        continue

    # --------------------------
    # 1. kopiera hela katalogen
    # --------------------------
    dst_dir = os.path.join(dst_root, sample)

    print(f"Copying directory {src_dir} -> {dst_dir}")

    shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)

    # --------------------------
    # 2. kopiera report.tsv
    # --------------------------
    report_src = os.path.join(src_dir, "report.tsv")
    report_dst = os.path.join(dst_root, f"{sample}.tsv")

    if os.path.exists(report_src):

        print(f"Copying report {report_src} -> {report_dst}")

        shutil.copy2(report_src, report_dst)

    else:
        print(f"No report.tsv found for {sample}")