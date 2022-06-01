

if __name__ == "__main__":
    his_annotated_file = "/Users/dhthutrang/Documents/BIOINFO/Episplicing/ENCODE_episplicing/majiq/temp.tsv"
    with open(his_annotated_file) as reader:
        lines = [line.rstrip().split('\t') for line in reader]

    LSV_info = defaultdict()
    for i, line in enumerate(lines):