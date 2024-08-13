from mindi import MindiTool
from termcolor import colored
import json
import numpy as np
import pandas as pd
from utils import parse_fasta
from scheduling import MiniBucketScheduler


configfile: 'config/config_IR.yaml'

buckets = config['buckets']
out = config['out']
mode = config['mode']

if mode not in {'IR', 'MR', 'DR'}:
    raise ValueError(f'Invalid mode {mode}.')

if mode == 'IR':
    extraction_type = 'Inverted Repeats'
    db_name = 'inverted_repeats'
elif mode == 'DR':
    extraction_type = 'Direct Repeats'
    db_name = 'direct_repeats'
elif mode == 'MR':
    extraction_type = 'Mirror Repeats'
    db_name = 'mirror_repeats'


def load_bucket(bucket):
    global buckets
    global out

    print(f"Reading schedule data from '{out}/schedule_{buckets}.json'...")
    with open(f"{out}/schedule_{buckets}.json", mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket)]

def extract_id(assembly):
    assembly = Path(assembly).name
    if '_' in assembly:
        return '_'.join(assembly.split('_')[:2])

    return assembly.split('.fna')[0]


def extract_name(accession):
    accession = Path(accession).name
    return accession.split(".fna")[0]

rule all:
    input:
        "%s/%s_completed/%s_db.parquet.snappy" % (out, mode, db_name)

rule schedule:
    output:
        "%s/schedule_%s.json" % (out, config["buckets"])
    params:
        out=Path(config["out"]).resolve(),
        buckets=config["buckets"],
        files=config["files"],
    run:
        print(colored(f"{len(params.files)} annotations detected.", "green"))
        assemblies = []
        with open(params.files, mode="r") as f:
            for line in f:
                line = line.strip()
                if line.count("\t") > 0:
                    line = line.split("\t")[0]

                assemblies.append(line)

        Path(params.out).mkdir(exist_ok=True)
        total_split = params.buckets

        mini_bucket_scheduler = MiniBucketScheduler()
        scheduled_files = mini_bucket_scheduler.schedule(assemblies, total_buckets=total_split)
        mini_bucket_scheduler.saveas(scheduled_files, output[0])

        # splitted_batches = {batch_id: job.tolist() for batch_id, job in enumerate(np.array_split(assemblies, total_split))}

        # with open(output[0], mode="w", encoding="UTF-8") as f:
        #    json.dump(splitted_batches, f, indent=4)


rule extractRepeats:
    input:
        '%s/schedule_%s.json' % (out, config["buckets"])
    output:
        touch("%s/%s_completed/bucket_{bucket}.%s.completed" % (out, mode, mode))
    params:
        out=Path(config['out']).resolve(),
        minrep=int(config['minrep']),
        max_spacer_length=int(config['max_spacer_length']),
        min_arm_length=int(config['min_arm_length']),
    run:
        # PREPARE DESTINATION DIRECTORIES
        Path(f"{params.out}/{mode}_completed").mkdir(exist_ok=True, parents=True)
        destination_dir = params.out.joinpath(f"{mode}_extracted_accessions")
        destination_dir.mkdir(exist_ok=True)

        
        
        bucket = load_bucket(wildcards.bucket)

        tempdir = Path(params.out).joinpath(f"{mode}_temp")
        destination_dir = params.out.joinpath(f"{mode}_extracted_accessions")

        mindi = MindiTool(tempdir=tempdir)


        def _extract_right_hand(left_arm: str, mode: str) -> str:
            if mode == 'IR':
                return ''.join({
                                'a': 't', 
                                't': 'a', 
                                'g': 'c', 
                                'c': 'g'
                                }[n] for n in left_arm)[::-1]
            elif mode == 'DR':
                return left_arm
            elif mode == 'MR':
                return left_arm[::-1]


        for accession in bucket:
            print(f'Processing accession {accession}.')
            accession_id = extract_name(accession)
            destination = destination_dir.joinpath(accession_id + f".{mode}.csv")

            if mode == 'IR':
                _ = mindi.extract_inverted(
                          accession, 
                          minrep=params.minrep,
                          min_arm_length=params.min_arm_length,
                          max_spacer_length=params.max_spacer_length
                    )
            elif mode == 'MR':
                _ = mindi.extract_mirror(
                        accession,
                        minrep=params.minrep,
                        min_arm_length=params.min_arm_length,
                        max_spacer_length=params.max_spacer_length
                        )
            elif mode == 'DR':
                raise NotImplementedYet('Direct Repeats extraction is not currently supported')
            

            repeat_df = mindi.to_dataframe()
            
            # accession_id_extracted = extract_id(mindi.fnp)
            # accession_id = extract_id(accession)
            # assert accession_id_extracted == accession_id

            total_chr = 0
            for seqID, sequence in parse_fasta(accession):

                temp = repeat_df[repeat_df['seqID'] == seqID]
                total_chr = 1

                total = temp.shape[0]
                total_validated = 0

                for _, row in temp.iterrows():

                    start = row['start']
                    end = row['end']
                    derived_seq = row['sequence']
                    arm_seq = row['sequenceOfArm']
                    arm_len = len(arm_seq)
                    spacer_seq = row['sequenceOfSpacer']
                    if spacer_seq == ".":
                        spacer_seq = ""

                    spacer_true_len = row['spacerLength']
                    arm_true_len = row['armLength']
                    seq_true_len = row['sequenceLength']

                    spacer_len = len(spacer_seq)

                    section = sequence[start: end]
                    
                    right_hand_side = _extract_right_hand(arm_seq, mode)
                    assert section == derived_seq == arm_seq + spacer_seq + right_hand_side, f"{accession} on {seqID} - start-end: {start}-{end}."
                    assert arm_len * 2 + spacer_len == len(section) == len(derived_seq) == arm_true_len * 2 + spacer_true_len == seq_true_len, f"{accession} on {seqID} - start-end: {start}-{end}."

                    total_validated += 1

                assert total_validated == total
            print(colored(f'Accession {accession} has passed all checks.', 'green'))
            assert total_chr == 1

            print(f'Accession {accession} has been processed succesfully.')
            mindi.moveto(destination)

rule reduceRepeats:
    input:
        expand("%s/%s_completed/bucket_{bucket}.%s.completed" % (out, mode, mode), bucket=range(buckets))
    output:
        "%s/%s_completed/%s_db.parquet.snappy" % (out, mode, db_name)
    params:
        out=Path(config["out"]).resolve(),
        buckets=config["buckets"]
    run:
        target = params.out.joinpath(f"{mode}_extracted_accessions")
        print(f"Fetching extracted accessions from target {target}.")
        files = [extracted_file for extracted_file in target.glob(f"*.{mode}.csv")]

        extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
        df_all = []
        empty_assemblies = []

        for file in files:
            print(f"Processing file {file}...")
            accession_id = extract_id(file)
            df = pd.read_table(file)
            if df.shape[0] == 0:
                empty_assemblies.append(file)
                continue

            df.loc[:, "#assembly_accession"] = accession_id
            df_all.append(df)
        
        with open(f"{params.out}/{mode}_completed/empty_accessions.{mode}.txt", mode="w", encoding="UTF-8") as f:
            for accession in empty_assemblies:
                f.write(str(accession) + "\n")

        df_all = pd.concat(df_all, axis=0)
        df_all.to_parquet(output[0], engine="fastparquet")
