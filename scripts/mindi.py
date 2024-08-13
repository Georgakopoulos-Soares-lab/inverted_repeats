# Mindi: A Non B-DNA extractor tool for data analysis

__author__ = "Nikol Chantzi"
__version__ = "1.0.1"
__email__ = "nmc6088@psu.edu"


### > IMPORTS BEGIN

import subprocess
import pandas as pd
import csv
import tempfile
from pathlib import Path
from termcolor import colored
import functools
import re
from typing import ClassVar, Optional, Callable, Iterator
from dotenv import load_dotenv
# from tailhunter import hunt_tail
import os
import stat
import gzip
from utils import parse_fasta
from dataclasses import dataclass
import shutil
import uuid
import sys

csv.field_size_limit(sys.maxsize)

### < IMPORTS END

@dataclass
class Defaults:

    minrep: int = 8
    spacer_length: int = 8
    min_arm_length: int = 8
    max_sru: int = 9
    min_consensus_repeats: int = 3

extract_id = lambda accession: '_'.join(Path(accession).name.split('_')[:2])
extract_name = lambda accession: Path(accession).name.split('.fna')[0]

def nonbdna_register(mode):
    def nonbdna(func):
        @functools.wraps(func)
        def wrapper(self, *args, **kwargs):
            global extract_name

            self.reset()
            if 'accession' in kwargs:
                accession = kwargs['accession']
            else:
                accession = args[0]

            accession = Path(accession).resolve()
            if not accession.is_file():
                raise FileNotFoundError(f'Unable to detect accession {accession}.')

            accession_name = extract_name(accession)
            accession_id = extract_id(accession)

            # accession_tmp_dir = self.tempdir.joinpath(accession_name)
            accession_tmp_dir = tempfile.TemporaryDirectory(prefix=f"{accession_id}.{mode}.")
            accession_tmp_dir_path = self.tempdir.joinpath(accession_tmp_dir.name)

            # if accession_tmp_dir.is_dir():
            #    shutil.rmtree(accession_tmp_dir)

            cur_dir = os.getcwd()
            
            if 'minrep' in kwargs:
                minrep = kwargs['minrep']
            elif len(args) > 1:
                minrep = args[1]
            else:
                # fallback to default
                minrep = int(Defaults.minrep)

            # handle zipped files
            if accession.name.endswith(".gz"):

                # gzip compression strategy
                with tempfile.NamedTemporaryFile(
                                                 dir=accession_tmp_dir_path, 
                                                 delete=False, 
                                                 suffix='.fna') as unzipped_tmp:

                    with gzip.open(accession, 'rb') as handler:
                        unzipped_tmp.write(handler.read())

                    accession = Path(unzipped_tmp.name).name
            
            elif accession.name.endswith(".xz"):

                # different compression strategy
                raise NotImplementedYet()

            os.chdir(accession_tmp_dir_path)

            out_tsv = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.tsv')
            out_gff = accession_tmp_dir_path.joinpath(accession_name + f'_{mode}.gff')
            rand_accession_name = str(uuid.uuid4())
            
            match mode:
                case 'IR':
                    command = f"{self.nonBDNA} -seq {accession} -out {rand_accession_name} -minIRrep {minrep} -skipAPR -skipSTR -skipMR -skipDR -skipGQ -skipZ -skipSlipped -skipCruciform -skipTriplex -skipWGET"
                case 'MR':
                    command = f"{self.nonBDNA} -seq {accession} -out {rand_accession_name} -minMRrep {minrep} -skipAPR -skipSTR -skipIR -skipDR -skipGQ -skipZ -skipSlipped -skipCruciform -skipTriplex -skipWGET"
                case 'DR':
                    raise NotImplementedYet('Direct Repeats Finder is not currently supported by this version of Mindi.')
                case 'STR':
                    raise NotImplementedYet('Short Tandem Repeats are not supported yet!')
                case _:
                    raise ValueError(f'Unknown mode {mode}.')
            

            _ = subprocess.run(command, shell=True, 
                                        check=True,
                                        stdout=subprocess.DEVNULL,
                                        stderr=subprocess.DEVNULL,
                                )

            # check operation was succesful
            if not Path(rand_accession_name + f"_{mode}.tsv").is_file():
                raise FileNotFoundError(f"Failed to extract {mode} for {accession}.")

            shutil.move(rand_accession_name + f"_{mode}.tsv", accession_name + f"_{mode}.tsv")
            shutil.move(rand_accession_name + f"_{mode}.gff", accession_name + f"_{mode}.gff")

            destination = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
            if destination.is_file():
                os.remove(destination)
            
            shutil.move(out_tsv, destination)

            # remove redundant files
            os.unlink(out_gff)
            os.chdir(cur_dir)

            self.fn = self.tempdir.joinpath(accession_name + f'_{mode}.tsv')
            result = func(self, *args, **kwargs)

            accession_tmp_dir.cleanup()

            return result
        return wrapper
    return nonbdna


class MindiTool:
    
    MINDI_FIELDS: ClassVar[list[str]] = [
                                         "seqID", 
                                         "start", 
                                         "end", 
                                         "sequenceOfArm", 
                                         "sequenceOfSpacer", 
                                         "sequence", 
                                         "armLength", 
                                         "spacerLength", 
                                         "length", 
                                         "arm_a", 
                                         "arm_c", 
                                         "arm_g", 
                                         "arm_t"
                                         ]

    REGEX_FIELDS: ClassVar[list[str]] = [
                    "seqID",
                    "start",
                    "end",
                    "sequence",
                    "strand",
                    "gc_content",
                    "length",
                    "stackerLength",
                    "spacerLength",
                    ]

    TAIL_FIELDS: ClassVar[list[str]] = [
                            "seqID", 
                            "start", 
                            "end", 
                            "sequence", 
                            "length", 
                            "consensus", 
                            "sru", 
                            "consensus_repeats"
                            ]

    FRAME_FIELDS: ClassVar[list[str]] = [
                             "seqID", 
                             "start", 
                             "end", 
                             "sequenceOfArm", 
                             "sequenceOfSpacer", 
                             "sequence", 
                             "armLength", 
                             "spacerLength", 
                             "sequenceLength",
                             "arm_a",
                             "arm_g",
                             "arm_c",
                             "arm_t",
                             "composition"
                          ]

    HDNA_max_at_content: ClassVar[float] = 0.8
    HDNA_min_pyrine: ClassVar[float] = 0.9
    HDNA_max_pyrine: ClassVar[float] = 0.9

    def __init__(self, tempdir: Optional[str] = None,
                       nonBDNA: Optional[str] = None,
                 RTRF: Optional[str] = None) -> None:
        load_dotenv()
        if nonBDNA is None:
            nonBDNA = os.getenv('nonBDNA')

        if RTRF is None:
            RTRF = os.getenv('RTRF')

        if tempdir is None:
            self.tempdir = Path().cwd()
        else:
            self.tempdir = Path(tempdir).resolve()
            self.tempdir.mkdir(exist_ok=True)

        self.nonBDNA = nonBDNA
        self.RTRF = RTRF
        self.cur_mode = None
        self.fn = None
        self.fnp = None

        if not Path(self.nonBDNA).is_file():
            raise FileNotFoundError(f"Invalid executable path {self.nonBDNA}.")

    def reset(self) -> None:
        self.fn = None
        self.fnp = None
        self.cur_mode = None

    def moveto(self, dest: os.PathLike[str]) -> None:
        shutil.move(self.fnp, dest)
    
    def to_dataframe(self, usecols: bool = True) -> pd.DataFrame:
        # STILL EXPERIMENTAL
        if self.cur_mode == 'RE' or self.cur_mode == "TAIL":
            tmp_file = self.fn
        else:
            tmp_file = self.fnp

        if usecols:
            mindi_frame = pd.read_table(tmp_file)

            if self.cur_mode == "IR" or self.cur_mode == "DI" or self.cur_mode == "MR":
                mindi_frame.loc[:, "sequenceOfSpacer"] = mindi_frame["sequenceOfSpacer"].fillna(".")

        else:
            mindi_frame = pd.read_table(tmp_file, header=None, skiprows=1)

        return mindi_frame

    def set_tempdir(self, tempdir: os.PathLike[str]) -> None:
        self.tempdir = Path(tempdir).resolve()
        self.tempdir.mkdir(exist_ok=True)
    
    @nonbdna_register(mode='IR')
    def extract_inverted(self, accession: os.PathLike[str], 
                        minrep: int = 8, 
                        min_arm_length: Optional[int] = None, 
                        max_spacer_length: Optional[int] = None,
               ) -> "MindiTool":

       mindi_table = self.process_table(
                    min_arm_length=min_arm_length, 
                    max_spacer_length=max_spacer_length
                    )
       self.cur_mode = 'IR'
       return self

    @nonbdna_register(mode='MR') 
    def extract_mirror(self, accession: os.PathLike[str], 
                             minrep: int = 8, 
                             min_arm_length: Optional[int] = None, 
                             max_spacer_length: Optional[int] = None,
                            ) -> "MindiTool":
        mindi_table = self.process_table(min_arm_length=min_arm_length, 
                                         max_spacer_length=max_spacer_length)
        self.cur_mode = 'MR'
        return self

    def filter_hdna(self) -> pd.DataFrame:
        # assert current extraction is of mode 'Mirror Repeat'
        if self.cur_mode != 'MR':
            raise ValueError(f"Invalid mode '{mode}' for hdna filtering.")

        mindi_table = pd.read_table(self.fnp)
        mindi_table.loc[:, "pyrine"] = mindi_table["sequence"].str.count("a|g")
        mindi_table.loc[:, "pyrimidine"] = mindi_table["sequence"].str.count("c|t")
        mindi_table.loc[:, "at_content"] = mindi_table["sequence"].str.count("a|t")

        mindi_table = mindi_table[(mindi_table["at_content"] <= 0.8) & ((mindi_table["pyrimidine"] >= 0.9) | (mindi_table["pyrine"] >= 0.9))]

        return mindi_table

    def extract_regex(self, accession: os.PathLike[str], 
                            stacker: str = "g", 
                            minrep: int = 3,
                            multiplicity: int = 3) -> "MindiTool":

        self.cur_mode = "RE"
        regex = RegexExtractor(
                               stacker=stacker, 
                               minrep=minrep, 
                               multiplicity=multiplicity
                               )


        # global extract_name--cluster-config config/cluster_overlap.yaml --cores $j
        # accession_name = extract_name(accession)
        # accession_tmp_dir = self.tempdir.joinpath(accession_name)

        with tempfile.NamedTemporaryFile(dir=self.tempdir, delete=False, suffix=".re", mode="w") as file:
            dict_writer = csv.DictWriter(file, delimiter="\t", fieldnames=MindiTool.REGEX_FIELDS)
            dict_writer.writeheader()

            for table in regex.parse_regex(accession):
                for row in table:
                    dict_writer.writerow(row)

            self.fn = file.name

        return self


    def cleanup(self) -> None:
        if self.fnp and Path(self.fnp).is_file():
            os.remove(self.fnp)
        else:
            print(colored(f"WARNING! Processed file {self.fnp} does not exist.", "red"))


    def filter_cruciform(self) -> pd.DataFrame:
        # assert current extraction is of mode 'Inverted Repeat'
        if self.cur_mode != 'IR':
            raise ValueError(f"Invalid mode '{mode}' for cruciform filtering.")

        mindi_table = pd.read_table(self.fnp)
        return mindi_table


    def extract_tails(self, accession: os.PathLike[str], minrepeat: int = 8) -> "MindiTool":
        self.cur_mode = "TAIL"
        with tempfile.NamedTemporaryFile(dir=self.tempdir, 
                                         delete=False, 
                                         suffix=".tail", 
                                         mode="w") as file:
            dict_writer = csv.DictWriter(
                                         file, 
                                         delimiter="\t", 
                                         fieldnames=MindiTool.TAIL_FIELDS
                                         )
            dict_writer.writeheader()
            for seqID, sequence in parse_fasta(accession):
                for mononucleotide_tail in hunt_tail(sequence=sequence, 
                                                     minrepeat=minrepeat,
                                                     seqID=seqID
                                                     ):
                    dict_writer.writerow(mononucleotide_tail)

            self.fn = file.name
        return self


    @nonbdna_register(mode='STR')
    def extract_tandem(self, accession: os.PathLike[str], 
                             min_sru: Optional[int] = None) -> "MindiTool":
        raise NotImplementedYet()
        self._process_RTRF()


    def _process_RTRF(self) -> "MindiTool":
        raise NotImplementedYet()


    @nonbdna_register(mode='DR')
    def extract_direct(self, accession: os.PathLike[str], 
                    minrep: int = 8,
                    min_arm_length: Optional[int] = None,
                    max_spacer_length: Optional[int] = None
                ) -> "MindiTool":
        raise NotImplementedYet()

    def process_table(self, min_arm_length: Optional[int] = None, 
                            max_spacer_length: Optional[int] = None, 
                            ) -> "MindiTool":
        to_drop = ["Source",
                   "Type",
                   "Score",
                   "Strand",
                   "Subset",
                   "Permutations",
                   "Sequence",
                   "Start",
                   "Stop"]
        tmp_file = tempfile.NamedTemporaryFile(
                                               dir=self.tempdir, 
                                               prefix=str(self.fn).replace(".tsv", "") + ".",
                                               suffix=".tsv", 
                                               delete=False, 
                                               mode="w"
                                            )
        self.fnp = tmp_file.name
        nucleotides = {'a', 'g', 'c', 't'}

        with tmp_file as fout:
            fout_writer = csv.DictWriter(
                                         fout, 
                                         delimiter="\t", 
                                         fieldnames=MindiTool.FRAME_FIELDS
                                         )

            fout_writer.writeheader()

            with open(self.fn, mode="r", encoding="UTF-8") as fh:
                reader = csv.DictReader(fh, delimiter="\t")

                for row in reader:
                    arm_length = int(row['Repeat'])
                    spacer_length = int(row['Spacer'])
                    sequence = row['Sequence']

                    if (isinstance(min_arm_length, int) and arm_length < min_arm_length) or (isinstance(max_spacer_length, int) and spacer_length > max_spacer_length):
                        continue

                    start = int(row['Start']) - 1
                    end = int(row['Stop'])
                    sequence_length = int(row['Length'])
                    total_coordinate_length = end - start

                    if sequence_length < total_coordinate_length:
                        sequence = sequence[:sequence_length]
                        end = end - (total_coordinate_length - sequence_length)
                        print(colored(f'Invalid sequence length detected for {self.fn} on (start,end)=({start},{end}) with sequence {sequence}.', 'red'))
                        raise ValueError(f'Invalid sequence length detected for {self.fn}.')


                    # find sequence of arm
                    repeat = int(row['Repeat'])
                    sequence_of_arm = sequence[:repeat]

                    # find spacer
                    del row['Spacer']
                    true_spacer_length = sequence_length - 2 * repeat 
                    right_arm = sequence[repeat+true_spacer_length:]
                

                    if any(n not in nucleotides for n in right_arm) or any(n not in nucleotides for n in sequence_of_arm):
                        continue

                    # spacer = sequence[repeat:repeat+spacer_length]
                    true_spacer = sequence[repeat:repeat+true_spacer_length]
                    if len(true_spacer) == 0:
                        true_spacer = "."

                    # process composition
                    composition = re.search(r"(\d+)A/(\d+)C/(\d+)G/(\d+)T", row['Composition'])

                    a_content = composition.group(1)
                    c_content = composition.group(2)
                    g_content = composition.group(3)
                    t_content = composition.group(4)
                    
                    # skip maximum spacer length
                    if (isinstance(max_spacer_length, int) and true_spacer_length > max_spacer_length):

                        # invalid record?
                        continue

                    row.update({
                            "start": start,
                            "end": end,
                            "sequenceOfArm": sequence_of_arm,
                            "sequenceOfSpacer": true_spacer,
                            "spacer": true_spacer_length,
                            "sequence": sequence,
                            "arm_a": a_content,
                            "arm_c": c_content,
                            "arm_g": g_content,
                            "arm_t": t_content,
                        })

                    for col in to_drop:
                        del row[col]

                    row["sequenceLength"] = row.pop("Length")
                    row["seqID"] = row.pop("Sequence_name")
                    row["armLength"] = row.pop("Repeat")
                    row["spacerLength"] = row.pop("spacer")
                    row["composition"] = row.pop("Composition")
                    fout_writer.writerow(row)

        return self


class RegexExtractor:

    def __init__(self, stacker: str = "g", multiplicity: int = 3, minrep: int = 3) -> None:
        self.stacker = stacker
        self.multiplicity = multiplicity
        self.complementary_stacker = RegexExtractor.invert(self.stacker)
        self.minrep = minrep
        self.motifs = [
                       "%s{%s,}[agct]{1,7}" % (self.stacker, self.minrep) * self.multiplicity + "%s{%s,}" % (self.stacker, self.minrep),
                       "%s{%s,}[agct]{1,7}" % (self.complementary_stacker, self.minrep) * self.multiplicity + "%s{%s,}" % ( self.complementary_stacker, self.minrep)
                       ]

        self.strict_motif = "%s" % (self.stacker * self.minrep) + "[agct]{1,7}"
        self.strict_motif = self.strict_motif * self.multiplicity + self.strict_motif

        # self.complementary_strict_motif 

    
    @staticmethod
    def complement(nucleotide: str) -> str:
        match nucleotide:
            case 'a':
                return 't'
            case 't':
                return 'a'
            case 'g':
                return 'c'
            case 'c':
                return 'g'
            case 'n':
                return 'n'
            case _:
                raise ValueError(f'Invalid nucleotide {nucleotide}.')

    @staticmethod
    def invert(kmer: str) -> str:
        return ''.join(RegexExtractor.complement(n) for n in kmer)[::-1]


    def parse_regex(self, accession: os.PathLike[str]) -> Iterator[list[dict]]:
        for seqID, seq in parse_fasta(accession):

            for motif in self.motifs:
                matches = re.finditer(motif, seq)
                table = []

                for match in matches:

                    sequence = match.group()
                    gc_content = sequence.count('g') + sequence.count('c')
                    strand = "+" if sequence.startswith(self.stacker) else "-"

                    if strand == "+":
                        stacker = self.stacker
                    else:
                        stacker = self.complementary_stacker

                    stackers = re.findall("%s{%s,}" % (stacker, self.minrep), sequence)
                    total_stacker_len = sum(map(len, stackers))

                    spacers = re.sub("%s{%s,}" % (stacker, self.minrep), "|", sequence)
                    total_spacer_len = sum(map(len, spacers.split("|")))

                    start = match.start()
                    end = match.end()

                    assert seq[start: end] == sequence, f"Invalid sequence detected at chromosome {seqID} at position ({start},{end})."

                    table.append(
                        {
                            "seqID": seqID,
                            "start": start,
                            "end": end,
                            "sequence": sequence,
                            "strand": strand,
                            "gc_content": gc_content,
                            "length": len(sequence),
                            "stackerLength": total_stacker_len,
                            "spacerLength": total_spacer_len,
                        }
                    )

                table = sorted(
                               table, 
                               key=lambda x: (x['seqID'], x['start']), 
                               reverse=False
                               )
                yield table


