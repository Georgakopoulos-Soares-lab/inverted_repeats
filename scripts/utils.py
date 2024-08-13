from Bio import SeqIO
import time
import gzip
import json
import os
from pathlib import Path
import threading
import logging
from typing import Optional

class ProgressTracker:

    def __init__(self, total_accessions: int, 
                       filename: str = "biologs/tracker.log", 
                       sleeping_time: int = 300,
                       bucket_id: Optional[int] = None) -> None:
        self.log_filename = Path(filename).resolve()
        self.log_filename.parent.mkdir(exist_ok=True)

        self.bucket_id = bucket_id

        self.sleeping_time = sleeping_time

        self.counter = 0
        self.total_accessions = total_accessions

        logging.basicConfig(
                            filename=self.log_filename,
                            level=logging.INFO, 
                            format="%(asctime)s:%(levelname)s:%(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S",
                        )
    
    def _get_progress(self) -> str:
        return f"{self.counter * 1e2 / self.total_accessions:.2f}"

    def track_progress(self) -> None:
        bucket_info = f"Bucket: {self.bucket_id};" if self.bucket_id else ""
        while True:
            logging.info(f"{bucket_info}Progress level: {self._get_progress()}%.")
            time.sleep(self.sleeping_time)



def parse_fasta(accession: os.PathLike[str]) -> tuple[str]:
    accession = Path(accession).resolve()
    if accession.name.endswith(".gz"):
        file = gzip.open(accession, 'rt')
    else:
        file = open(accession, mode='r', encoding='UTF-8')
        
    for record in SeqIO.parse(file, 'fasta'):
        yield str(record.id), str(record.seq).lower()

    file.close()

def load_bucket(bucket_id: int, schedule_path: os.PathLike[str]) -> list[str]:
    with open(schedule_path, mode="r", encoding="UTF-8") as f:
        return json.load(f)[str(bucket_id)]
