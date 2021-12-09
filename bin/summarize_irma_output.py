#!/usr/bin/env python

"""
Summarize results from IRMA.
"""

import logging
import os
import re
# from collections import defaultdict
from dataclasses import dataclass, fields
from functools import partial
from itertools import repeat
# from multiprocessing import Pool
# from typing import Dict, List, Optional, Tuple

import click
import pandas as pd
# import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from rich.console import Console
from rich.logging import RichHandler

LOG_FORMAT = "%(asctime)s %(levelname)s: %(message)s [in %(filename)s:%(lineno)d]"
logging.basicConfig(format=LOG_FORMAT, level=logging.INFO)

@dataclass
class SpecimenStatsEntry:
    ID: str
    Initials: int
    Assembled: int
    WidowsPairs: int
    QCFailed: int
    Chimeric: int
    No_Match: int
    Comment: str

def parse_irma_specimen_stats(irma_result: str):
    logging.info(f"Parsing IRMA specimen statistics results for {irma_result}")
    specimen_stats_snippet = pd.DataFrame()
    try:
        read_counts = pd.read_csv(f"{irma_result}/tables/READ_COUNTS.txt", sep="\t", index_col="Record")
        specimen_stats_entry = SpecimenStatsEntry(
            ID=irma_result,
            Initials=read_counts["Reads"].get("1-initial"),
            Assembled=read_counts["Reads"].get("3-match"),
            WidowsPairs=read_counts["PairsAndWidows"].get("3-match"),
            QCFailed=read_counts["Reads"].get("2-failQC"),
            Chimeric=read_counts["Reads"].get("3-chimeric"),
            No_Match=read_counts["Reads"].get("3-nomatch"),
            Comment = ""
        )
    except FileNotFoundError:
        try:
            sorted_read_stats = pd.read_csv(f"{irma_result}/sorted_read_stats.txt", sep="\t")
            specimen_stats_entry = SpecimenStatsEntry(
                ID=irma_result,
                Initials=int(sorted_read_stats["Read Patterns"].sum()),
                Assembled=0,
                WidowsPairs=0,
                QCFailed=0,
                Chimeric=0,
                No_Match=0,
                Comment = "FluReads"
            )
        except FileNotFoundError:
            if os.path.exists(f"{irma_result}/R0.fa"):
                specimen_stats_entry = SpecimenStatsEntry(
                    ID=irma_result,
                    Initials=0,
                    Assembled=0,
                    WidowsPairs=0,
                    QCFailed=0,
                    Chimeric=0,
                    No_Match=0,
                    Comment = "NID"
                )
            else:
                specimen_stats_entry = SpecimenStatsEntry(
                    ID=irma_result,
                    Initials=0,
                    Assembled=0,
                    WidowsPairs=0,
                    QCFailed=0,
                    Chimeric=0,
                    No_Match=0,
                    Comment = "NO `tables/READ_COUNTS.txt` OR `sorted_read_stats.txt` OR `R0.fa` WAS FOUND!!!!!"
                )
    return pd.DataFrame([specimen_stats_entry]).astype({field.name: field.type for field in fields(SpecimenStatsEntry)})


@dataclass
class LdtSummaryEntry:
    ID: str
    Segment: pd.CategoricalDtype()
    Reads: int
    Align_Score: float
    Depth: int
    Sites: int
    Perc_Cover: float
    AvgDepth: float
    SDDepth: float

@dataclass
class LdtReadEntry:
    ID: str
    Segment: pd.CategoricalDtype()
    Reads: int

@dataclass
class LdtAlignEntry:
    ID: str
    Segment: pd.CategoricalDtype()
    Align_Score: float

@dataclass
class LdtCoverageEntry:
    ID: str
    Segment: pd.CategoricalDtype()
    Depth: int
    Sites: int
    Perc_Cover: float
    AvgDepth: float
    SDDepth: float

def parse_irma_ldt_fasta_result(irma_result: str, coverage_depth_threshold: int):
    logging.info(f"Parsing IRMA LDT and FASTA results for {irma_result}")
    fasta_dictionary = dict()
    # ldt_snippet = pd.DataFrame(columns=[field.name for field in fields(LdtSummaryEntry)]) # .astype({field.name: field.type for field in fields(LdtSummaryEntry)})
    ldt_snippet = pd.DataFrame()

    fasta_files = list(filter(lambda x: x.endswith(".fasta"), os.listdir(irma_result)))
    for fasta_file in fasta_files:
        for seq_record in SeqIO.parse(f"{irma_result}/{fasta_file}", "fasta"):
            fasta_dictionary[seq_record.id] = seq_record.seq

    record_name_pattern = re.compile("|".join(fasta_dictionary.keys()))

    with open(f"{irma_result}/logs/READ_log.txt") as f:
        ldt_read_snippet = pd.DataFrame(columns=[field.name for field in fields(LdtReadEntry)]) # .astype({field.name: field.type for field in fields(LdtReadEntry)})
        for line in f.readlines():
            if (record_name_match := record_name_pattern.search(line)) is not None:
                this_record = LdtReadEntry(ID=irma_result, Segment=record_name_match.group(0), Reads=int(line.split(":")[1].strip()))
                ldt_read_snippet = ldt_read_snippet.append([this_record])

        ldt_snippet = ldt_snippet.append(ldt_read_snippet)

    with open(f"{irma_result}/logs/ASSEMBLY_log.txt") as f:
        ldt_align_snippet = pd.DataFrame(columns=[field.name for field in fields(LdtAlignEntry)]) # .astype({field.name: field.type for field in fields(LdtAlignEntry)})
        for line in f.readlines():
            if (record_name_match := record_name_pattern.search(line)) is not None:
                this_record = LdtAlignEntry(ID=irma_result, Segment=record_name_match.group(0), Align_Score=float(line.split()[0].strip()))
                ldt_align_snippet = ldt_align_snippet.append([this_record])

        ldt_snippet = ldt_snippet.merge(ldt_align_snippet.groupby(["ID", "Segment"], as_index=False).agg({"Align_Score": "mean"}), how="right")

    ldt_coverage_snippet = pd.DataFrame(columns=[field.name for field in fields(LdtCoverageEntry)]) # .astype({field.name: field.type for field in fields(LdtCoverageEntry)})
    for segment_name in sorted(fasta_dictionary.keys()):
        segment_coverage = pd.read_csv(f"{irma_result}/tables/{segment_name}-coverage.txt", sep="\t")
        this_record = LdtCoverageEntry(
            ID=irma_result,
            Segment=segment_name,
            Depth=float((segment_coverage["Coverage Depth"] > coverage_depth_threshold).sum()),
            Sites=len(segment_coverage),
            Perc_Cover=float((segment_coverage["Coverage Depth"] > coverage_depth_threshold).value_counts(normalize=True)),
            AvgDepth=float(segment_coverage["Coverage Depth"].mean()),
            SDDepth=float(segment_coverage["Coverage Depth"].std(ddof=0))
        )
        ldt_coverage_snippet = ldt_coverage_snippet.append([this_record])

    ldt_snippet = ldt_snippet.merge(ldt_coverage_snippet, how="right")

    fasta_snippet = pd.DataFrame(columns=["ID", "Segment", "Sequence"])
    for segment, sequence in fasta_dictionary.items():
        fasta_snippet = fasta_snippet.append({"ID": irma_result, "Segment": segment, "Sequence": sequence}, ignore_index=True) 

    return ldt_snippet.astype({field.name: field.type for field in fields(LdtSummaryEntry)}), fasta_snippet


@click.command()
@click.option("--coverage-depth-threshold", default=25)
@click.argument("irma_results", nargs=-1)
def main(irma_results, coverage_depth_threshold):
    from rich.traceback import install

    install(show_locals=True, width=120, word_wrap=True)
    logging.basicConfig(
        format="%(message)s",
        datefmt="[%Y-%m-%d %X]",
        level=logging.DEBUG,
        handlers=[RichHandler(rich_tracebacks=True, tracebacks_show_locals=True)],
    )

    df_specimen_stats = pd.concat(map(parse_irma_specimen_stats, irma_results))
    df_specimen_stats.to_csv("Specimen_Stats.txt", sep="\t", index=False)

    df_ldt, df_fasta = map(partial(pd.concat, ignore_index=True), zip(*map(parse_irma_ldt_fasta_result, df_specimen_stats[df_specimen_stats["Comment"] == ""]["ID"], repeat(coverage_depth_threshold))))
    df_ldt.to_csv("Ldt_Summary.txt", sep="\t", index=False)

    os.mkdir("summary_fasta")
    os.chdir("summary_fasta")
    segment_names = df_fasta["Segment"].unique()
    list(map(partial(os.makedirs), map(lambda x: x.replace("_", "/", 1), segment_names)))
    # genera = {x.split("_")[0] for x in segment_names}
    # for genus in genera:
    for segment_name in segment_names:
        segment_records = []
        for _, record in df_fasta[df_fasta["Segment"] == segment_name].iterrows():
            seq_record = SeqRecord(record["Sequence"], id=record["ID"], description="")
            SeqIO.write(seq_record, f"{segment_name.replace('_', '/', 1)}/{record['ID']}.fasta", "fasta")
            segment_records.append(seq_record)
        SeqIO.write(segment_records, f"{segment_name.replace('_', '/', 1)}_all.fasta", "fasta")


if __name__ == "__main__":
    main()
