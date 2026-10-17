#!/usr/bin/env python3
"""
BOLD Barcode Profiler
======================
Automated pipeline for identifying and profiling host and vector DNA barcodes
(e.g., COI, COII, 16S, cytb) from metagenomic contigs/assemblies against BOLD Systems data.

Author: Alexandre F. da Silva
License: MIT
"""

import os
import sys
import csv
import argparse
import subprocess
import shutil
import re
from typing import Dict, List, Optional, Set
import pandas as pd
import numpy as np


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="BOLD Barcode Profiler: End-to-end pipeline to screen metagenomic contigs against BOLD Systems database.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-q", "--query", required=True,
        help="Path to input assembled contigs/sequences in FASTA format."
    )
    parser.add_argument(
        "-b", "--bold", required=True,
        help="Path to raw BOLD Systems TSV file (downloaded from BOLD Systems for target taxon)."
    )
    parser.add_argument(
        "-o", "--outdir", default="bold_results",
        help="Output directory for all generated databases, alignments, tables, and figures."
    )
    parser.add_argument(
        "-p", "--prefix", default=None,
        help="Prefix for output result files. (Default: query file basename without extension)."
    )
    parser.add_argument(
        "-t", "--threads", type=int, default=4,
        help="Number of CPU threads to use for BLASTn search."
    )
    parser.add_argument(
        "-m", "--marker", default="all",
        help="Target marker gene(s) to filter (e.g., 'COI-5P', 'COI', 'COII', '16S', 'cytb' or 'all'). "
             "Multiple markers can be comma-separated, e.g. 'COI-5P,COII'."
    )
    parser.add_argument(
        "--min-identity", type=float, default=80.0,
        help="Minimum percentage nucleotide identity cutoff (pident)."
    )
    parser.add_argument(
        "--min-coverage", type=float, default=50.0,
        help="Minimum query alignment coverage cutoff (qcovs, in percent)."
    )
    parser.add_argument(
        "--evalue", type=float, default=1e-5,
        help="Maximum BLASTn E-value threshold."
    )
    parser.add_argument(
        "--max-target-seqs", type=int, default=5,
        help="Maximum number of aligned sequences to keep per query in BLASTn."
    )
    parser.add_argument(
        "--skip-blast", action="store_true",
        help="Skip database generation and BLASTn execution if raw BLAST results already exist in output directory."
    )
    parser.add_argument(
        "--no-plot", action="store_true",
        help="Disable automatic figure plotting (only generate tables)."
    )
    parser.add_argument(
        "--no-extract-fastas", action="store_true",
        help="Disable automatic extraction of identified contig sequences into marker-specific FASTA files."
    )
    parser.add_argument(
        "--remote-ncbi", action="store_true",
        help="Execute remote NCBI BLASTn against GenBank (db: nt) for the extracted marker contigs to cross-validate results."
    )
    parser.add_argument(
        "--ncbi-db", default="nt",
        help="NCBI database for remote BLASTn search (default: 'nt')."
    )
    parser.add_argument(
        "--ncbi-evalue", type=float, default=1e-5,
        help="Maximum E-value threshold for remote NCBI BLASTn search."
    )
    parser.add_argument(
        "--ncbi-max-target-seqs", type=int, default=5,
        help="Maximum target sequences to report per query in remote NCBI BLASTn."
    )
    parser.add_argument(
        "--ncbi-entrez-query", default=None,
        help="Optional Entrez query filter for remote NCBI BLASTn (e.g. 'Culicidae[Organism]' or 'Diptera[Organism]')."
    )
    return parser.parse_args()


def check_tool_dependency(tool_name: str) -> str:
    """Verifies that an external binary (e.g. makeblastdb, blastn) exists in PATH."""
    path = shutil.which(tool_name)
    if not path:
        sys.exit(
            f"[-] Error: '{tool_name}' was not found in PATH.\n"
            f"    Please activate the required Conda/Micromamba environment (e.g., 'conda activate bold-profiler')\n"
            f"    or install NCBI BLAST+ (e.g., 'conda install -c bioconda blast')."
        )
    return path


def validate_inputs(query_path: str, bold_path: str):
    """Checks that the input files exist and are not empty."""
    if not os.path.isfile(query_path):
        sys.exit(f"[-] Error: Query FASTA file not found: {query_path}")
    if os.path.getsize(query_path) == 0:
        sys.exit(f"[-] Error: Query FASTA file is empty: {query_path}")

    if not os.path.isfile(bold_path):
        sys.exit(f"[-] Error: BOLD TSV file not found: {bold_path}")
    if os.path.getsize(bold_path) == 0:
        sys.exit(f"[-] Error: BOLD TSV file is empty: {bold_path}")


def prepare_bold_fasta(bold_tsv: str, out_fasta: str) -> int:
    """
    Parses the BOLD TSV file, extracts valid nucleotide sequences,
    cleans up gaps/whitespace, and outputs a formatted FASTA file.
    Header is kept strictly as >processid to avoid BLAST+ 50-character limit errors.
    """
    print(f"[*] Parsing BOLD TSV and extracting FASTA: {bold_tsv}")
    seq_count = 0
    with open(bold_tsv, "r", encoding="utf-8", errors="ignore") as in_f, \
         open(out_fasta, "w", encoding="utf-8") as out_f:
        
        reader = csv.reader(in_f, delimiter="\t")
        header = next(reader, None)
        if not header:
            sys.exit(f"[-] Error: BOLD TSV appears to have no header line: {bold_tsv}")

        col_map = {col.strip().lower(): idx for idx, col in enumerate(header)}
        if "processid" not in col_map or "nuc" not in col_map:
            sys.exit("[-] Error: BOLD TSV must contain 'processid' and 'nuc' columns.")

        pid_idx = col_map["processid"]
        nuc_idx = col_map["nuc"]

        for row in reader:
            if len(row) <= max(pid_idx, nuc_idx):
                continue
            pid = row[pid_idx].strip()
            raw_seq = row[nuc_idx].strip()

            if not pid or not raw_seq:
                continue

            # Clean nucleotide sequence
            clean_seq = "".join([c for c in raw_seq.upper() if c in "ACGTUNRYKMSWBDHV"])
            if len(clean_seq) < 50:  # Minimum valid barcode sequence length
                continue

            out_f.write(f">{pid}\n{clean_seq}\n")
            seq_count += 1

    print(f"[+] Successfully extracted {seq_count:,} sequences into: {out_fasta}")
    return seq_count


def make_blast_database(fasta_file: str, db_prefix: str):
    """Builds a nucleotide BLAST database using makeblastdb."""
    check_tool_dependency("makeblastdb")
    print(f"[*] Building BLAST database '{db_prefix}' from: {fasta_file}")
    cmd = [
        "makeblastdb",
        "-in", fasta_file,
        "-dbtype", "nucl",
        "-out", db_prefix,
        "-title", "BOLD_DB"
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.exit(f"[-] Error building BLAST database:\n{res.stderr}")
    print("[+] BLAST database built successfully.")


def run_blastn(query_fasta: str, db_prefix: str, output_tsv: str, threads: int,
               evalue: float, max_target_seqs: int):
    """Executes BLASTn search with tabular format including query coverage."""
    check_tool_dependency("blastn")
    print(f"[*] Running BLASTn search: {query_fasta} vs {db_prefix} (threads: {threads})")
    
    outfmt = "6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"
    cmd = [
        "blastn",
        "-query", query_fasta,
        "-db", db_prefix,
        "-outfmt", outfmt,
        "-max_target_seqs", str(max_target_seqs),
        "-evalue", str(evalue),
        "-num_threads", str(threads),
        "-out", output_tsv
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if res.returncode != 0:
        sys.exit(f"[-] Error running blastn:\n{res.stderr}")
    print(f"[+] BLASTn search completed: {output_tsv}")


def load_bold_metadata(bold_tsv: str, target_process_ids: Set[str]) -> Dict[str, Dict[str, str]]:
    """
    Extracts metadata from the BOLD TSV file only for the processids that had BLAST hits.
    This saves significant memory when handling multi-gigabyte BOLD dumps.
    """
    print(f"[*] Loading BOLD metadata for {len(target_process_ids):,} matched IDs...")
    metadata = {}
    with open(bold_tsv, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader, None)
        if not header:
            return metadata

        col_map = {col.strip().lower(): idx for idx, col in enumerate(header)}
        pid_idx = col_map.get("processid")
        if pid_idx is None:
            return metadata

        # Select relevant metadata columns
        fields_to_extract = {
            "sampleid": col_map.get("sampleid"),
            "marker_code": col_map.get("marker_code"),
            "class": col_map.get("class"),
            "order": col_map.get("order"),
            "family": col_map.get("family"),
            "subfamily": col_map.get("subfamily"),
            "genus": col_map.get("genus"),
            "species": col_map.get("species"),
            "subspecies": col_map.get("subspecies"),
            "identification": col_map.get("identification"),
            "country": col_map.get("country"),
            "province_state": col_map.get("province_state"),
        }

        for row in reader:
            if len(row) <= pid_idx:
                continue
            pid = row[pid_idx].strip()
            if pid in target_process_ids:
                row_meta = {}
                for field, idx in fields_to_extract.items():
                    val = row[idx].strip() if (idx is not None and len(row) > idx) else ""
                    row_meta[field] = val

                # Fallback for species if empty: check identification or genus
                if not row_meta.get("species"):
                    ident_val = row[col_map["identification"]].strip() if ("identification" in col_map and len(row) > col_map["identification"]) else ""
                    if ident_val:
                        row_meta["species"] = ident_val
                    elif row_meta.get("genus"):
                        row_meta["species"] = f"{row_meta['genus']} sp."

                metadata[pid] = row_meta

    return metadata


def annotate_and_filter_hits(blast_tsv: str, bold_tsv: str, out_prefix: str,
                             min_identity: float, min_coverage: float,
                             target_markers: Optional[List[str]]) -> Optional[pd.DataFrame]:
    """
    Annotates BLAST hits with BOLD metadata, applies user thresholds,
    determines best hits per query contig, and exports tabular summaries.
    """
    cols = [
        "qseqid", "qlen", "sseqid", "slen", "pident", "length",
        "mismatch", "gapopen", "qstart", "qend", "sstart", "send",
        "evalue", "bitscore", "qcovs", "qcovhsp"
    ]

    if not os.path.exists(blast_tsv) or os.path.getsize(blast_tsv) == 0:
        print("[-] Warning: BLAST output is empty. No hits were found.")
        return None

    df_blast = pd.read_csv(blast_tsv, sep="\t", names=cols)
    if df_blast.empty:
        print("[-] Warning: BLAST output table contains no records.")
        return None

    # Clean processid in case sseqid contains pipes (|) from legacy databases
    df_blast["clean_processid"] = df_blast["sseqid"].astype(str).apply(lambda x: x.split("|")[0].strip())
    hit_ids = set(df_blast["clean_processid"].unique())
    meta_dict = load_bold_metadata(bold_tsv, hit_ids)

    # Annotate dataframe with metadata from BOLD TSV (with fallback to sseqid tokens)
    meta_rows = []
    for _, row in df_blast.iterrows():
        pid = row["clean_processid"]
        meta = meta_dict.get(pid, {}).copy()

        # Fallback to tokens inside pipe-separated sseqid if any metadata is empty
        sseq_parts = [p.strip() for p in str(row["sseqid"]).split("|")]
        if not meta.get("marker_code") and len(sseq_parts) > 2 and sseq_parts[2] not in ["", "NA"]:
            meta["marker_code"] = sseq_parts[2]
        if not meta.get("species") and len(sseq_parts) > 7 and sseq_parts[7] not in ["", "NA"]:
            meta["species"] = sseq_parts[7].replace("_", " ")
        if not meta.get("genus") and len(sseq_parts) > 6 and sseq_parts[6] not in ["", "NA"]:
            meta["genus"] = sseq_parts[6]
        if not meta.get("family") and len(sseq_parts) > 5 and sseq_parts[5] not in ["", "NA"]:
            meta["family"] = sseq_parts[5]

        meta_rows.append(meta)

    df_meta = pd.DataFrame(meta_rows)
    df_annotated = pd.concat([df_blast.reset_index(drop=True), df_meta.reset_index(drop=True)], axis=1)
    df_annotated["processid"] = df_annotated["clean_processid"]

    # Apply quality filters
    initial_count = len(df_annotated)
    df_filtered = df_annotated[
        (df_annotated["pident"] >= min_identity) &
        (df_annotated["qcovs"] >= min_coverage)
    ].copy()

    # Apply marker filtering if requested
    if target_markers and "all" not in [m.lower() for m in target_markers]:
        target_norm = [m.lower().strip() for m in target_markers]
        df_filtered = df_filtered[
            df_filtered["marker_code"].astype(str).str.lower().isin(target_norm)
        ]

    print(f"[*] Hits before filtering: {initial_count:,} | Passing thresholds: {len(df_filtered):,}")

    if df_filtered.empty:
        print("[-] Warning: No hits passed the specified quality/marker filters.")
        return None

    # Save all passing hits
    all_hits_file = f"{out_prefix}_all_filtered_hits.tsv"
    df_filtered.to_csv(all_hits_file, sep="\t", index=False)
    print(f"[+] All filtered hits saved: {all_hits_file}")

    # Best hit selection: sort by bitscore desc, pident desc, length desc, then drop duplicates per qseqid
    df_besthits = df_filtered.sort_values(
        by=["bitscore", "pident", "length"],
        ascending=[False, False, False]
    ).drop_duplicates(subset=["qseqid"]).reset_index(drop=True)

    # Save best hits (TSV and CSV)
    besthits_tsv = f"{out_prefix}_besthits.tsv"
    besthits_csv = f"{out_prefix}_besthits.csv"
    df_besthits.to_csv(besthits_tsv, sep="\t", index=False)
    df_besthits.to_csv(besthits_csv, index=False)
    print(f"[+] Best hits table saved: {besthits_tsv}")
    print(f"[+] Best hits CSV saved: {besthits_csv}")

    # Generate Google Docs friendly 2-column format
    two_col_file = f"{out_prefix}_summary_2columns.tsv"
    two_col_rows = []
    for idx, r in df_besthits.iterrows():
        two_col_rows.append([f"=== CONTIG {idx+1}: {r['qseqid']} ===", ""])
        two_col_rows.append(["Contig Length (bp)", str(r["qlen"])])
        two_col_rows.append(["Marker Gene", str(r.get("marker_code", "NA"))])
        two_col_rows.append(["Assigned Species", str(r.get("species", "NA"))])
        two_col_rows.append(["Assigned Genus", str(r.get("genus", "NA"))])
        two_col_rows.append(["Family", str(r.get("family", "NA"))])
        two_col_rows.append(["Order / Class", f"{r.get('order', 'NA')} / {r.get('class', 'NA')}"])
        two_col_rows.append(["Country / Locality", f"{r.get('country', '')} {r.get('province_state', '')}".strip() or "NA"])
        two_col_rows.append(["Percent Identity (% pident)", f"{r['pident']:.2f}%"])
        two_col_rows.append(["Query Coverage (% qcovs)", f"{r['qcovs']:.1f}%"])
        two_col_rows.append(["Aligned Length (bp)", str(r["length"])])
        two_col_rows.append(["E-value", str(r["evalue"])])
        two_col_rows.append(["Bit Score", str(r["bitscore"])])
        two_col_rows.append(["BOLD Process ID", str(r.get("processid", "NA"))])
        two_col_rows.append(["BOLD Sample ID", str(r.get("sampleid", "NA"))])
        two_col_rows.append(["", ""])

    with open(two_col_file, "w", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Parameter / Metadata", "Value / Description"])
        writer.writerows(two_col_rows)

    print(f"[+] Two-column summary saved (Google Docs friendly): {two_col_file}")
    return df_besthits


def extract_marker_fastas(query_fasta: str, df_besthits: pd.DataFrame, out_prefix: str) -> Dict[str, str]:
    """
    Extracts contig sequences identified in the best hits table from the query FASTA,
    grouped and saved into separate FASTA files by marker gene (e.g. COI-5P, 18S-3P),
    as well as a consolidated FASTA file with all identified marker contigs.
    """
    print(f"\n[*] Extracting identified marker contigs from: {query_fasta}")

    # Read query FASTA sequences into memory
    seq_records: Dict[str, str] = {}
    current_id = None
    current_seq = []
    with open(query_fasta, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith(">"):
                if current_id:
                    seq_records[current_id] = "".join(current_seq)
                current_id = line_str[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line_str)
        if current_id:
            seq_records[current_id] = "".join(current_seq)

    # Group hits by sanitized marker name
    markers_dict: Dict[str, List[Dict[str, str]]] = {}
    all_extracted_records = []

    for _, row in df_besthits.iterrows():
        qid = str(row["qseqid"]).strip()
        seq = seq_records.get(qid)
        if not seq:
            # Fallback search if query ID has extra tokens or differences
            for k, v in seq_records.items():
                if k == qid or k.startswith(qid) or qid.startswith(k):
                    seq = v
                    break

        if not seq:
            print(f"[-] Warning: Contig sequence '{qid}' not found in query FASTA.")
            continue

        marker = str(row.get("marker_code", "")).strip()
        if not marker or marker.lower() in ["na", "nan", "none", ""]:
            marker = "Unassigned"

        clean_marker = "".join([c if c.isalnum() or c in "-_." else "_" for c in marker])
        species_name = str(row.get("species", "")).strip() or f"{row.get('genus', 'Taxon')} sp."
        pident = float(row.get("pident", 0.0))
        qcovs = float(row.get("qcovs", 0.0))
        bitscore = int(row.get("bitscore", 0))
        qlen = int(row.get("qlen", len(seq)))
        processid = str(row.get("processid", "NA")).strip()

        # FASTA header: first token is clean contig ID, followed by descriptive annotations
        annotated_header = (
            f">{qid} marker={marker} species=\"{species_name}\" "
            f"identity={pident:.2f}% coverage={qcovs:.1f}% length={qlen}bp "
            f"bitscore={bitscore} bold_id={processid}"
        )

        rec = {
            "qseqid": qid,
            "header": annotated_header,
            "seq": seq,
            "marker": marker,
            "clean_marker": clean_marker,
            "species": species_name,
            "pident": pident,
            "qcovs": qcovs,
            "bitscore": bitscore,
            "processid": processid
        }
        all_extracted_records.append(rec)
        markers_dict.setdefault(clean_marker, []).append(rec)

    output_fastas: Dict[str, str] = {}

    # Write separate FASTA per marker
    for clean_marker, records in markers_dict.items():
        marker_fasta = f"{out_prefix}_contigs_{clean_marker}.fasta"
        with open(marker_fasta, "w", encoding="utf-8") as out_f:
            for r in records:
                out_f.write(f"{r['header']}\n")
                s = r["seq"]
                for i in range(0, len(s), 80):
                    out_f.write(s[i:i+80] + "\n")
        output_fastas[clean_marker] = marker_fasta
        print(f"[+] Marker FASTA saved ({len(records)} sequences): {marker_fasta}")

    # Write combined FASTA containing all identified marker contigs
    if all_extracted_records:
        all_fasta = f"{out_prefix}_contigs_all_markers.fasta"
        with open(all_fasta, "w", encoding="utf-8") as out_f:
            for r in all_extracted_records:
                out_f.write(f"{r['header']}\n")
                s = r["seq"]
                for i in range(0, len(s), 80):
                    out_f.write(s[i:i+80] + "\n")
        output_fastas["ALL"] = all_fasta
        print(f"[+] Consolidated marker FASTA saved ({len(all_extracted_records)} sequences): {all_fasta}")

    return output_fastas


def run_remote_ncbi_blast(marker_fastas: Dict[str, str], df_besthits: pd.DataFrame,
                          out_prefix: str, ncbi_db: str = "nt", evalue: float = 1e-5,
                          max_target_seqs: int = 5, entrez_query: Optional[str] = None):
    """
    Executes remote NCBI BLASTn against GenBank (e.g. 'nt') for the extracted marker contigs
    using NCBI's fast interactive QBLAST API and produces a cross-validation comparison table
    between BOLD and NCBI GenBank identifications.
    """
    import urllib.request
    import urllib.parse
    import xml.etree.ElementTree as ET
    import time
    import re

    print(f"\n[*] Initiating remote NCBI BLASTn search against '{ncbi_db}' database...")
    print(f"[*] Note: Queries are sent to NCBI servers via HTTPS QBLAST API.")
    if entrez_query:
        print(f"[*] Entrez Query Filter: {entrez_query}")

    target_fasta = marker_fastas.get("ALL")
    if not target_fasta or not os.path.exists(target_fasta):
        non_all = [p for k, p in marker_fastas.items() if k != "ALL"]
        if not non_all:
            print("[-] No marker FASTA available for remote BLASTn.")
            return
        target_fasta = non_all[0]

    # Read queries to submit
    queries: Dict[str, str] = {}
    cur_id = None
    cur_seq = []
    with open(target_fasta, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            l = line.strip()
            if not l:
                continue
            if l.startswith(">"):
                if cur_id:
                    queries[cur_id] = "".join(cur_seq)
                cur_id = l[1:].split()[0]
                cur_seq = []
            else:
                cur_seq.append(l)
        if cur_id:
            queries[cur_id] = "".join(cur_seq)

    all_ncbi_hits = []
    for qid, seq in queries.items():
        print(f"\n[*] Submitting query '{qid}' ({len(seq)} bp) to NCBI...")
        post_data = {
            'CMD': 'Put',
            'PROGRAM': 'blastn',
            'DATABASE': ncbi_db,
            'QUERY': seq,
            'EXPECT': str(evalue),
            'HITLIST_SIZE': str(max_target_seqs)
        }
        if entrez_query:
            post_data['ENTREZ_QUERY'] = entrez_query

        data = urllib.parse.urlencode(post_data).encode('utf-8')
        http_headers = {'User-Agent': 'BOLDBarcodeProfiler/1.0 (bioinformatics tool; NCBI QBLAST client)'}
        req = urllib.request.Request('https://blast.ncbi.nlm.nih.gov/Blast.cgi', data=data, headers=http_headers)

        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                html = resp.read().decode('utf-8')
        except Exception as e:
            print(f"[-] Network error submitting query to NCBI: {e}")
            continue

        rid_match = re.search(r'RID = (\w+)', html)
        rtoe_match = re.search(r'RTOE = (\d+)', html)
        if not rid_match:
            print(f"[-] Could not obtain RID from NCBI for '{qid}'.")
            continue

        rid = rid_match.group(1)
        rtoe = int(rtoe_match.group(1)) if rtoe_match else 30
        print(f"[*] Request ID (RID): {rid} | Estimated wait: {rtoe}s...")

        time.sleep(rtoe)
        poll_url = f'https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID={rid}'
        poll_req = urllib.request.Request(poll_url, headers=http_headers)

        q_hits = []
        for attempt in range(15):
            try:
                with urllib.request.urlopen(poll_req, timeout=30) as resp:
                    txt = resp.read().decode('utf-8')
            except Exception:
                time.sleep(8)
                continue

            if 'Status=READY' in txt:
                print(f"[+] NCBI search completed for '{qid}'!")
                for line in txt.splitlines():
                    if not line.startswith('#') and line.strip():
                        parts = line.split('\t')
                        if len(parts) >= 12:
                            q_hits.append({
                                'qseqid': qid,
                                'sseqid': parts[1],
                                'pident': float(parts[2]),
                                'length': int(parts[3]),
                                'mismatch': int(parts[4]),
                                'gapopen': int(parts[5]),
                                'qstart': int(parts[6]),
                                'qend': int(parts[7]),
                                'sstart': int(parts[8]),
                                'send': int(parts[9]),
                                'evalue': parts[10],
                                'bitscore': float(parts[11])
                            })
                            if len(q_hits) >= max_target_seqs:
                                break
                break
            elif 'Status=WAITING' in txt:
                print(f"[*] Query '{qid}' in NCBI queue... checking in 10s (attempt {attempt+1})")
                time.sleep(10)
            else:
                break

        # Query E-utilities for sequence titles
        for h in q_hits:
            acc = h['sseqid']
            try:
                esum_url = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={acc}&retmode=xml'
                esum_req = urllib.request.Request(esum_url, headers=http_headers)
                with urllib.request.urlopen(esum_req, timeout=15) as r:
                    root = ET.fromstring(r.read())
                    elem = root.find('.//Item[@Name="Title"]')
                    h['stitle'] = elem.text if elem is not None else ""
            except Exception:
                h['stitle'] = ""
            all_ncbi_hits.append(h)

    remote_out_tsv = f"{out_prefix}_ncbi_blastn_raw.tsv"
    if all_ncbi_hits:
        df_ncbi = pd.DataFrame(all_ncbi_hits)
        df_ncbi.to_csv(remote_out_tsv, sep="\t", index=False)
        print(f"\n[+] Raw NCBI BLASTn results saved: {remote_out_tsv}")

        # Top hit per contig for NCBI
        df_ncbi_besthits = df_ncbi.sort_values(
            by=["bitscore", "pident", "length"],
            ascending=[False, False, False]
        ).drop_duplicates(subset=["qseqid"]).reset_index(drop=True)

        # Merge with BOLD best hits
        comparison_rows = []
        for _, b_row in df_besthits.iterrows():
            qid = str(b_row["qseqid"]).strip()
            n_match = df_ncbi_besthits[df_ncbi_besthits["qseqid"].astype(str) == qid]

            comp = {
                "Contig_ID": qid,
                "Contig_Length_bp": b_row.get("qlen", "NA"),
                "Marker_Gene": b_row.get("marker_code", "NA"),
                "BOLD_Species": b_row.get("species", "NA"),
                "BOLD_Identity_pct": b_row.get("pident", "NA"),
                "BOLD_Coverage_pct": b_row.get("qcovs", "NA"),
                "BOLD_Process_ID": b_row.get("processid", "NA"),
                "NCBI_Accession": n_match["sseqid"].values[0] if not n_match.empty else "No hit",
                "NCBI_Identity_pct": n_match["pident"].values[0] if not n_match.empty else "NA",
                "NCBI_BitScore": n_match["bitscore"].values[0] if not n_match.empty else "NA",
                "NCBI_Evalue": n_match["evalue"].values[0] if not n_match.empty else "NA",
                "NCBI_Title": n_match["stitle"].values[0] if not n_match.empty else "NA"
            }
            comparison_rows.append(comp)

        df_comp = pd.DataFrame(comparison_rows)
        comp_tsv = f"{out_prefix}_bold_vs_ncbi_comparison.tsv"
        comp_csv = f"{out_prefix}_bold_vs_ncbi_comparison.csv"
        df_comp.to_csv(comp_tsv, sep="\t", index=False)
        df_comp.to_csv(comp_csv, index=False)
        print(f"[+] BOLD vs NCBI cross-validation table saved: {comp_tsv}")
        print(f"[+] BOLD vs NCBI cross-validation CSV saved: {comp_csv}")
    else:
        print("[-] No hits returned from NCBI BLASTn search.")


def generate_publication_figures(df_besthits: pd.DataFrame, out_prefix: str):
    """
    Generates a publication-grade multi-panel figure summarizing:
      - Panel A: Contig Total Length vs Alignment Length
      - Panel B: Percentage Identity & BLASTn Bit Score
      - Panel C: Structured Summary Table with key taxonomic and BOLD metadata
    Exports both 300 DPI PNG and vector SVG.
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec
    except ImportError:
        print("[-] Matplotlib not available. Skipping figure generation.")
        return

    print("[*] Generating publication-ready summary figures (PNG @ 300 DPI & SVG)...")

    df = df_besthits.copy()
    if len(df) == 0:
        return

    # Limit to top 15 contigs if there are many, sorted by length
    if len(df) > 15:
        print(f"[*] Note: Plotting top 15 contigs by length out of {len(df):,} total hits.")
        df = df.sort_values(by="qlen", ascending=False).head(15)

    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig = plt.figure(figsize=(13.5, 9.2), dpi=300)
    gs = gridspec.GridSpec(2, 2, height_ratios=[1.15, 1.1], hspace=0.38, wspace=0.28)

    x = np.arange(len(df))
    width = 0.35

    # Panel A: Contig Total Length & BOLD Alignment Length
    ax1 = fig.add_subplot(gs[0, 0])
    rects1 = ax1.bar(x - width/2, df["qlen"], width, label="Contig Total Length (bp)",
                     color="#1d3557", alpha=0.9, edgecolor="black", linewidth=0.8)
    rects2 = ax1.bar(x + width/2, df["length"], width, label="BOLD Alignment Length (bp)",
                     color="#457b9d", alpha=0.9, edgecolor="black", linewidth=0.8)

    ax1.set_ylabel("Length (bp)", fontsize=11, fontweight="bold")
    ax1.set_title("A. Contig vs. BOLD Alignment Length", fontsize=12, fontweight="bold", loc="left", pad=10)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{row['qseqid']}\n({row.get('marker_code', 'Marker')})" for _, row in df.iterrows()],
                        fontsize=9.5, fontweight="normal")
    ax1.legend(frameon=True, facecolor="white", framealpha=0.9, fontsize=9)
    ax1.grid(axis="y", linestyle="--", alpha=0.5)

    for rect in rects1:
        h = rect.get_height()
        ax1.annotate(f"{int(h)} bp", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8, fontweight="bold")
    for rect in rects2:
        h = rect.get_height()
        ax1.annotate(f"{int(h)} bp", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8)

    # Panel B: Identity % and Bit Score
    ax2 = fig.add_subplot(gs[0, 1])
    ax2_twin = ax2.twinx()

    rects_id = ax2.bar(x - width/2, df["pident"], width, label="% Identity",
                       color="#2a9d8f", alpha=0.9, edgecolor="black", linewidth=0.8)
    rects_bit = ax2_twin.bar(x + width/2, df["bitscore"], width, label="Bit Score",
                             color="#e76f51", alpha=0.9, edgecolor="black", linewidth=0.8)

    ax2.set_ylabel("Identity (%)", fontsize=11, fontweight="bold", color="#1b6d63")
    ax2_twin.set_ylabel("BLASTn Bit Score", fontsize=11, fontweight="bold", color="#b24329")
    ax2.set_title("B. Taxonomic Alignment Metrics (BOLD)", fontsize=12, fontweight="bold", loc="left", pad=10)
    ax2.set_xticks(x)
    
    # Species label
    labels = []
    for _, row in df.iterrows():
        sp = row.get("species", "")
        if pd.isna(sp) or not str(sp).strip():
            sp = f"{row.get('genus', 'Taxon')} sp."
        labels.append(f"{row['qseqid']}\n({sp})")
        
    ax2.set_xticklabels(labels, fontsize=9.5, fontweight="normal")
    min_id = max(0, df["pident"].min() - 5)
    ax2.set_ylim(min(min_id, 85), 102)
    ax2.grid(axis="y", linestyle="--", alpha=0.5)
    ax2_twin.grid(False)

    for rect in rects_id:
        h = rect.get_height()
        ax2.annotate(f"{h:.1f}%", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                     textcoords="offset points", ha="center", va="bottom", fontsize=8, fontweight="bold", color="#14524a")
    for rect in rects_bit:
        h = rect.get_height()
        ax2_twin.annotate(f"{int(h)}", xy=(rect.get_x() + rect.get_width()/2, h), xytext=(0, 3),
                          textcoords="offset points", ha="center", va="bottom", fontsize=8, fontweight="bold", color="#802611")

    # Panel C: Summary Table Card
    ax3 = fig.add_subplot(gs[1, :])
    ax3.axis("off")
    ax3.set_title("C. Summary of BOLD Marker Barcode Identifications", fontsize=12, fontweight="bold", loc="left", pad=6)

    headers = [
        "Contig\nID",
        "Contig\nLength",
        "Marker\nGene",
        "Assigned Taxon\n(BOLD Systems)",
        "Nucleotide\nIdentity",
        "Aligned\nLength",
        "Query\nCoverage",
        "Bit\nScore",
        "BOLD\nProcess ID"
    ]
    col_widths = [0.11, 0.09, 0.09, 0.23, 0.09, 0.09, 0.09, 0.08, 0.13]

    table_data = []
    for _, r in df.iterrows():
        tax_name = r.get("species", "")
        if pd.isna(tax_name) or not str(tax_name).strip():
            tax_name = f"{r.get('genus', 'Unclassified')} sp."
        country_str = r.get("country", "")
        if pd.notna(country_str) and str(country_str).strip():
            taxon_cell = f"{tax_name}\n({str(country_str).strip()})"
        else:
            taxon_cell = tax_name

        table_data.append([
            r["qseqid"],
            f"{int(r['qlen']):,} bp",
            str(r.get("marker_code", "NA")),
            taxon_cell,
            f"{float(r['pident']):.2f}%",
            f"{int(r['length']):,} bp",
            f"{float(r['qcovs']):.1f}%",
            f"{int(r['bitscore']):,}",
            str(r.get("processid", "NA"))
        ])

    table = ax3.table(
        cellText=table_data,
        colLabels=headers,
        colWidths=col_widths,
        loc="center",
        cellLoc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9.5)
    table.scale(1.0, 1.85)

    # Styling table
    for (i, j), cell in table.get_celld().items():
        if i == 0:
            cell.set_facecolor("#1d3557")
            cell.set_text_props(color="white", fontweight="bold", va="center")
        else:
            cell.set_text_props(va="center")
            if i % 2 == 0:
                cell.set_facecolor("#f8f9fa")
            else:
                cell.set_facecolor("#ffffff")
        cell.set_edgecolor("#dee2e6")

    # Remove incompatible tight_layout when using manual GridSpec margins
    out_png = f"{out_prefix}_summary_figure.png"
    out_svg = f"{out_prefix}_summary_figure.svg"
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)

    print(f"[+] Publication figure (300 DPI PNG) saved: {out_png}")
    print(f"[+] Vector figure (SVG) saved: {out_svg}")


def main():
    args = parse_arguments()

    # Determine prefix
    if not args.prefix:
        query_base = os.path.basename(args.query)
        prefix = os.path.splitext(query_base)[0]
        if prefix.endswith(".fa") or prefix.endswith(".fasta"):
            prefix = os.path.splitext(prefix)[0]
    else:
        prefix = args.prefix

    os.makedirs(args.outdir, exist_ok=True)
    validate_inputs(args.query, args.bold)

    db_fasta = os.path.join(args.outdir, "bold_db.fasta")
    db_prefix = os.path.join(args.outdir, "bold_db")
    blast_out = os.path.join(args.outdir, f"{prefix}_blastn_raw.tsv")
    final_prefix = os.path.join(args.outdir, prefix)

    target_markers = [m.strip() for m in args.marker.split(",")] if args.marker else ["all"]

    print("=" * 60)
    print(" BOLD Barcode Profiler Pipeline")
    print(f" - Query FASTA:    {args.query}")
    print(f" - BOLD TSV:       {args.bold}")
    print(f" - Output Dir:     {args.outdir}")
    print(f" - Prefix:         {prefix}")
    print(f" - Target Markers: {', '.join(target_markers)}")
    print(f" - Min Identity:   {args.min_identity}%")
    print(f" - Min Coverage:   {args.min_coverage}%")
    print("=" * 60)

    # Step 1 & 2: Build BOLD Database (unless skipping)
    if not args.skip_blast or not os.path.exists(blast_out):
        prepare_bold_fasta(args.bold, db_fasta)
        make_blast_database(db_fasta, db_prefix)
        # Step 3: Run BLASTn
        run_blastn(
            query_fasta=args.query,
            db_prefix=db_prefix,
            output_tsv=blast_out,
            threads=args.threads,
            evalue=args.evalue,
            max_target_seqs=args.max_target_seqs
        )
    else:
        print(f"[*] Reusing existing BLASTn results: {blast_out}")

    # Step 4: Annotate and Filter Hits
    df_besthits = annotate_and_filter_hits(
        blast_tsv=blast_out,
        bold_tsv=args.bold,
        out_prefix=final_prefix,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        target_markers=target_markers
    )

    # Step 5: Extract Marker-Specific Contig FASTAs
    marker_fastas = {}
    if df_besthits is not None and not args.no_extract_fastas:
        marker_fastas = extract_marker_fastas(
            query_fasta=args.query,
            df_besthits=df_besthits,
            out_prefix=final_prefix
        )

    # Step 6: Publication Summary Figures
    if df_besthits is not None and not args.no_plot:
        generate_publication_figures(df_besthits, final_prefix)

    # Step 7: Remote NCBI BLASTn Cross-Validation (Optional)
    if df_besthits is not None and args.remote_ncbi:
        run_remote_ncbi_blast(
            marker_fastas=marker_fastas,
            df_besthits=df_besthits,
            out_prefix=final_prefix,
            ncbi_db=args.ncbi_db,
            evalue=args.ncbi_evalue,
            max_target_seqs=args.ncbi_max_target_seqs,
            entrez_query=args.ncbi_entrez_query
        )

    print("\n[✓] Analysis completed successfully!")


if __name__ == "__main__":
    main()
