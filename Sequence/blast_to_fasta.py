import json
import argparse
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO

# Function to create a FASTA record from BLAST hit data


def create_fasta_record(hit):
    hit_def = hit["hit_def"]
    hit_id = hit["hit_id"]
    hit_acc = hit["hit_acc"]
    hit_identity = hit["hit_hsps"][0]["hsp_identity"]
    hit_seq = hit["hit_hsps"][0]["hsp_hseq"]

    # Create a header for the FASTA file
    fasta_header = f"{hit_id}"

    # Create a SeqRecord for the FASTA format
    record = SeqRecord(
        Seq(hit_seq),
        id=fasta_header,
        description=f"Protein: {hit_def} | Identity: {hit_identity}%",
    )
    return record


# Main function to read JSON and write FASTA


def process_blast_json(input_path, output_path):
    # Load BLAST results from JSON file
    with open(input_path, "r") as json_file:
        blast_data = json.load(json_file)

    # List to store FASTA records
    fasta_records = []

    # Process each BLAST hit and create a FASTA record
    for hit in blast_data["hits"]:
        record = create_fasta_record(hit)
        fasta_records.append(record)

    # Output the FASTA file
    SeqIO.write(fasta_records, output_path, "fasta")

    print(f"FASTA file generated: {output_path}")


# Set up argparse for command-line arguments
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process BLAST JSON to FASTA")
    parser.add_argument(
        "-i", "--input", required=True, help="Path to input BLAST JSON file"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Path to output FASTA file"
    )

    args = parser.parse_args()

    # Run the BLAST processing
    process_blast_json(args.input, args.output)
