import pysam
import argparse
import math

def split_and_write_bam(input_bam, output_bam, chunk_size):
    in_bam = pysam.AlignmentFile(input_bam, "rb")
    out_bam = pysam.AlignmentFile(output_bam, "wb", template=in_bam)

    for read in in_bam.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        seq = read.query_sequence
        qual = read.query_qualities

        if not seq or not qual:
            continue

        read_len = len(seq)

        # Skip reads shorter than chunk_size
        if read_len <= chunk_size:
            out_bam.write(read)  # Optionally write as-is
            continue

        # Calculate number of chunks and adjusted chunk length
        num_chunks = math.ceil(read_len / chunk_size)
        true_chunk_size = math.ceil(read_len / num_chunks)

        for i in range(num_chunks):
            start = i * true_chunk_size
            end = min(start + true_chunk_size, read_len)

            chunk_seq = seq[start:end]
            chunk_qual = qual[start:end]

            chunk_read = pysam.AlignedSegment()
            chunk_read.query_name = f"{read.query_name}_chunk{i}"
            chunk_read.query_sequence = chunk_seq
            chunk_read.flag = read.flag
            chunk_read.reference_id = read.reference_id
            chunk_read.reference_start = read.reference_start + start  # Approximate position
            chunk_read.mapping_quality = read.mapping_quality
            chunk_read.cigar = [(0, len(chunk_seq))]  # 'M' for all bases
            chunk_read.next_reference_id = -1
            chunk_read.next_reference_start = -1
            chunk_read.template_length = 0
            chunk_read.query_qualities = chunk_qual

            out_bam.write(chunk_read)

    in_bam.close()
    out_bam.close()
    print(f"✅ Done! Output written to: {output_bam}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split long Nanopore BAM reads into ~400bp synthetic short reads (BAM output).")
    parser.add_argument("input_bam", help="Input aligned BAM file")
    parser.add_argument("output_bam", help="Output BAM file with synthetic short reads")
    parser.add_argument("--chunk", type=int, default=400, help="Target chunk size (default: 400 bp)")
    args = parser.parse_args()

    split_and_write_bam(args.input_bam, args.output_bam, args.chunk)
