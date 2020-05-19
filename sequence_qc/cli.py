import click

from sequence_qc import sequence_qc


@click.command()
@click.option("--ref_fasta", help="Path to reference fasta, containing all regions in bed_file")
@click.option("--bam_file", help="Path to BAM file for calculating noise")
@click.option("--bed_file", help="Path to BED file containing regions over which to calculate noise")
@click.option("--threshold", default=0.02, help="Alt allele frequency past which to ignore positions from the calculation")
@click.option("--add_indels", default=False, help="Include bases from insertions / deletions in noise calculation")
@click.option("--truncate", default=True, help="Whether to exclude trailing bases from reads that only partially overlap the bed file")
@click.option("--ignore_overlaps", default=True, help="Don't double count overlapping reads (fragment aware)")
@click.option("--flag_filter", default=0, help="Reads with any of these flags set will be excluded from the calculation")
def calculate_noise(ref_fasta, bam_file, bed_file, threshold, add_indels, truncate, ignore_overlaps, flag_filter):
    """
    Calculate noise level of given bam file, across the given positions in `bed_file`.

    :param count:
    :param name:
    :return:
    """
    noise_level = sequence_qc.calculate_noise(
        ref_fasta=ref_fasta,
        bam_path=bam_file,
        bed_file_path=bed_file,
        noise_threshold=threshold,
        add_indels=add_indels,
        truncate=truncate,
        ignore_overlaps=ignore_overlaps,
        flag_filter=flag_filter,
    )

    # todo: add parameter -o for output file and print to there
    print(noise_level)



# for AlignmentFile?
# @click.option("--min_mapq", default=False, help="Exclude reads with a lower mapping quality")
# @click.option("--min_basq", default=False, help="Exclude bases with a lower base quality")
