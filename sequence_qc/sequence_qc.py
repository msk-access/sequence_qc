from pysam import AlignmentFile
from pybedtools import BedTool
from pyfaidx import Fasta


def calculate_noise(ref_fasta, bam_path, bed_file_path, noise_threshold):
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysam

    :param ref_fasta: string - Path to reference fasta which will be used to compare to
        `bam_path`.
    :param bam_path: string - Path to a BAM file which will be interrogated for noise
        level.
    :param bed_file_path: string - Path to a bed-format file with positions
        over which to calculate the noise level.
    :noise_threshold: float - Upper limit to consider a position "noisy".
        This should be a fraction between 0 and 1, which will determine whether a
        position in the supplied bed_file with any alt reads has too many to be considered
        just noise.

    :return:
    """
    # todo: why do we need to use check_seq?
    ref = Fasta(ref_fasta)
    af = AlignmentFile(bam_path, check_sq=False)
    bed_file = BedTool(bed_file_path)

    alt_count = 0
    total_count = 1e-9

    for region in bed_file.intervals:
        pileup = af.pileup(region.chrom.replace('chr', ''), region.start, region.stop)

        for p in pileup:
            refbase = ref[region.chrom][p.pos:p.pos + 1]
            # todo: need this?
            refbase = str(refbase)

            bases = p.get_query_sequences()
            mismatches = list(filter((refbase).__ne__, bases))
            alt_count += len(mismatches)
            total_count += len(bases)

    return alt_count / total_count
