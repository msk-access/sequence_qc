"""Main module."""


from pysam import AlignmentFile
from pybedtools import BedTool


def calculate_noise(bam_path, bed_file_path, noise_threshold):
    """
    Create file of noise across specified regions in `bed_file` using pybedtools and pysam

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
    af = AlignmentFile(bam_path, check_sq=False)
    bed_file = BedTool(bed_file_path)

    alt_count = 0
    total_count = 0

    for region in bed_file.intervals:
        pileup = af.pileup(region.chrom.replace('chr', ''), region.start, region.stop)

        # todo: get this part working
        for position in pileup:
            alt_count += position.A if position.A < noise_threshold else 0
            alt_count += position.C if position.C < noise_threshold else 0
            alt_count += position.G if position.G < noise_threshold else 0
            alt_count += position.T if position.T < noise_threshold else 0
            total_count += position.depth

    return alt_count / total_count
