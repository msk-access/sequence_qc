import pysam
import pandas as pd

from collections import defaultdict


def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.

    Reads are added to read_dict until a pair is found.

    :param: bam pysam.AlignmentFile
    :param: region str - region to search, with format chr:start-end
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(region=region_string):
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            try:
                if read.is_read1 and read.tlen == -read_dict[qname][1].tlen:
                    yield read, read_dict[qname][1]
                else:
                    if read.tlen == -read_dict[qname][0].tlen:
                        yield read_dict[qname][0], read
            except AttributeError:
                pass
            del read_dict[qname]


def get_fragment_size_for_sample(sample_id, bam_file_path, tag, noise_df, mifs, mafs):
    """
    Search through positions in `noise_df` for reads from `bam_file_path`
    and write tlen information for reads that represent either noise or
    the genotype of the sample at each position

    :param: sample_id str - Sample ID to be used in output files
    :param: bam_file_path str - Path to bam file
    :param: tag str -
    :param: noise_df pd.DataFrame - Data frame with columns "chrom" and "pos"
    :param: include_readname bool -
    :param: mifs int - Minimum tlen of reads to include in calculation
    :param: mafs int - Maximum tlen of reads to include in calculation
    """
    filename = sample_id + '_noise_by_tlen.tsv'
    out_fh = open(filename, 'w')
    bamfile = pysam.AlignmentFile(bam_file_path, "rb")
    if len(noise_df.index) == 0:
        return
    for i, noise_pos in enumerate(noise_df.itertuples()):
        get_fragment_size_for_noisy_position(sample_id, bamfile, noise_pos, tag, out_fh, mifs, mafs)
    out_fh.close()
    # Read file back in so it can be returned for plotting
    noisy_tlen_df = pd.read_csv(
        filename,
        sep="\t",
        header=None,
        names=["Sample", "Type", "read_id", "Var", "Size", "Chr", "Pos", "geno_not_geno"],
        dtype="object",
    )
    return noisy_tlen_df


def get_fragment_size_for_noisy_position(sample_id, bamfile, noise_pos, tag, out_fh, mifs, mafs):
    """
    The TLEN of each of the paired reads will count as either GENOTYPE or NOISE,
    depending on the base at the query position.

    Called once for each position.

    :param sample_id: str
    :param bamfile: pysam.AlignmentFile
    :param noise_pos: pd.Series - Query position, with 'chrom' and 'pos' columns
    :param tag: str
    :param out_fh: filehandle
    :param mifs: int
    :param mafs: int
    """
    start = int(noise_pos.pos) - 300
    if start < 1:
        start = '1'
    start = str(start)
    end = str(int(noise_pos.pos) + 300)
    region_string = str(noise_pos.chrom) + ":" + start + "-" + end

    for read1, read2 in read_pair_generator(bamfile, region_string=region_string):
        base_counts = {'A': noise_pos.A, 'C': noise_pos.C, 'G': noise_pos.G, 'T': noise_pos.T}
        genotype = max(base_counts, key=base_counts.get)
        non_geno_bases = ['A', 'C', 'G', 'T']
        non_geno_bases.remove(genotype)

        geno_coordinate = "\t".join(
            [str(noise_pos.chrom), str(noise_pos.pos), genotype])

        non_geno_coordinate = "\t".join(
            [str(noise_pos.chrom), str(noise_pos.pos), 'not_' + genotype])

        try:
            tlen_is_within_range = (mifs <= abs(read1.tlen) <= mafs)
            if not tlen_is_within_range: continue

            position_is_within_read_1 = (read1.pos <= int(noise_pos.pos) <= read1.aend)
            position_is_within_read_2 = (read2.pos <= int(noise_pos.pos) <= read2.aend)

            read_1_base = read1.seq[int(noise_pos.pos) - read1.pos]
            read_2_base = read2.seq[int(noise_pos.pos) - read2.pos]
            read_1_has_genotype_at_position = (read_1_base == genotype)
            read_2_has_genotype_at_position = (read_2_base == genotype)
            read_1_has_non_geno_at_position = (read_1_base in non_geno_bases)
            read_2_has_non_geno_at_position = (read_2_base in non_geno_bases)
            read_1_has_n_at_position = (read_1_base == 'N')
            read_2_has_n_at_position = (read_2_base == 'N')

            if ((position_is_within_read_1 and read_1_has_non_geno_at_position)
                or
                (position_is_within_read_2 and read_2_has_non_geno_at_position)):

                # Todo: Handle overlapping reads with different bases
                sub_type = genotype + '>' + (read_1_base if position_is_within_read_1 else read_2_base)
                line = [sample_id, tag, read1.qname, sub_type, str(abs(read1.tlen)), non_geno_coordinate]
                out_fh.write("\t".join(line) + "\n")
                continue

            elif ((position_is_within_read_1 and read_1_has_genotype_at_position)
                or
                (position_is_within_read_2 and read_2_has_genotype_at_position)):

                line = [sample_id, tag, read1.qname, "GENOTYPE", str(abs(read1.tlen)), geno_coordinate]
                out_fh.write("\t".join(line)+ "\n")
                continue

            elif ((position_is_within_read_1 and read_1_has_n_at_position)
                or
                (position_is_within_read_2 and read_2_has_n_at_position)):

                line = [sample_id, tag, read1.qname, "N", str(abs(read1.tlen)), geno_coordinate]
                out_fh.write("\t".join(line)+ "\n")
                continue

        except (IndexError, AttributeError, TypeError) as e:
            # todo: investigate why TypeError is thrown and read1.aend is sometimes undefined
            print(e)
            pass
