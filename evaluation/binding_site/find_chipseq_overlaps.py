from pybedtools import BedTool

 
def peak_file_for_mf(mf_symbol, chipseq_dir, gene_symbol_aliases):
    for gene_symbol in [mf_symbol, *gene_symbol_aliases.get(mf_symbol, [])]: # [] is the default return when the mf_symbol is not a dictionary key; if the key exists its alias list is returned
        chipseq_peak_file = chipseq_dir / f"{gene_symbol}.bed.gz"
        if chipseq_peak_file.exists():
            return chipseq_peak_file 
    return None


def cpgs_with_peak_overlap(cpg_coordinates, chipseq_peak_file, window_around_cpg = 250):
    """Return CpGs when their [start-window_around_cpg, end+window_around_cpg) interval overlaps a peak."""
    cpg_as_bed_interval = []
    for cpg_id, cpg_coordinate in cpg_coordinates.iterrows():
        window_start = max(0, cpg_coordinate["CpG_beg"] - window_around_cpg)
        window_end = cpg_coordinate["CpG_end"] + window_around_cpg
        cpg_as_bed_interval.append([cpg_coordinate["CpG_chrm"], window_start, window_end, cpg_id])
    cpg_bed = BedTool(cpg_as_bed_interval)
    chipseq_peaks = BedTool(chipseq_peak_file)
    overlapping_intervals = cpg_bed.intersect(chipseq_peaks, u = True) # u = unique: count each CpG once, also if it overlaps more peaks
    overlapping_cpg_ids = set()
    for interval in overlapping_intervals:
        overlapping_cpg_ids.add(interval.name)
    return overlapping_cpg_ids
 