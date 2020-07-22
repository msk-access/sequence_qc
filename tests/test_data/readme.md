SeraCare_0-5.sam
- only 10 reads from SeraCare 0.5% bam
- introduced test changes in reads:
    - N at first position
    - 5bp insertion at 2nd position
    - 1bp deletion at 8th position
    - low mapping quality
    - low base quality
    - reversed insert size (todo: correct way to do this?)
    - todo: each of 11 possible flags (+ combinations?)


ref.fa
- first 92 bases of region from chr1 that matches reads from SeraCare_0-5.sam
- represents a "fake" chromosome of only 92 bases
