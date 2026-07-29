#!/usr/bin/env bash

### The script identifies gene sub-blocks that have a single random-id feature
### label but propagated coordinates in more than one genome. These are not
### "original singletons" in a pre-panda sense. They are post-panda rows where a
### lone random-id label marks a sequence now found in multiple genomes.
###
### For each matching row, the script extracts every propagated interval and its
### length. For example, AAB#0#chrXVI:207883-209683 has length 1800.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
file="${script_dir}/../sts/gene/sts-by-genomes.txt"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}

### identify relevant columns
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "Class_id")
            classCol = i

        if ($i == "Features_id")
            featCol = i

        if ($i == "N_feats")
            nFeatsCol = i

        if ($i == "N_feats_rid")
            nFeatsRidCol = i

        if ($i == "Ν_pres")
            nPresCol = i

        if ($i == "F_pres")
            fPresCol = i
    }

    if (!classCol || !featCol || !nFeatsCol || !nFeatsRidCol ||
        !nPresCol || !fPresCol) {
        print "Error: missing required column in " FILENAME > "/dev/stderr"
        exit 1
    }

    nGenomeCols = 0
    for (i = 3; i < nPresCol; i++) {
        nGenomeCols++
        genomeCols[nGenomeCols] = i
        genomeNames[i] = $i
    }

    print "Features_id", "N_pres", "F_pres", "Genome_id",
          "Interval_id", "Size"

    next
}

### consider only gene rows
$classCol == "gene" {

    ### keep rows with a single random-id feature label and propagated
    ### coordinates in more than one genome
    if (($nFeatsCol + 0) != 1 ||
        ($nFeatsRidCol + 0) != 1 ||
        ($nPresCol + 0) <= 1)
        next

    recoveredRidRows++
    nPres = $nPresCol + 0
    fPres = $fPresCol + 0

    for (g = 1; g <= nGenomeCols; g++) {
        col = genomeCols[g]

        if ($col == "")
            continue

        nIntervals = split($col, intervals, ";")

        for (i = 1; i <= nIntervals; i++) {
            intervalId = intervals[i]

            if (intervalId == "")
                continue

            coord = intervalId
            sub(/^.*:/, "", coord)

            if (split(coord, pos, "-") != 2)
                continue

            size = (pos[2] + 0) - (pos[1] + 0)
            recoveredRidIntervals++
            totalSize += size

            if (recoveredRidIntervals == 1 ||
                size < minSize)
                minSize = size

            if (size > maxSize)
                maxSize = size

            print $featCol, nPres, fPres, genomeNames[col],
                  intervalId, size
        }
    }
}

END {
    print "summary", "candidate_rows", recoveredRidRows + 0,
          "", "", "" > "/dev/stderr"
    print "summary", "candidate_intervals", recoveredRidIntervals + 0,
          "", "", "" > "/dev/stderr"

    if (recoveredRidIntervals > 0) {
        printf "summary\tmin_size\t%d\t\t\t\n", minSize > "/dev/stderr"
        printf "summary\tmean_size\t%.2f\t\t\t\n",
               totalSize / recoveredRidIntervals > "/dev/stderr"
        printf "summary\tmax_size\t%d\t\t\t\n", maxSize > "/dev/stderr"
    }
}
' "$file"
