#!/usr/bin/env bash

### The script identifies genes that were originally annotated as private 
### and determines whether graph-based feature propagation recovered homologous 
### coordinates in additional strains. It dynamically identifies the relevant 
### columns from the header, counts non-empty strain columns for each 
### original private gene, and reports how many remain private or become 
### shared accessory genes, together with their percentages.

file="sts-by-genomes.txt"

awk -F '\t' '
BEGIN {
    OFS = "\t"
}

### identify columns from the header
NR == 1 {
    for (i = 1; i <= NF; i++) {
        if ($i == "Class_id")
            classCol = i

        if ($i == "Features_id")
            featCol = i

        /*
         * The column immediately before F_pres is Ν_pres,
         * and the preceding column is the last strain.
         */
        if ($i == "F_pres")
            lastStrainCol = i - 2
    }

    next
}

### consider only gene rows
$classCol == "gene" {

    ### count the original annotation labels
    nLabels = split($featCol, labels, ",")

    ### keep only genes originally classified as private
    if (nLabels != 1)
        next

    originalSingletons++

    ### count strains with at least one propagated coordinate
    nStrains = 0

    for (i = featCol + 1; i <= lastStrainCol; i++) {
        if ($i != "")
            nStrains++
    }

    if (nStrains == 1)
        remainingSingletons++
    else if (nStrains >= 2)
        sharedAccessory++
    else
        noCoordinates++
}

END {
    print "category", "count"

    print "original_singletons",
          originalSingletons + 0

    print "remain_singletons",
          remainingSingletons + 0

    print "shared_accessory",
          sharedAccessory + 0

    print "no_coordinates",
          noCoordinates + 0

    if (originalSingletons > 0) {
        printf "remain_singletons_pct\t%.2f\n",
               100 * remainingSingletons / originalSingletons

        printf "shared_accessory_pct\t%.2f\n",
               100 * sharedAccessory / originalSingletons
    }
}
' "$file"
