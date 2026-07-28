#!/usr/bin/env bash

### The script identifies genes that were originally annotated as private 
### and determines whether graph-based feature propagation recovered homologous 
### coordinates in additional strains. It dynamically identifies the relevant 
### columns from the header, counts non-empty strain columns for each 
### original private gene, and reports how many remain private, become 
### shared accessory genes or near core, together with their percentages.

file="../sts/gene/sts-by-genomes.txt"

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

        if ($i == "Ν_pres")
            nPresCol = i

        if ($i == "F_pres")
            fPresCol = i
    }

    next
}

### consider only gene rows
$classCol == "gene" {

    ### count the original annotation labels
    nLabels = split($featCol, labels, ",")

    ### keep only genes originally classified as singletons
    if (nLabels != 1)
        next

    originalSingletons++

    nPres = $nPresCol + 0
    fPres = $fPresCol + 0

    if (nPres == 0) {
    	### sanity check, this should be always 0
        noCoordinates++
    }
    else if (nPres == 1) {
        remainingSingletons++
    }
    else {
        recoveredSingletons++

        if (fPres > 0.976)
            core++
        else
            sharedAccessory++
    }
}

END {
    print "category", "count"

    print "original_singletons",
          originalSingletons + 0

    print "remain_singletons",
          remainingSingletons + 0

    print "recovered_singletons",
          recoveredSingletons + 0

    print "core",
          core + 0

    print "shared_accessory",
          sharedAccessory + 0

    print "no_coordinates",
          noCoordinates + 0

    if (originalSingletons > 0) {
        printf "remain_singletons_pct\t%.2f\n",
               100 * remainingSingletons / originalSingletons

        printf "recovered_singletons_pct\t%.2f\n",
               100 * recoveredSingletons / originalSingletons

        printf "core_pct\t%.2f\n",
               100 * core / originalSingletons

        printf "shared_accessory_pct\t%.2f\n",
               100 * sharedAccessory / originalSingletons
    }

    if (recoveredSingletons > 0) {
        printf "core_among_recovered_pct\t%.2f\n",
               100 * core / recoveredSingletons

        printf "shared_accessory_among_recovered_pct\t%.2f\n",
               100 * sharedAccessory / recoveredSingletons
    }
}
' "$file"
