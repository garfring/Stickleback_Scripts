#Program used to find genic SNPs or SNPs within a specified distance up and downstream of the SNP. Returns list of genes and relevant SNP. 
#The fourth element of the within_range function specifies the distance: Eg: '0' for genic SNPs, 500 for genes within 500 bps of SNP.
function within_range(val, lower, upper, proximity) {
    # you can specify the "proximity" as required
    return val > lower - proximity && val < upper + proximity
}

BEGIN {
    OFS="\t"
}

$1 == id && within_range(pos, $4, $5, 0) {
    name = gensub(/.*Name=([^\t]*).*/, "\\1", 1)
    if (name ~ /[^[:space:]]+/)
        print id, pos, name
}
