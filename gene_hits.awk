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