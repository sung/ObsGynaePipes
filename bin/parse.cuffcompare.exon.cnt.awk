BEGIN{FS="\t"}{split($9,fields,";"); split(fields[2], d, "transcript_id ");split(d[2], f, "\""); Cnt[f[2]]++}END{for(c in Cnt)print c,Cnt[c]}
