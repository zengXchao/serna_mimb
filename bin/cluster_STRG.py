import os
import sys

STRG_CLUSTER = sys.argv[1]
MERGED_GTF_WO_STR = sys.argv[2]
OUT_GTF = sys.argv[3]

strg_cluster = {}
with open(STRG_CLUSTER) as fh:
	for line in fh:
		a = line.rstrip().split()
		chrm,start,end,strand,cluster = a[0],int(a[1]),int(a[2]),a[-2],int(a[-1])
		if cluster not in strg_cluster:
			strg_cluster[cluster] = [chrm, strand, [], []]
		strg_cluster[cluster][2].append(start)
		strg_cluster[cluster][3].append(end)
for x in strg_cluster:
	strg_cluster[x][2] = min(strg_cluster[x][2])
	strg_cluster[x][3] = max(strg_cluster[x][3])

with open(OUT_GTF,"w") as fhW:
	for x in sorted(strg_cluster):
		chrm,strand,start,end = strg_cluster[x][0],strg_cluster[x][1],str(strg_cluster[x][2]+1),str(strg_cluster[x][3])
		gid,tid = "STRG."+str(x), "STRG."+str(x)+".1"
		fhW.write("\t".join([chrm,"StringTie","transcript",start,end,".",strand,".","gene_id \"%s\"; transcript_id \"%s\";"%(gid,tid)])+"\n")
		fhW.write("\t".join([chrm,"StringTie","exon",start,end,".",strand,".","gene_id \"%s\"; transcript_id \"%s\"; exon_number \"1\";"%(gid,tid)])+"\n")
	with open(MERGED_GTF_WO_STR) as fh:
		for line in fh:
			fhW.write(line)