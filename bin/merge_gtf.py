import os
import sys
import tempfile
import math

def elog(msg):
	print(msg, file=sys.stderr)

def get_chr_len():
	Params["chrLen"] = {}
	with open(Params["fpathChrNameLen"]) as fh:
		for line in fh:
			a = line.rstrip().split()
			Params["chrLen"][a[0]] = int(a[1])

def parse_GTF(fpath):
	elog("Parsing: " + fpath)
	gtf = []
	with open(fpath) as fh:
		for line in fh:
			if line[0] == "#": continue
			line = line.rstrip()
			a = line.split("\t")
			b = a[-1].split()
			att = {}
			for i in range(0,len(b),2):
				att[b[i]] = b[i+1][1:-2]
			att["seqname"], att["source"], att["feature"], att["start"], att["end"], att["strand"] = a[0], a[1], a[2], a[3], a[4], a[6]
			att["startInt"], att["endInt"] = int(a[3]),int(a[4])
			att["line"] = line
			att["main"] = "\t".join(a[:-1])
			att["attributes"] = a[-1]
			gtf.append(att)
	return gtf

def gtf_to_bed(fpathGtf, fpathBed, feature="gene", nameField="gene_name", flanking = 0, nameFieldKey=None):
	bedList = []
	for x in parse_GTF(fpathGtf):
		a, b = x["startInt"]-flanking, x["endInt"]+flanking
		if a < 0 : continue
		if x["seqname"] not in Params["chrLen"]: continue
		if b > Params["chrLen"][x["seqname"]]: continue
		if nameFieldKey:
			if x["feature"] == feature and nameFieldKey in x[nameField]: bedList.append( (x["seqname"], a, b, x[nameField], "0", x["strand"]) )
		elif x["feature"] == feature: 
			bedList.append( (x["seqname"], a, b, x[nameField], "0", x["strand"]) )
	bedList.sort(key=lambda x: (x[0], x[1]))
	OUT = open(fpathBed,"w") 
	for bed in bedList:
		OUT.writelines("%s\t%s\t%s\t%s\t%s\t%s\n" % tuple(bed))
	OUT.close()

def intersect_bed(inBed1, inBed2, outTxt):
	elog("Intersection: " + outTxt)
	os.system("bedtools intersect -nonamecheck -a %s -b %s -wa -wb > %s" % (inBed1, inBed2, outTxt))

def get_unannoated_gene():
	
	gencodeGeneBedFile = Params["tmpDir"] + "/genecodeGene.bed"
	stringtieGene5kBedFile = Params["tmpDir"] + "/stringtieGene5k.bed"
	stringtieGene5kXgencodeGeneFile = Params["tmpDir"] + "/stringtieGene5k_x_genecodeGene.txt"
	
	gtf_to_bed(Params["fpathGencodeCollapseGtf"], gencodeGeneBedFile, "transcript", "transcript_id", 0)
	gtf_to_bed(Params["fpathStringtieGtf"], stringtieGene5kBedFile, "exon", "gene_id", 5000)
	intersect_bed(stringtieGene5kBedFile, gencodeGeneBedFile, stringtieGene5kXgencodeGeneFile)

	stringtieGene = {}
	with open(stringtieGene5kBedFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			stringtieGene[a[3]] = a[0]+":"+a[1]+"-"+a[2]+"\t"+a[-1]
	stringtieGeneNum = len(stringtieGene)
	
	with open(stringtieGene5kXgencodeGeneFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			if a[3] in stringtieGene and a[5] == a[11]: del stringtieGene[a[3]]

	antisense = {}
	with open(stringtieGene5kXgencodeGeneFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			if a[3] in stringtieGene and a[5] != a[11]: antisense[a[3]] = a[-3]

	intergenic = {}
	for x in stringtieGene:
		if x not in antisense: intergenic[x] = stringtieGene[x]
	
	return antisense,intergenic
	
def get_intron_rs():
	retainedIntronByRS = {}
	with open(Params["fpathIntronRS"]) as fh:
		# fpkmThreshold = float(fh.readline().split()[-3])
		fh.readline()
		fh.readline()
		for line in fh:
			a = line.rstrip().split()
			iID, iFPKM, eFPKM, RS = a[0], float(a[1]), float(a[2]), float(a[3])
			try:
				# if math.log(eFPKM,10) >= fpkmThreshold and math.log(iFPKM,10) >= fpkmThreshold and RS >= 0.05:
				# if math.log(iFPKM,10) >= fpkmThreshold and RS >= 0.05:
				if iFPKM >= 0.1 and RS >= 0.1:
					retainedIntronByRS[iID] = " intron_FPKM \"%.2f\"; flanking_exon_FPKM \"%.2f\"; retention_score \"%.2f\";" % (iFPKM, eFPKM, RS)
			except:
				pass
			
	return retainedIntronByRS

def filter_overlapping_intron():
	retainedIntronByRS = get_intron_rs()
	gencodeExonBedFile = Params["tmpDir"] + "/gencodeExon.bed"
	gencodeIntronBedFile = Params["tmpDir"] + "/gencodeIntron.bed"
	gencodeIntronXgencodeExonFile = Params["tmpDir"] + "/gencodeIntron_x_gencodeExon.txt"
	gtf_to_bed(Params["fpathStringtieRefGtf"], gencodeExonBedFile, "transcript", "transcript_id", 0, "_exon_")
	gtf_to_bed(Params["fpathStringtieRefGtf"], gencodeIntronBedFile, "transcript", "transcript_id", 0, "_intron_")
	intersect_bed(gencodeIntronBedFile, gencodeExonBedFile, gencodeIntronXgencodeExonFile)
	with open(gencodeIntronXgencodeExonFile) as fh:
		for line in fh:
			a = line.rstrip().split()
			if a[3] in retainedIntronByRS and a[5] == a[11]:	
				del retainedIntronByRS[a[3]]
	return retainedIntronByRS

def get_nested_gene():
	gencodeGeneBedFile = Params["tmpDir"] + "/gencodeGene.bed"
	nestGeneFile = Params["tmpDir"] + "/nestGene.bed"
	gtf_to_bed(Params["fpathGencodeCollapseGtf"], gencodeGeneBedFile, "transcript","gene_id")
	os.system("bedtools intersect -nonamecheck -a %s -b %s -wao -f 1 -s | awk '{if($4!=$10){print $4}}' | sort | uniq > %s"%(gencodeGeneBedFile,gencodeGeneBedFile,nestGeneFile))
	nestGene = {}
	with open(nestGeneFile) as fh:
		for line in fh:
			nestGene[line.rstrip()] = 0
	return nestGene

def rebuild_GTF():
	fpathMergedGTF = Params["tmpDir"] + "/merged.gtf"
	mergedGTFfh = open(fpathMergedGTF, "w")
	stringtieGTF = parse_GTF(Params["fpathStringtieGtf"])
	antisense,intergenic = get_unannoated_gene()
	antiGeneNum, intergenicGeneNum = len(antisense),len(intergenic)
	annotatedGeneNum, irGeneNum = 0, 0
	for x in stringtieGTF:
		tag = None
		if x["gene_id"] in antisense: 
			tag = " transcript_category \"antisense\"; associated_gene \"%s\"; " % antisense[x["gene_id"]]
		elif x["gene_id"] in intergenic: 
			tag = " transcript_category \"intergenic\";"
	
		if tag:
			if x["feature"] == "transcript":
				mergedGTFfh.write("%s\n" % (x["line"]+tag))
			else:
				mergedGTFfh.write("%s\n" % x["line"])

	gencodeGTF = parse_GTF(Params["fpathGencodeCollapseGtf"])
	stringtieRefGTF = parse_GTF(Params["fpathStringtieRefGtf"])
	retainedIntronByRS = filter_overlapping_intron()
	retainedIntronTx = {}
	for x in stringtieRefGTF:
		if x["feature"] == "exon" and x["transcript_id"] in retainedIntronByRS:
			txID = x["transcript_id"].split("_intron_")[0]
			if txID not in retainedIntronTx: retainedIntronTx[txID] = []
			retainedIntronTx[txID].append(x["main"] + "\tgene_id \"%s\"; transcript_id \"%s_RI\"; exon_id \"%s\";" % (x["gene_id"].split("/")[0], x["transcript_id"].split("_intron_")[0], x["transcript_id"]) + retainedIntronByRS[x["transcript_id"]])
	irGeneNum = len(retainedIntronTx)
	nestGene = get_nested_gene()
	preTid, preTxLine, preExons = None, None, []
	for x in gencodeGTF:
		if x["gene_id"] in nestGene: continue
		if x["feature"] == "transcript":
			if preTid and preTid != x["transcript_id"] and preTid in retainedIntronTx:
				mergedGTFfh.write(preTxLine)
				for xx in preExons:
					mergedGTFfh.write("%s\n" % xx)
				for xx in retainedIntronTx[preTid]:
					mergedGTFfh.write("%s\n" % xx)
		
			mergedGTFfh.write(x["main"]+"\tgene_id \"%s\"; transcript_id \"%s\"; transcript_biotype \"%s\"; transcript_category \"collapsed\"; \n" % (x["gene_id"],x["transcript_id"],x["transcript_biotype"]))

			preTxLine = x["main"]+"\tgene_id \"%s\"; transcript_id \"%s_RI\"; transcript_biotype \"%s\"; transcript_category \"intron_retained\"; \n" % (x["gene_id"],x["transcript_id"],x["transcript_biotype"])
			preTid = x["transcript_id"]
			preExons[:] = []
			annotatedGeneNum += 1
		elif x["feature"] == "exon":
			mergedGTFfh.write("%s\n" % x["line"])
			preExons.append(x["main"]+"\tgene_id \"%s\"; transcript_id \"%s_RI\"; exon_id \"%s\";" % (x["gene_id"],x["transcript_id"],x["exon_id"]))

	if preTid in retainedIntronTx:
		mergedGTFfh.write(preTxLine)
		for xx in preExons:
			mergedGTFfh.write("%s\n" % xx)
		for xx in retainedIntronTx[preTid]:
			mergedGTFfh.write("%s\n" % xx)
			
		
	mergedGTFfh.close()
	os.system("mv %s %s" % (fpathMergedGTF,Params["fpathMergedGTF"]))

	# allGeneNum = annotatedGeneNum + antiGeneNum + intergenicGeneNum + irGeneNum
	# print("\n# Statistics of merged transcripts")
	# print("All transcripts: %d" % allGeneNum)
	# print(" - collapsed: %d" % annotatedGeneNum)
	# print(" - intron_retained: %d" % irGeneNum)
	# print(" - antisense: %d" % antiGeneNum)
	# print(" - intergenic: %d" % intergenicGeneNum)

Params = {
	"fpathMergedGTF":sys.argv[1],
	"fpathChrNameLen":sys.argv[2],
	"fpathGencodeCollapseGtf":sys.argv[3],
	"fpathStringtieRefGtf":sys.argv[4],
	"fpathStringtieGtf":sys.argv[5],
	"fpathIntronRS":sys.argv[6],
	"tmpDir":tempfile.mkdtemp(dir=os.getcwd()),
}
get_chr_len()
rebuild_GTF()
os.system("rm -rf " + Params["tmpDir"])





