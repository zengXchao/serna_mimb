#!/usr/bin/env python
#-*- coding:utf-8 -*-

import sys
import os
import getopt
import tempfile

Usage = """
Merge exons to create a unique set of exons for each gene based on a GENCODE GTF file
=====================================================================================
\x1b[1mUSAGE:\x1b[0m
  %s --gtf FILE [--includeIntron -o FILE]

\x1b[1mHELP:\x1b[0m
  -g,--gtf            <String>
                        Input a GENCODE GTF file
\x1b[1mOptions\x1b[0m
  -i,--includeIntron  <None>
                        Include intron features
  -o,--out            <String>
                        Output a filtered or modified GTF file (default: stdout)
  -h,--help           <None>
                        Print usage

\x1b[1mAUTHOR:\x1b[0m
  Chao Zeng
""" % (sys.argv[0])

Params = { 'outGTF': None, 'inGTF': None, 'includeIntron': False, 'tmpDir': None, "tidDuplicateCheck":{}}

def elog(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def init():	
	opts, args = getopt.getopt(sys.argv[1:], 'hig:o:a', ['out=', 'gtf=', 'includeIntron', 'help'])

	for op, val in opts:
		if op in ["-h", "--help"]:
			elog(Usage)
			exit(-1)
		elif op in ["-g", "--gtf"]:
			Params['inGTF'] = val
		elif op in ["-o", "--out"]:
			Params['outGTF'] = val
		elif op in ["-i", "--includeIntron"]:
			Params['includeIntron'] = True
	
	if Params['inGTF'] == None:
		elog("[ERROR] GTF file is not specified!")
		elog(Usage)
		exit(-1)
	if Params['outGTF'] == None:
		Params['outGTF'] = sys.stdout
	else:
		Params['outGTF'] = open(Params['outGTF'], "w")

	Params['tmpDir'] = tempfile.mkdtemp(dir=os.getcwd())
	return Params

def cleanup():
	os.system('rm -rf %s' % (Params['tmpDir']))
	if Params['outGTF'] != sys.stdout:
		Params['outGTF'].close()

def parseGTF(fpath):
	gtf_record = []
	with open(fpath) as fh:
		for line in fh:	
			if line[0] == "#": continue
			line = line.rstrip()
			a = line.split("\t")
			b = a[-1].split()
			att = {}
			for i in range(0,len(b),2):
				att[b[i]] = b[i+1][1:-2]
			record = {}
			record["seqname"], record["source"], record["feature"], record["start"], record["end"], record["strand"], record["attributes"] = a[0], a[1], a[2], int(a[3]), int(a[4]), a[6], att
			record["line"] = line
			# att["startInt"], att["endInt"] = int(a[3]),int(a[4])
			# att["line"] = line
			# att["main"] = "\t".join(a[:-1])
			# att["attributes"] = a[-1]
			gtf_record.append(record)
	return gtf_record

def isAbundantRNA(biotype_type):
	# Reference to https://www.gencodegenes.org/pages/biotypes.html
	if biotype_type in ["Mt_rRNA","Mt_tRNA","miRNA","rRNA","snoRNA",
		"Mt_tRNA_pseudogene","tRNA_pseudogene","snoRNA_pseudogene","rRNA_pseudogene","miRNA_pseudogene",]:
		return True
	else:
		return False

def calTranscriptLength(gtf_record):
	tid2Len = {}
	for record in gtf_record:
		if record["feature"] == "exon":
			tid, exonLen = record["attributes"]["transcript_id"], record["end"]-record["start"]+1
			if exonLen < 1:
				elog("[ERROR] Invaliaded length of exon")
				elog(record)
				cleanup()
				exit(-1)
			if tid not in tid2Len: tid2Len[tid] = 0
			tid2Len[tid] += exonLen
	return tid2Len

def get_exons_and_introns(exons, strand):
	mergedExons, introns = [], []
	tmpExonBed = "%s/exons.bed" % Params["tmpDir"]
	tmpSortedExonBed = "%s/exons.sorted.bed" % Params["tmpDir"]
	tmpMergedExonBed = "%s/exons.merged.bed" % Params["tmpDir"]
	with open(tmpExonBed, "w") as fh:
		for x in exons:
			fh.write(x + "\n")
	os.system("sort -k1,1 -k2,2n %s > %s"%(tmpExonBed, tmpSortedExonBed))
	os.system("bedtools merge -d 1 -i %s > %s"%(tmpSortedExonBed, tmpMergedExonBed))

	with open(tmpMergedExonBed) as fh:
		for line in fh:
			line = line.strip()
			if len(line) == 0: continue
			x = line.split()
			
			try:
				s, e = int(mergedExons[-1][1]) + 1, int(x[1]) - 1
				introns.append([str(s), str(e)])
			except:
				pass

			mergedExons.append(x[1:])

	if strand == "-":
		mergedExons.reverse()
		introns.reverse()

	return mergedExons, introns

def get_exons_and_introns(exons, strand):
	mergedExons, introns = [], []
	tmpExonBed = "%s/exons.bed" % Params["tmpDir"]
	tmpSortedExonBed = "%s/exons.sorted.bed" % Params["tmpDir"]
	tmpMergedExonBed = "%s/exons.merged.bed" % Params["tmpDir"]
	with open(tmpExonBed, "w") as fh:
		for x in exons:
			fh.write(x + "\n")
	os.system("sort -k1,1 -k2,2n %s > %s"%(tmpExonBed, tmpSortedExonBed))
	os.system("bedtools merge -d 1 -i %s > %s"%(tmpSortedExonBed, tmpMergedExonBed))

	with open(tmpMergedExonBed) as fh:
		for line in fh:
			line = line.strip()
			if len(line) == 0: continue
			x = line.split()
			
			try:
				s, e = int(mergedExons[-1][1]) + 1, int(x[1]) - 1
				introns.append([str(s), str(e)])
			except:
				pass

			mergedExons.append(x[1:])

	if strand == "-":
		mergedExons.reverse()
		introns.reverse()

	return mergedExons, introns

def print_buf(buf):
	txExons, strand, chrom, gid, gname, tid = [], None, None, None, None, None
	for x in buf:
		a, attr = x[0], x[1]
		if a[2] == "gene":
			chrom, strand = a[0], a[6]
			gid, gname = attr["gene_id"], attr["gene_name"]
			tid = gname
			if tid in Params["tidDuplicateCheck"]: tid = gname + "_" + gid
			Params["tidDuplicateCheck"][tid] = 1
			a[1], a[2] = "COLLAPSE", "transcript"
			Params["outGTF"].write("\t".join(a) + "\tgene_id \"%s\"; transcript_id \"%s\"; transcript_biotype \"%s\";\n" % (gid, tid, attr["gene_type"]))
		elif a[2] == "exon":
			txExons.append("\t".join([a[0], a[3], a[4]]))

	exons, introns = get_exons_and_introns(txExons, strand)
	
	for i,x in enumerate(exons):
		Params["outGTF"].write("\t".join([chrom, "COLLAPSE", "exon", x[0], x[1], ".", strand, "."]) + "\tgene_id \"%s\"; transcript_id \"%s\"; exon_id \"%s_exon_%d\";\n" % (gid, tid, tid, (i+1)))
		if Params["includeIntron"]:
			try:
				Params["outGTF"].write("\t".join([chrom, "COLLAPSE", "exon", introns[i][0], introns[i][1], ".", strand, "."]) + "\tgene_id \"%s\"; transcript_id \"%s\"; intron_id \"%s_intron_%d\";\n" % (gid, tid, tid, (i+1)))
			except:
				pass		
	buf[:] = []


def buf_per_gene():
	gid, buf = None, []
	with open(Params['inGTF']) as fh:
		for line in fh:
			if line[0] == "#": continue
			line = line.strip()
			x = line.split("\t")
			
			attr = {}
			for y in x[-1].split(";"):
				y = y.strip()
				blankIdx = y.find(" ")
				k, v = y[:blankIdx], y[blankIdx+2:-1]
				if len(k) == 0: continue
				attr[k] = v

			if gid and gid != attr["gene_id"]: print_buf(buf)

			gid = attr["gene_id"]
			if not isAbundantRNA(attr["gene_type"]):
				buf.append([x[:-1], attr])

	print_buf(buf)


if __name__ == "__main__":	
	init()
	buf_per_gene()
	cleanup()






