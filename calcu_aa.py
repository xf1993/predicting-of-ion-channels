#!/usr/bin/python
# -*- coding: utf-8 -*-

from Bio import SeqIO
a=SeqIO.parse(r"D:\yezo\课件\TEAM\456.txt","fasta")#打开要处理的氨基酸序列文件
outfile=open(r"D:\yezo\课件\TEAM\夏凤数据挖掘.txt","w")#输入要输出的氨基酸频率文件地址
outfile.write("\t\t\tAla\tArg\tAsn\tAsp\tCys\tGln\tGlu\tGly\tHis\tIle\tLeu\tLys\tMet\tPhe\tPro\tSer\tTrp\tTyr\tVal\tThr"+"\n")
ID=""
for seq_record in a:
	if len(seq_record.seq)>=30:
		ID=seq_record.id
		Ala = str(seq_record.seq.count("A")/len(seq_record.seq))
		Arg = str(seq_record.seq.count("R")/len(seq_record.seq))
		Asn = str(seq_record.seq.count("N")/len(seq_record.seq))
		Asp = str(seq_record.seq.count("D")/len(seq_record.seq))
		Cys = str(seq_record.seq.count("C")/len(seq_record.seq))
		Gln = str(seq_record.seq.count("Q")/len(seq_record.seq))
		Glu = str(seq_record.seq.count("E")/len(seq_record.seq))
		Gly = str(seq_record.seq.count("G")/len(seq_record.seq))
		His = str(seq_record.seq.count("H")/len(seq_record.seq))
		Ile = str(seq_record.seq.count("I")/len(seq_record.seq))
		Leu = str(seq_record.seq.count("L")/len(seq_record.seq))
		Lys = str(seq_record.seq.count("K")/len(seq_record.seq))
		Met = str(seq_record.seq.count("M")/len(seq_record.seq))
		Phe = str(seq_record.seq.count("F")/len(seq_record.seq))
		Pro = str(seq_record.seq.count("P")/len(seq_record.seq))
		Ser = str(seq_record.seq.count("S")/len(seq_record.seq))
		Trp = str(seq_record.seq.count("W")/len(seq_record.seq))
		Tyr = str(seq_record.seq.count("Y")/len(seq_record.seq))
		Val = str(seq_record.seq.count("V")/len(seq_record.seq))
		Thr = str(seq_record.seq.count("T")/len(seq_record.seq))
		outfile.write("\t"+Ala+"\t"+Cys+"\t"+Asp+"\t"+Glu+"\t"+Phe+"\t"+Gly+"\t"+His+"\t"+Ile+"\t"+Lys+"\t"+Leu+"\t"+Met+"\t"+Asn+"\t"+Pro+"\t"+Gln+"\t"+Arg+"\t"+Ser+"\t"+Thr+"\t"+Val+"\t"+Trp+"\t"+Tyr+"\n")
		#输出的是每个序列的氨基酸频率，求整个基因组的，需要调整
		
