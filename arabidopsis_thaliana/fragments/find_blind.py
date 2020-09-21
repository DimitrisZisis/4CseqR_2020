pattern="CATG"
flag=""
fastaFileName="/media/dimitris/3AEC27A3EC275881/4CseqR-master/4CseqR/Library_Fragments/fragments_by_fouCseq/FourCseq_fragments_arab.fa"
fastaOutputFileName="FourCseq_fragments_arab_blindNonBlind.txt"
fastaFile=open(fastaFileName,"r")
fastaOutputFile=open(fastaOutputFileName,"w")
for aLine in fastaFile:
    aLine=aLine.rstrip()
    if aLine[0]==">":
        headerLine=aLine
        
        continue
    else:
        if pattern in aLine:
            flag="NonBlind"
            #print "1"
        else:
            flag="Blind"
            #print "0"
        #print headerLine+flag+"\n"+aLine+"\n"
        print headerLine+"\t"+flag+"\n"
        fastaOutputFile.write(headerLine+"\t"+flag+"\n")
fastaOutputFile.close()
fastaFile.close()
