pattern="CATG"
flag=""
fastaFileName="/media/dimitris/3AEC27A3EC275881/4CseqR-master/4CseqR/Library_Fragments/mm10_AAGCTT_fragments.fa"
fastaOutputFileName="mm10_AAGCTT_blindNonBlind_frend_len.txt"
fastaFile=open(fastaFileName,"r")
fastaOutputFile=open(fastaOutputFileName,"w")
for aLine in fastaFile:
    aLine=aLine.rstrip()
    if aLine[0]==">":
        headerLine=aLine
        
        continue
    else:
        firstPos="NA"
        lastPos="NA"
        if pattern in aLine:
            length=len(aLine)
            firstPos=aLine.find(pattern)
            lastPos=aLine.rfind(pattern)
            fiveLen=firstPos
            threeLen=length-lastPos-len(pattern)
            flag="NonBlind"
            print "-->",length,firstPos,lastPos
        else:
            flag="Blind"
            fiveLen="NA"
            threeLen="NA"
            #print "0"
        #print headerLine+flag+"\n"+aLine+"\n"
        print headerLine+"\t"+flag+"\t"+str(fiveLen)+"\t"+str(threeLen)+"\n"
        fastaOutputFile.write(headerLine+"\t"+flag+"\t"+str(fiveLen)+"\t"+str(threeLen)+"\n")
fastaOutputFile.close()
fastaFile.close()
