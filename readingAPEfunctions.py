__author__ = 'victoriatorrance'

import re
import Primers


def AddPreamble(seq, seqNAME = '', comment = ''):
    lenDNA = len(seq)
    blankAPE = 'LOCUS       %s               %i bp ds-DNA     linear       05-JAN-2016\nDEFINITION  .\nACCESSION   \nVERSION     \nSOURCE      .\n  ORGANISM  .\nCOMMENT     \nCOMMENT     %s\nCOMMENT     lApEinfo:methylated:1\nORIGIN\n'%(seqNAME, lenDNA, comment)
    if seq.find('LOCUS') != -1:
        return seq
    else :
        new = blankAPE + seq +'\n//'
        return new


def returnSeq(file_contents):
    ''' Read an APE file and return the sequence only as a string'''
    seqRaw = file_contents.split('ORIGIN')[1]
    seqRaw = seqRaw.replace("//", "")
    seqRaw = re.sub('[\s+]', '', seqRaw)
    seqRaw = ''.join([i for i in seqRaw if not i.isdigit()])
    return seqRaw


def changeSeqNumberName(file_contents, DNAlen, DNAname):
    tag = "LOCUS       %s                 %i"%(DNAname, DNAlen)
    t = file_contents.split(' bp ds-DNA')[1]
    return tag +'bp ds-DNA'+ t


def insertSeq(file_contents, insert):
    ''' insert a DNA seq into APE file- Repaces it really'''
    beforeSeq = file_contents.split('ORIGIN')[0]
    afterSeq = file_contents.split('ORIGIN')[1]
    new = beforeSeq + 'ORIGIN\n' + insert + '\n//'
    return new


def replaceAPEseq(file_contents, insert, homologyF,  homologyR, geneName):
    rawSeq = returnSeq(file_contents)
    newSeq = insertSeqIntoSeq(rawSeq, insert, homologyF,  homologyR)
    newApe = insertSeq(file_contents, newSeq)
    newApe = changeSeqNumberName(newApe, len(newSeq), geneName)
    return newApe


def insertSeqIntoSeq(seq, insert, homologyF,  homologyR):
    '''Takes a DNA string as argument and inserts a sequence by homologous recombination'''
    seq = seq.lower()
    homologyF  = homologyF.lower()
    homologyR  = homologyR.lower()
    if seq.find(homologyF) != -1:
        beforeHomology = seq.split(homologyF)[0]
    else:
        print 'Forward homology region does not exist in plasmid'
        exit()
    homologyR_rc = Primers.reverseComp(homologyR).lower()
    if seq.find(homologyR_rc) != -1:
        afterHomology = seq.split(homologyR_rc)[1]
    else:
        print 'Reverse homology region does not exist in plasmid'
        exit()
    LenSequenceRepacedByInsert =  len(seq) - (len(beforeHomology)+len(afterHomology)+len(homologyF) + len(homologyR))
    newPlamidSeq = beforeHomology + homologyF + insert + homologyR_rc + afterHomology
    before =  len(beforeHomology)
    return newPlamidSeq, before, LenSequenceRepacedByInsert

def replaceAPEseq(file_contents, insert, homologyF,  homologyR, geneName ):
    '''Takes an ape file contents and insert (seq of DNA) as arguments and updates the ape file to include the insert'''
    rawSeq = returnSeq(file_contents)
    newSeq = insertSeqIntoSeq(rawSeq, insert, homologyF,  homologyR)[0]
    before =  insertSeqIntoSeq(rawSeq, insert, homologyF,  homologyR)[1]
    LenSequenceRepacedByInsert = insertSeqIntoSeq(rawSeq, insert, homologyF,  homologyR)[2]
    x = len(insert) -LenSequenceRepacedByInsert
    newApe = insertSeq(file_contents, newSeq)
    newApe = changeSeqNumberName(newApe, len(newSeq), geneName)
    newApe = updateMiscFeatureLocs(newApe, x, before)
    return newApe

def insertFeatureAtLoc(file_contents, start, stop, label = 'no_label', colour = 'green'):
    '''Creates highlighted sequence from locations start to stop'''
    misc ='misc_feature    %i..%i\n                     /label=%s\n                     /ApEinfo_fwdcolor=%s\n                     /ApEinfo_revcolor=%s\n                     /ApEinfo_graphicformat=arrow_data {{0 1 2 0 0 -1} {} 0}\n                     width 5 offset 0\n'%(start, stop, label, colour, colour)
    commentTag = 'COMMENT     ApEinfo:methylated:1' or 'COMMENT     ApEinfo:methylated:0'
    aboveORIGIN = file_contents.split('ORIGIN')[0]
    belowORIGIN = 'ORIGIN'+file_contents.split('ORIGIN')[1]
    aboveComment = aboveORIGIN.split(commentTag)[0] + 'COMMENT     ApEinfo:methylated:1'
    belowComment = aboveORIGIN.split(commentTag)[1]
    features = belowComment
    if len(features)>1:
        featsSplit = features.split('FEATURES             Location/Qualifiers')
        before_featureTag = featsSplit[0]+'FEATURES             Location/Qualifiers'
        after_featureTag = featsSplit[1]
        beforeMisc = after_featureTag.split('misc_feature', 1)[0]
        afterMisc = '     misc_feature' + after_featureTag.split('misc_feature', 1)[1]
        misc_and_insert = before_featureTag + beforeMisc + misc + afterMisc
    else:
        misc_and_insert = 'FEATURES             Location/Qualifiers\n' + misc
    new_contents = aboveComment + misc_and_insert + belowORIGIN
    return new_contents

def insertFeature(file_contents, seq, label = 'no_label', colour = 'cyan'):
    '''take an APE file and highlight a region of DNA (seq)'''
    seqfromfile = returnSeq(file_contents)
    x = seqfromfile.split(seq)
    first = len(x[0])+1
    last = len(x[0])+len(seq)
    newFile = insertFeatureAtLoc(file_contents, first, last, label, colour)
    return newFile


def updateMiscFeatureLocs(file_contents, lenDNA, startPos):
    '''updates the highlighting of a DNA region after insertion of DNA'''
    startPos = startPos+1
    listfindalls = re.findall("misc_feature(.*?)\n", file_contents)
    oldints1, oldints2 = [], []
    newints1, newints2 = [], []
    for i in listfindalls:
        x = re.sub('complement', '', i)
        x= re.sub('[(){}<>]', '', x)
        xsplit =  x.split('..')
        num1 = int(xsplit[0])
        oldints1.append(str(num1))
        num2 = int(xsplit[1])
        oldints2.append(str(str(num2)))
        if num1 > startPos:
            olnum1updated = num1+lenDNA
            num2updated = num2+lenDNA
        else:
            olnum1updated = num1
            num2updated = num2
        newints1.append(str(olnum1updated))
        newints2.append(str(num2updated))
    oldints=zip(oldints1, oldints2)
    newints = zip(newints1, newints2)
    myoldlines, mynewlines = [],[]
    for t in oldints:
        myline = t[0]+'..'+t[1]
        myoldlines.append(myline)
    for t in newints:
        myline = t[0]+'..'+t[1]
        mynewlines.append(myline)
    OldandNew = zip(myoldlines, mynewlines)

    for t in OldandNew:
        file_contents = re.sub(t[0], t[1], file_contents)
    return file_contents

def SaveToApe(newFile, fileName):
    with open(r"%s.ape"%(fileName), "w") as fp:
        fp.write(newFile)

