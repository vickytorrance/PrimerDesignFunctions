__author__ = 'victoriatorrance'


import Primers
import pickle
from intermine.webservice import Service
import readingAPEfunctions as AP


################ USER INPUT REQUIRED HERE################

#If constructing a new plasmid and inserting a sequence into a vector by homologous recombination read in the plasmid here
#MyPlasmid = open(r"pDL1728.ape").read()

# The tm of the primers can be changed here.
PrimerTM = 55

# IF the PCR product will contain 'homology tails' insert them here
homology_F = ''
homology_R = ''


# If homology regions for deletion or tagging a DNA seq are being designed
# the part of the primer used for amplification can be inserted here
amplify_F = ''
amplify_R = ''

#########################################################

PrimerTypeAns = raw_input('Would you like to... \n a: design primers to be used to delete a region of DNA (homology tails) \n b: just amplify the DNA\n c: incorporate a TAG\n').lower()

if PrimerTypeAns == 'a':
    PrimerType = 'homology'
elif PrimerTypeAns == 'b':
    PrimerType = 'amplification'
elif PrimerTypeAns == 'c':
    PrimerPosition = 'TAG'
    TAG = raw_input('a: N-Terminal (flanking the ATG)\nb: C-Terminal (flanking the stop codon)\n').lower()

if PrimerTypeAns == 'a' or PrimerTypeAns == 'b' :
    print 'Where would you like the primers to anneal?\n'

    ans_PrimerPosition = raw_input(' a: The promoter of a given ORF \n b: Specific distances flanking the ORF\n c: ORF\n').lower()

    if ans_PrimerPosition == 'a':
        PrimerPosition = 'Promoter'
    elif ans_PrimerPosition == 'b':
        PrimerPosition = 'Flanking'
        upstream = int(raw_input('Upstream flanking distance?\n'))
        downstream = int(raw_input('Downstream flanking distance?\n'))
    elif ans_PrimerPosition == 'c':
        PrimerPosition = 'ORF'

if PrimerTypeAns == 'b':
    writeAPE = raw_input('Create an APE file of the PCR product? y/n\n').lower()
    plas = raw_input('If cloning the product into a vector would you like an APE file of the PCR product in your plasmid? y/n\n').lower()

direction = raw_input('Do you want the insert in the opposite direction? y/n\n').lower()

genes = open(r'input')
promoters = open(r"yeast promoter sizes from YPA.txt")

print 'Please wait a wee moment for your sequences to appear...'
service = Service("http://yeastmine.yeastgenome.org/yeastmine/service")

# Get a new query on the class (table) you will be querying:
query = service.new_query("Gene")

# The view specifies the output columns
query.add_view("secondaryIdentifier", "symbol", "length", "flankingRegions.direction",
    "flankingRegions.sequence.length", "flankingRegions.sequence.residues")

# Uncomment and edit the line below (the default) to select a custom sort order:
# query.add_sort_order("Gene.secondaryIdentifier", "ASC")

# You can edit the constraint values below
query.add_constraint("Gene", "IN", "Verified_ORFs", code = "B")
query.add_constraint("flankingRegions.direction", "=", "both", code = "C")
query.add_constraint("flankingRegions.distance", "=", "1.0kb", code = "A")
query.add_constraint("flankingRegions.includeGene", "=", "true", code = "D")


mydict2 = {}
# created a dict which contains the systematic gene name as the key and then
# item[0] = short gene name, item[1] = the upstream sequence, item[2] = the ORF sequence
for row in query.rows():
    seq = row["flankingRegions.sequence.residues"]
    upstreamSeq = seq[0:1000]
    downstreamSeq = seq[-1000:]
    interseq = seq[1000:]
    ORFseq = interseq[:-1000]
    mydict2[row["secondaryIdentifier"]] = row["symbol"],upstreamSeq,downstreamSeq, ORFseq, seq
    x = row["flankingRegions.sequence.residues"]

with open('tempthings.pickle', 'wb') as handle:
  pickle.dump(mydict2, handle)

promoterdict = {}

with open('tempthings.pickle', 'rb') as handle:
  mydict = pickle.load(handle)

g = open(r"output_sequences.txt", "w")

for line in promoters:
    linestring = line.rstrip()
    linelist = linestring.split('\t')
    genename = linelist[0]
    genename2 = genename.split (' ')
    genename3 = genename2[0]
    promoter = linelist[3]
    promoterdict[genename3] = promoter

if PrimerPosition =='Promoter':
    for gene1 in genes:
        gene = gene1.rstrip('\n')
        gene = gene.rstrip('\r')
        upstreamseq = str(mydict[gene][1])
        try:
            promoter_len =  int(promoterdict[gene])
            promoter_seq = upstreamseq [-promoter_len:]
            if direction == 'y':
                promoter_seq = Primers.reverseComp(promoter_seq)
        #    Remember to include if sequence is to be reverse orientation 
        #  promoter_seq = Primers.reverseComp(promoter_seq)
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(promoter_seq, PrimerTM)
                primers =  homology_F + primers[0], homology_R + primers[1]
                if writeAPE == 'y':
                    AP.SaveToApe(promoter_seq, str(mydict[gene][0])+'   '+gene)
                if plas == 'y':
                    newPlasmid = AP.replaceAPEseq(MyPlasmid, promoter_seq, homology_F , homology_R ,str(mydict[gene][0]) )
                    newPlasmid = AP.insertFeature(newPlasmid, promoter_seq, label = str(mydict[gene][0]), colour = 'cyan')
                    AP.SaveToApe(newPlasmid, str(mydict[gene][0]))
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(promoter_seq)
                primers =  primers[0] + amplify_F, primers[1] + amplify_R

            print primers[0], primers[1]
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')
        except KeyError:
            print "we don't know the promoter length of %s"%gene


elif PrimerPosition =='Flanking':
        for gene1 in genes:
            gene = gene1.rstrip('\n')
            gene = gene.rstrip('\r')
            up1000_seq =  str(mydict[gene][1])
            up1000_seq = up1000_seq[-upstream:]
            down1000_seq =  str(mydict[gene][2])
            down1000_seq = down1000_seq[:downstream]
            mySeq = up1000_seq +str(mydict[gene][3])+  down1000_seq
            if direction == 'y':
                mySeq = Primers.reverseComp(mySeq)
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(mySeq, PrimerTM)
                primers = homology_F + primers[0], homology_R + primers[1]
                if writeAPE == 'y':
                    AP.SaveToApe(mySeq, str(mydict[gene][0])+'   '+gene)
                if plas == 'y':
                    newPlasmid = AP.replaceAPEseq(MyPlasmid, mySeq, homology_F , homology_R, str(mydict[gene][0]) )
                    newPlasmid = AP.insertFeature(newPlasmid, mySeq, label = str(mydict[gene][0]), colour = 'cyan')
                    AP.SaveToApe(newPlasmid, str(mydict[gene][0]))
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(mySeq)
                primers =  primers[0] + amplify_F, primers[1] + amplify_R
            print str(mydict[gene][0]), gene, primers
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')


elif PrimerPosition =='ORF':
        for gene1 in genes:
            gene = gene1.rstrip('\n')
            gene = gene.rstrip('\r')
            SEQorf = str(mydict[gene][3])
            if direction == 'y':
                SEQorf = Primers.reverseComp(SEQorf)
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(SEQorf, PrimerTM)
                primers = homology_F + primers[0], homology_R + primers[1]
                if writeAPE == 'y':
                    AP.SaveToApe(SEQorf, str(mydict[gene][0])+'   '+gene)
                if plas == 'y':
                    newPlasmid = AP.replaceAPEseq(MyPlasmid, SEQorf, homology_F , homology_R ,str(mydict[gene][0]) )
                    newPlasmid = AP.insertFeature(newPlasmid, SEQorf, label = str(mydict[gene][0]), colour = 'cyan')
                    AP.SaveToApe(newPlasmid, str(mydict[gene][0]))
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(SEQorf)
                primers =   primers[0] + amplify_F, primers[1] + amplify_R
            print str(mydict[gene][0]), gene,primers
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')


elif PrimerPosition == 'TAG':
        for gene1 in genes:
            gene = gene1.rstrip('\n')
            gene = gene.rstrip('\r')
            down1000_seq =  str(mydict[gene][2])
            down1000_seq =  down1000_seq[:40]
            up1000_seq =  str(mydict[gene][1])
            up1000_seq = up1000_seq[-40:]
            if TAG == 'b' or TAG == 'c':
                SEQorf = str(mydict[gene][3])
                SEQorf = SEQorf[-43:]
                mySeq = SEQorf + down1000_seq
            elif TAG == 'a' or TAG == 'n':
                SEQorf = str(mydict[gene][3])
                SEQorf = SEQorf[:43]
                mySeq = up1000_seq + SEQorf
            primers = Primers.designHomologyPair(mySeq)
            primers =  primers[0] + amplify_F, primers[1] + amplify_R

            print str(mydict[gene][0]),gene,primers
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')

