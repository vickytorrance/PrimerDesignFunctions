__author__ = 'victoriatorrance'

import Primers
import pickle
from intermine.webservice import Service

homology_A = ''
homology_B = ''
amplify_A = ''
amplify_B = ''

homology_vector = 'ccctcactaaagggaacaaaagctggagctccaccgcgt'
homology_renilla = 'tccgtttcctttgttctggatcataaactttcgaagtcat'

PrimerTypeAns = raw_input('Would you like to... \n a: design primers to be used to delete a region of DNA \n b: just amplify the DNA\n c: incorporate a TAG\n').lower()

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

genes = open(r'input')

#promoters = open (r'C:\Users\b2041087\Dropbox\PrimerDesignFunctions\yeast promoter sizes from YPA.txt')
#promoters = open(r"/Users/victoriatorrance/Desktop/yeast promoter sizes from YPA.txt")
promoters = open(r"yeast promoter sizes from YPA.txt")

print 'Please wait a wee moment for your sequences to apper...'
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
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(promoter_seq, 55)
                primers =  homology_A + primers[0], homology_B + primers[1]
                if writeAPE == 'y':
                    with open(r"PCR Product %s  %s.ape"%(gene, str(mydict[gene][0])), "w") as fp:
                        fp.write(promoter_seq)
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(promoter_seq)
                primers =  primers[0] + amplify_A, primers[1] + amplify_B
            print homology_vector + primers[0], homology_renilla + primers[1]
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')

            #with open(r"C:\Users\b2041087\Dropbox\PrimerDesignFunctions\promoter_seqs_%s.ape"%(gene), "w") as fp:
             #   fp.write(promoter_seq)

        except KeyError:
            print "we don't know the promoter length of %s :("%gene



elif PrimerPosition =='Flanking':
        for gene1 in genes:
            gene = gene1.rstrip('\n')
            gene = gene.rstrip('\r')
            up1000_seq =  str(mydict[gene][1])
            up1000_seq = up1000_seq[-upstream:]
            down1000_seq =  str(mydict[gene][2])
            down1000_seq = down1000_seq[:downstream]
            mySeq = up1000_seq +str(mydict[gene][3])+  down1000_seq
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(mySeq, 55)
                primers = homology_A + primers[0], homology_B + primers[1]
                if writeAPE == 'y':
                    with open(r"PCR Product %s  %s.ape"%(gene, str(mydict[gene][0])), "w") as fp:
                        fp.write(mySeq)
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(mySeq)
                primers =  primers[0] + amplify_A, primers[1] + amplify_B
            print str(mydict[gene][0]), gene, primers
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')


elif PrimerPosition =='ORF':
        for gene1 in genes:
            gene = gene1.rstrip('\n')
            gene = gene.rstrip('\r')
            SEQorf = str(mydict[gene][3])
            if PrimerType == 'amplification':
                primers = Primers.designPrimerpair(SEQorf, 55)
                primers = homology_A + primers[0], homology_B + primers[1]
                if writeAPE == 'y':
                    with open(r"PCR Product %s  %s.ape"%(gene, str(mydict[gene][0])), "w") as fp:
                        fp.write(SEQorf)
            elif PrimerType =='homology':
                primers = Primers.designHomologyPair(SEQorf)
                primers =   primers[0] + amplify_A, primers[1] + amplify_B
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
            primers =  primers[0] + amplify_A, primers[1] + amplify_B

            print str(mydict[gene][0]),gene,primers
            g.write(str(mydict[gene][0])+'\t'+ gene +'\t'+ str(primers) + '\n')


# also...
# the promoter and a specific distance along the ORF
# the promoter and stop codon
# the promoter and a downstream region of DNA
# a specific region of the ORF and another specific region of the ORF
# a specific region of the ORF and a downstream flanking distance
# an upstream flanking region and a specific region in the DNA