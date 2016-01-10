# Vicky Torrance
#
# Library: VickyModules
# Module: Primers
#
# Description: Module for primer-related scripts

def closest(mylist, Number):
    ''' returns the index of the integer in mylist which is closest to Number '''
    aux = [abs(Number-num) for num in mylist]
    return aux.index(min(aux))

def calculateTm(g):

    '''calculates TM value from input string.'''

    numA= g.count('A')
    numT= g.count('T')
    numC= g.count('C')
    numG= g.count('G')

    numAT= numA+ numT
    numGC= numG+ numC

    tm = 4*numGC + 2*numAT

    return tm

def reverseComp(seq):

    '''Function to calculate reverse complement.'''

    seq2 = seq[::-1]
    seq2 = seq2.upper()
    seq3 = seq2.replace('A','t')
    seq4 = seq3.replace('T', 'A')
    seq5 = seq4.replace('G', 'c')
    seq6 = seq5.replace('C', 'G')
    seq7 = seq6.replace('t', 'T')
    seq8 = seq7.replace('c', 'C')

    return seq8


def designPrimerpair(seq, TM, nmin=15, nmax=30):

    '''returns forward and reverse primers. takes arguments, seq and tm'''

    seqRC= reverseComp(seq)
    xf = [seq[:a] for a in xrange(nmin, nmax)]
    xr = [seqRC[:a] for a in xrange(nmin, nmax)]

    xTMs = [calculateTm(x) for x in xf]
    xrTMs = [calculateTm(x) for x in xr]

    closestT_f = closest(xTMs, TM)
    Fprimer = xf[closestT_f]
    tmF = xTMs[closestT_f]
    closestT = closest(xrTMs, tmF)
    Rprimer = xr[closestT]
    tmR = xrTMs[closestT]

    return Fprimer, Rprimer

def designHomologyPair(seq, dist=40):

    seqRC= reverseComp(seq)
    forHomology = seq[:dist]
    revHomology = seqRC[:dist]

    return forHomology, revHomology

def designPrimerR(seq, TM, dist, nmin=15, nmax=35):
    ''' Returns a reverse primer exactly dist nucleotides away from start'''
    new = seq[:dist]
    xf = [new[-a:] for a in xrange(nmin, nmax)]
    xTMs = [calculateTm(x) for x in xf]
    closestT_f = closest(xTMs, TM)
    Fprimer = reverseComp(xf[closestT_f])
    tmF = xTMs[closestT_f]

    return Fprimer

