# License: GNU Affero General Public License v3 or later

from numpy import log2

class RiPP:

    def __init__(self, sequence):
        self.sequence = sequence
        self.cys30    = ''
        self.cys20    = ''
        self.cys_ser30 = ''
        self.cys_ser20 = ''
        self.aafreq   = dict([(i,0.) for i in \
                ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", \
                 "M", "F", "P", "S", "T", "W", "Y", "V"]])
        self.clfreq   = {'RHK':0., 'DE':0., 'STNQ':0., 'CGP':0., 'AVIL':0., 'MFYW':0.}
        self.charge   = ''
        self.avgcharge= ''
        self.avghydrop= ''
        self.length   = len(self.sequence)
        self.entropy = ''
        self.entropyratio = ''

    def calculate_features(self):

        '''
        Wimley-White whole residue hydrophobicity interface scale
        '''
        hydrophobicity_dict = {'A':0.17,
        'R':0.81,
        'N':0.42,
        'D':1.23,
        'C':-0.24,
        'Q':0.58,
        'E':2.02,
        'G':0.01,
        'H':0.96,
        'I':-0.31,
        'L':-0.56,
        'K':0.99,
        'M':-0.23,
        'F':-1.13,
        'P':0.45,
        'S':0.13,
        'T':0.14,
        'W':-1.85,
        'Y':-0.94,
        'V':0.07,
        'X':0.}


        charge_dict = {'A':0.,
        'R':1.,
        'N':0.,
        'D':-1.,
        'C':0.,
        'Q':0.,
        'E':-1.,
        'G':0.,
        'H':0.5,
        'I':0.,
        'L':0.,
        'K':1.,
        'M':0.,
        'F':0.,
        'P':0.,
        'S':0.,
        'T':0.,
        'W':0.,
        'Y':0.,
        'V':0.,
        'X':0.}

        
        # --- cys20/30; includes normalization in case a peptide is shorter than 30 aa
        if self.length < 20:
            self.cys30 = self.sequence.count('C') / float(self.length)
            self.cys20 = self.sequence.count('C') / float(self.length)
            self.cys_ser30 = (self.sequence.count('C') + self.sequence.count('S')) / float(self.length)
            self.cys_ser20 = (self.sequence.count('C') + self.sequence.count('S')) / float(self.length)

        elif 20<= self.length <30:
            self.cys30 = self.sequence.count('C') / float(self.length)
            self.cys_ser30 = (self.sequence.count('C') + self.sequence.count('S')) / float(self.length)
            self.cys20 = max([self.sequence[rng:rng+20].count('C') / float(len(self.sequence[rng:rng+20])) \
                      for rng in xrange(0,self.length,1) if len(self.sequence[rng:rng+20])==20])
            self.cys_ser20 = max([(self.sequence[rng:rng+20].count('C') + self.sequence[rng:rng+20].count('S')) \
                      / float(len(self.sequence[rng:rng+20])) for rng in xrange(0,self.length,1) if \
                      len(self.sequence[rng:rng+20])==20])
            
        else:
            self.cys30 = max([self.sequence[rng:rng+30].count('C') / float(len(self.sequence[rng:rng+30])) \
                      for rng in xrange(0,self.length,1) if len(self.sequence[rng:rng+30])==30])
            self.cys20 = max([self.sequence[rng:rng+20].count('C') / float(len(self.sequence[rng:rng+20])) \
                      for rng in xrange(0,self.length,1) if len(self.sequence[rng:rng+20])==20])
            self.cys_ser30 = max([(self.sequence[rng:rng+30].count('C') + self.sequence[rng:rng+30].count('S')) \
                      / float(len(self.sequence[rng:rng+30])) for rng in xrange(0,self.length,1) if \
                      len(self.sequence[rng:rng+30])==30])
            self.cys_ser20 = max([(self.sequence[rng:rng+20].count('C') + self.sequence[rng:rng+20].count('S')) \
                      / float(len(self.sequence[rng:rng+20])) for rng in xrange(0,self.length,1) if \
                      len(self.sequence[rng:rng+20])==20])

        # --- aafreq
        for aa in self.aafreq.keys():
            self.aafreq[aa] = self.sequence.count(aa) / float(self.length)

        # --- clfreq
        for aas in self.clfreq.keys():
            self.clfreq[aas] = sum( [self.sequence.count(aa) for aa in aas] ) / float(self.length)

        # --- charge & avgcharge
        self.charge = sum([charge_dict[aa] for aa in self.sequence if aa in self.aafreq])
        self.avgcharge = self.charge / float(self.length)

        # --- avghydrop
        self.avghydrop = sum([hydrophobicity_dict[aa] for aa in self.sequence if aa in self.aafreq]) / float(self.length)

        # --- k-tuplet entropy
        Es = []
        for rng in xrange(self.length):
            s = 0.
            seqtmp = self.sequence[rng:rng+10]
            if 1:#len(seqtmp)==10:
                for i in xrange(1,len(seqtmp)):
                    try: s -= self.aafreq[seqtmp[i]]*self.aafreq[seqtmp[i]]*self.aafreq[seqtmp[i-1]] \
                         *log2(self.aafreq[seqtmp[i]]*self.aafreq[seqtmp[i-1]])
                    except KeyError: pass
            Es.append(s)
        self.entropy = max(Es)


        # --- entropy ratio
        Es = []
        for rng in xrange(self.length):
            s = 0.
            seqtmp = self.sequence[rng:rng+10]
            if len(seqtmp)<5: continue
            aatmp = set(seqtmp)
            stotal = -sum([self.aafreq[i]*log2(self.aafreq[i]) for i in aatmp if i in self.aafreq])
            stemp  = -sum([(seqtmp.count(i)/10.)*log2((seqtmp.count(i)/10.)) for i in aatmp if i in self.aafreq])
            try: Es.append( stemp/stotal )
            except ZeroDivisionError: Es.append(2.)
        try:
            self.entropyratio = min(Es)
        except:
            raise ValueError('Error with min(Es) on smORF with sequence %s' %self.sequence)

    def make_list(self):
        
        aalist = ["A", "R", "N", "D", "C", "E", "Q", "G", "H", "I", "L", "K", \
                 "M", "F", "P", "S", "T", "W", "Y", "V"]
        cllist = ['RHK', 'DE', 'STNQ', 'CGP', 'AVIL', 'MFYW']

        L = []
        L += [self.aafreq[aa] for aa in aalist]
        L += [self.clfreq[aa] for aa in cllist]
        L += [self.cys30, self.cys20, self.charge, self.avgcharge, self.avghydrop, self.entropy, self.entropyratio]
            
        return L
    
    def get_features(self):
        feature_list = [self.sequence, self.cys30, self.cys20, self.cys_ser30, self.cys_ser20, self.charge, self.avgcharge, self.avghydrop, self.length, self.entropy, self.entropyratio]
        amino_acids = self.aafreq.keys()
        clusters = self.clfreq.keys()
        clusters.sort()
        amino_acids.sort()
        feature_list += [self.aafreq[aa] for aa in amino_acids]
        feature_list += [self.clfreq[cl] for cl in clusters]
        return feature_list


'''
ripp = RiPP('METEKYLQVVEDEEIEQLVGGAGPGWVETLTKDCPWNVPVACVTIMGQRICKKCY')
#ripp = RiPP('MNKDIDLSAIEISDLISETEQSDDALSQVMAASCTTTGCACSSSSSST')
ripp.calculate_features()

print ripp.sequence, ripp.cys30, ripp.cys20, ripp.entropy
print ripp.aafreq
print ripp.clfreq
print ripp.charge, ripp.avgcharge, ripp.avghydrop, ripp.length
'''









