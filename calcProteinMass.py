protein = 'SKADYEK'
fileName = 'rosalind_prtm.txt'
sigfigs = 3



mass = {
    'A': 71.03711,
    'C': 103.00919,
    'D': 115.02694,
    'E': 129.04259,
    'F': 147.06841,
    'G': 57.02146,
    'H': 137.05891,
    'I': 113.08406,
    'K': 128.09496,
    'L': 113.08406, 
    'M': 131.04049, 
    'N': 114.04293, 
    'P': 97.05276, 
    'Q': 128.05858, 
    'R': 156.10111,
    'S': 87.03203,
    'T': 101.04768,
    'V': 99.06841,
    'W': 186.07931,
    'Y': 163.06333
}

def calculateMass(proteinSeq, seqTag):
    MW = sum(mass[AA] for AA in proteinSeq)
    print(f"Monoisotopic Mass of {seqTag}: {MW:.{sigfigs}f} g/mol\n")

def loadSeq(path, concatSeqs):
    printData = True
    with open(path, 'r') as file:
        if concatSeqs:
            # Load sequeces as one long string
            loadedSeq = [file.read().replace('\n', '')]
            loadedSeqLen = len(loadedSeq[0])
            if printData:
                print(f'Loading File: {path}\n'
                      f'     Sequence length: {loadedSeqLen:,}\n')
            return loadedSeq, loadedSeqLen
        else:
            # Load sequeces as individual strings
            loadedSeq = file.read().strip().split('\n')
            if printData:
                print(f'Loading File: {path}')
            return loadedSeq


calculateMass(proteinSeq=protein, seqTag=protein)
seqs = loadSeq(path=fileName, concatSeqs=False)
calculateMass(proteinSeq=seqs[0], seqTag=fileName.replace('.txt', ''))
