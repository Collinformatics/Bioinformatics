from Bio.Seq import Seq
import re


DNA = ('AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGG'
       'ATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG')
fileName = 'rosalind_orf.txt'


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


def findORFs(dnaSeq):
    startCodon = "ATG"
    stopCodons = {"TAA", "TAG", "TGA"}
    ORFs = []

    # Analyze both forward and reverse complement strands
    sequences = {
        "Forward": dnaSeq,
        "Reverse Complement": str(Seq(dnaSeq).reverse_complement())
    }

    for strand, seq in sequences.items():
        print(f"\nScanning {strand} Strand:")
        startPositions = [match.start() for match in re.finditer(startCodon, seq)]

        for start in startPositions:
            for i in range(start + 3, len(seq), 3):
                codon = seq[i:i + 3]
                if codon in stopCodons:
                    orf = seq[start:i]  # Extract ORF sequence
                    ORFs.append((strand, orf))
                    print(f"  ORF: {orf}")  # Print each ORF found
                    break  # Stop at the first valid stop codon

    expressDNA(ORFs)


def expressDNA(seqs):
    substrates = {}

    print(f'\nExpressed ORFs:')
    for read, seq in seqs:
        # Express substrate
        substrate = Seq.translate(seq)
        if substrate in substrates.keys():
            substrates[substrate] += 1
        else:
            substrates[substrate] = 1

    for substrate in substrates.keys():
        print(f'     {substrate}')
    print('\n')


# Evaluate the test seq
findORFs(DNA)

# Evaluate sequence in a text file
seqRosa, countBP = loadSeq(fileName, True)
seqRosa = seqRosa[0].replace('>Rosalind_6530', '')
findORFs(seqRosa)
