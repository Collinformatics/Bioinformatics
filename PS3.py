from loadFile import loadFasta
import pandas as pd
import sys


# Specify the file path
pathProteins = 'ProteinsHuman.fasta.txt'
printSeqs = True


def frequency(sequences):
    aminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    printData = True
    counts = pd.DataFrame(0, index=aminoAcids, columns=['Counts'])


    for ID, sequence in sequences.items():
        seqLength = len(sequence)

        # Count the AAs
        for AA in sequence:
            if AA in aminoAcids:
                counts.loc[AA, 'Counts'] += 1

        # Calculate relative frequency of the AAs
        rf = pd.DataFrame(0.0, index=aminoAcids, columns=['RF'])
        for AA in counts.index:
            rf.loc[AA, 'RF'] = counts.loc[AA, 'Counts'] / seqLength
        rf = rf['RF'].tolist()

        # Add RF list to the dictionary
        sequences[ID] = {
            'sequence': sequence,
            'RF': rf
        }
        if printData:
            print(f'\n{ID}\n{sequences[ID]}\n\n')
            printData = False

    return sequences



# Load sequences
data = loadFasta(path=pathProteins, printData=printSeqs)

# Evaluate proteins
data = frequency(sequences=data)
