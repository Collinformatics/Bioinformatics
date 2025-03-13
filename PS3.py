from loadFile import loadFasta
import pandas as pd
import sys


# Specify the file path
pathProteins = 'ProteinsHuman.fasta.txt'
printSeqs = True


def frequency(sequences, saveFreq, savePath):
    aminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    printData = True
    counts = pd.DataFrame(0, index=aminoAcids, columns=['Counts'])
    vectorsRF = []


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
        rf = rf['RF'].tolist() # Convert to a list
        vectorsRF.append(f'{ID} '.join(map(str, rf)))


        # Add RF list to the dictionary
        sequences[ID] = {
            'sequence': sequence,
            'RF': rf
        }
        if printData:
            print(f'Processed Data:\n'
                  f'     {ID}\n'
                  f'     {sequences[ID]}\n')
            printData = False

    if saveFreq:
        # Save vectorsRF to a text file
        with open(savePath, 'w') as file:
            for vector in vectorsRF:
                file.write(vector + '\n')  # Write each vector on a new line

        print(f"AA Frequency data saved at:\n"
            f"     {savePath}\n")

    return sequences


# Load sequences
data = loadFasta(path=pathProteins, printData=printSeqs)

# Evaluate proteins
saveData = True
data = frequency(sequences=data, saveFreq=saveData,
                 savePath= pathProteins.replace('.fasta.txt', 'Freq.txt'))
