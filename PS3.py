from loadFile import loadFasta
import pandas as pd



# Specify the file path
pathProteins = 'ProteinsHuman.fasta.txt'
printSeqs = True
saveData = True
printData = True


def frequency(sequences, printData, saveFreq, savePath):
    # Evaluate a fasta file and determine the relative frequency of AAs in the sequence
    print(f'Calculating AA Frequency:')
    aminoAcids = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L',
                  'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']
    vectorsRF = []
    printNumber = 100 # Print every nth number datapoint

    # Evaluate the data
    iteration = 0
    for ID, sequence in sequences.items():
        seqLength = len(sequence)

        # Count the AAs
        counts = pd.DataFrame(0, index=aminoAcids, columns=['Counts'])
        for AA in sequence:
            if AA in aminoAcids:
                counts.loc[AA, 'Counts'] += 1

        # Calculate relative frequency of the AAs
        rf = pd.DataFrame(0.0, index=aminoAcids, columns=['RF'])
        for AA in counts.index:
            rf.loc[AA, 'RF'] = counts.loc[AA, 'Counts'] / seqLength
        rf = rf['RF'].tolist() # Convert to a list
        rfVector = f"{ID} " + " ".join(map(str, rf))
        vectorsRF.append(rfVector)

        # Add RF list to the dictionary
        sequences[ID] = {
            'sequence': sequence,
            'RF': rf
        }

        if printData:
            if iteration % printNumber == 0:
                print(f'Index {iteration:,}: {ID}\n'
                      f'{rf}\n')
            iteration += 1

    if saveFreq:
        # Save vectors to a text file
        with open(savePath, 'w') as file:
            file.write(" " + " ".join(aminoAcids) + '\n')
            for vector in vectorsRF:
                file.write(vector + '\n')

    return sequences


# Load sequences
data = loadFasta(path=pathProteins, printData=printSeqs)

# Evaluate proteins
data = frequency(sequences=data, printData=printData, saveFreq=saveData,
                 savePath=pathProteins.replace('.fasta.txt', 'FreqAA.txt'))
