# Use a Hidden Markov Model to compute the most likely state sequence of a protein
    # States: Soluble (S) or Transmembrane (T)


import numpy as np
import pandas as pd


protein = 'KNSFFFFFFFLIII'
fileSeqState = 'seqState.txt'
fileSeqSoluble = 'seqSoluble.txt'
fileSeqMembrane = 'seqTransmembrane.txt'
printData = True


def loadSeq(path, concatSeqs):
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


def transitionProb(fileName):
    # Load sequences
    seq, totalAA = loadSeq(path=fileName, concatSeqs=True)
    seq = seq[0]

    # Calculate: Counts AA
    counts = {}
    for AA in seq:
        if AA in counts.keys():
            counts[AA] += 1
        else:
            counts[AA] = 1

    # Calculate: Frequency AA
    frequency = {AA: count / totalAA for AA, count in counts.items()}
    frequency = sorted(frequency.items())
    if printData:
        print('Amino Acid probability:')
        for index in range(0, len(frequency), 2):
            AA1, freq1 = frequency[index]
            if index + 1 < len(frequency):
                AA2, freq2 = frequency[index + 1]
                print(f'     {AA1}: {freq1:.6f}     {AA2:}: {freq2:.6f}')
            else:
                print(f'     {AA1}: {freq1:.6f}')
    frequency = dict(frequency)

    # Calculate: Counts digram
    countsDigram = {}
    for index in range(totalAA):
        if index == totalAA-1:
            break
        digram = seq[index:index+2]
        if digram in countsDigram.keys():
            countsDigram[digram] += 1
        else:
            countsDigram[digram] = 1
    printCounts = True
    for digram in countsDigram.keys():
        if digram[0] not in ['S', 'T'] or digram[1] not in ['S', 'T']:
            printCounts = False
            break

    # Calculate: Frequency AA
    frequencyDigram = {digram: count / totalAA for digram, count in countsDigram.items()}
    if printData:
        if fileName == fileSeqState:
            print(f'\nDigram: Frequency')
            for digram, freq in frequencyDigram.items():
                print(f'     {digram}: {freq}')

    # Calculate: Transmission probabilites
    transitionMatrix = pd.DataFrame()
    totalTransitionFreq = {}
    for digram, freq in frequencyDigram.items():
        AA1 = digram[0]
        AA2 = digram[1]
        if AA1 in totalTransitionFreq.keys():
            totalTransitionFreq[AA1] += freq
        else:
            totalTransitionFreq[AA1] = freq
    for digram, freq in frequencyDigram.items():
        AA1 = digram[0]
        AA2 = digram[1]
        transitionMatrix.loc[AA1, AA2] = freq / totalTransitionFreq[AA1]
    if printData:
        print(f'\nTransition Probabilities:\n'
              f'{transitionMatrix}\n\n')
    # sys.exit()

    return transitionMatrix, frequency


def emissionProb(filePathSoluble, filePathMembrane):
    # Load sequences
    seqSolu = loadSeq(path=filePathSoluble, concatSeqs=False)
    seqMem = loadSeq(path=filePathMembrane, concatSeqs=False)
    print(type(seqSolu))

    # Calculate: Counts AA
    countsSolu = {}
    totalAASolu = 0
    countsMem = {}
    totalAAMem = 0
    for seq in seqSolu:
        for AA in seq:
            totalAASolu += 1
            if AA in countsSolu.keys():
                countsSolu[AA] += 1
            else:
                countsSolu[AA] = 1
    for seq in seqMem:
        for AA in seq:
            totalAAMem += 1
            if AA in countsMem.keys():
                countsMem[AA] += 1
            else:
                countsMem[AA] = 1

    # Calculate: Frequency AA
    frequencySolu = {AA: count / totalAASolu for AA, count in countsSolu.items()}
    frequencySolu = sorted(frequencySolu.items())
    frequencyMem = {AA: count / totalAAMem for AA, count in countsMem.items()}
    frequencyMem = sorted(frequencyMem.items())
    if printData:
        print('Amino Acid Frequency: Soluble proteins')
        for index in range(0, len(frequencySolu), 2):
            AA1, freq1 = frequencySolu[index]
            if index + 1 < len(frequencySolu):
                AA2, freq2 = frequencySolu[index + 1]
                print(f'     {AA1}: {freq1:.6f}     {AA2:}: {freq2:.6f}')
            else:
                print(f'     {AA1}: {freq1:.6f}')
        print('\nAmino Acid Frequency: Transmembrane proteins')
        for index in range(0, len(frequencyMem), 2):
            AA1, freq1 = frequencyMem[index]
            if index + 1 < len(frequencyMem):
                AA2, freq2 = frequencyMem[index + 1]
                print(f'     {AA1}: {freq1:.6f}     {AA2:}: {freq2:.6f}')
            else:
                print(f'     {AA1}: {freq1:.6f}')

    frequencySolu = dict(frequencySolu)
    frequencyMem = dict(frequencyMem)

    aminoAcids = set(frequencySolu.keys()).union(set(frequencyMem.keys()))
    emissionProb = {AA: {
        "S": frequencySolu.get(AA, 0), # P(AA | S)
        "T": frequencyMem.get(AA, 0) # P(AA | T)
    } for AA in aminoAcids}

    # Convert to df
    emissionProbDF = pd.DataFrame(emissionProb).T.sort_index()
    if printData:
        print(f'\nEmission Probability Matrix:\n'
              f'{emissionProbDF.map(lambda x: f'{x:.3f}')}\n\n'
              f'Sum: {sum(emissionProbDF.iloc[:,0])}'
              f'    {sum(emissionProbDF.iloc[:,1])}\n\n')

    return emissionProbDF


def viterbi(proteins, initialStateProb, transitionMatrix, emissionMatrix, initialRun):
    states = ['S', 'T']

    # if printData
    print(f'Viterbi Algorithm:\n'
          f'     States: {states}\n\n'
          f'Transition Probabilities:\n{transitionMatrix}\n\n'
          f'Emission Probabilities:\n{emissionMatrix}\n')

    seqStates = []
    for seq in proteins:
        seqLen = len(seq)

        # viterbiMatrix stores the highest probability at position i with each state
        viterbiMatrix = pd.DataFrame(0.0,
                                     index=range(seqLen),
                                     columns=states)
        backpointer = np.zeros((seqLen, len(states)),
                               dtype=int) # Backpointer for backtracking

        # Initial probabilities for the first amino acid
        for indexState, state in enumerate(states):
            viterbiMatrix.iloc[0, indexState] = (emissionMatrix[state].get(seq[0], 0) *
                                initialStateProb[state])

        # Fill in the path matrix (viterbiMatrix)
        for index in range(1, seqLen):
            for indexState, state in enumerate(states):
                probMax = -1  # Maximum probability
                stateMax = -1  # State with the maximum probability
                for indexStatePrev, statePrev in enumerate(states):
                    prob = (viterbiMatrix.iloc[index - 1, indexStatePrev] *
                            transitionMatrix[statePrev][state] *
                            emissionMatrix[state].get(seq[index], 0))
                    if prob > probMax:
                        probMax = prob
                        stateMax = indexStatePrev

                viterbiMatrix.iloc[index, indexState] = probMax
                backpointer[index, indexState] = stateMax
        if initialRun:
            print(f'Viterbi Matrix:\n{viterbiMatrix}\n')

        # Find the most likely sequence of states
        sequence = ''
        statePrev = np.argmax(viterbiMatrix.iloc[seqLen - 1])
        sequence += states[statePrev]

        # Backtrack through the Viterbi matrix
        for index in range(seqLen - 2, -1, -1):
            statePrev = backpointer[index + 1, statePrev]
            sequence += states[statePrev]
        sequence = sequence[::-1] # Reverse the sequence since we backtracked
        seqStates.append(sequence)

    if initialRun:
        for index in range(len(proteins)):
            print(f'Sequence: {proteins[index]}\n'
                  f'State Sequence: {seqStates[index]}\n')
        print('')
    else:
        print('')

    return seqStates


def evaluateTContent(stateSeqs, fileName):
    print(f'Evaluating the prediceted transmembrane content in: {fileName}')
    totalProteins = len(stateSeqs)
    countsTSegments = {'0': 0, '1': 0}
    for seq in stateSeqs:
        if 'T' not in seq:
            countsTSegments['0'] += 1
        elif 'S' not in seq:
            countsTSegments['1'] += 1
        else:
            count = 0
            for index in range(len(seq)):
                if index == 0:
                    if seq[index] == 'T':
                        count += 1
                elif index == len(seq) - 1:
                    break  # Terminate if at the last index
                else:
                    digram = seq[index:index + 2]
                    if digram == 'ST':
                        count += 1
            countStr = str(count)
            if countStr not in countsTSegments.keys():
                countsTSegments[countStr] = 1
            else:
                countsTSegments[countStr] += 1
    print(f'Total proteins: {totalProteins:,}')
    for numTSeg, count in countsTSegments.items():
        avg = 100 * count / totalProteins
        print(f'Total proteins with {numTSeg} transmembrane segments: {count:,}\n'
              f'     Average proteins in this set: '
              f'{count:,} / {totalProteins:,} = {round(avg, 2)} %')



# Evaluate training and experimantal data
freqTransition, probStates = transitionProb(fileName=fileSeqState)
emissionProbabilities = emissionProb(filePathSoluble=fileSeqSoluble,
                                     filePathMembrane=fileSeqMembrane)

# Predict states
viterbi(proteins=[protein],
        initialStateProb=probStates,
        transitionMatrix=freqTransition,
        emissionMatrix=emissionProbabilities,
        initialRun=True)


# Problem: 2
seqMembrane = loadSeq(path=fileSeqMembrane, concatSeqs=False)
seqStatesMembrane = viterbi(proteins=seqMembrane,
                            initialStateProb=probStates,
                            transitionMatrix=freqTransition,
                            emissionMatrix=emissionProbabilities,
                            initialRun=False)
evaluateTContent(stateSeqs=seqStatesMembrane, fileName=fileSeqMembrane)
