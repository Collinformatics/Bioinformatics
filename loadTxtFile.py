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
