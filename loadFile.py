def loadSeq(path, concatSeqs):
    printData = True
    with open(path, 'r') as file:
        if concatSeqs:
            # Load sequences as one long string
            loadedSeq = [file.read().replace('\n', '')]
            loadedSeqLen = len(loadedSeq[0])
            if printData:
                print(f'Loading File: {path}\n'
                      f'     Sequence length: {loadedSeqLen:,}\n')
            return loadedSeq
        else:
            # Load sequences as individual strings
            loadedSeq = file.read().strip().split('\n')
            if printData:
                print(f'Loading File: {path}')
            return loadedSeq


def loadFasta(path, printData):
    from Bio import SeqIO

    data = {record.id: record.seq for record in SeqIO.parse(path, "fasta")}

    if printData:
        for recordID, sequence in data.items():
            print(f'ID: {recordID} (name|id|description)\n'
                  f'Sequence: {sequence}\n')
        print(f'Loaded File:\n     {path}\n\n')
    return data
    
