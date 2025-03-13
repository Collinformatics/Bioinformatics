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

    proteins = SeqIO.parse(path, "fasta")
    if printData:
        # Parse the file and iterate over each record
        for record in proteins:
            print(f'ID: {record.id}\n'
                  f'Description: {record.description}\n'
                  f'Sequence: {record.seq}\n')
    return proteins
    
