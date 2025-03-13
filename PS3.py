from loadFile import loadFasta

# Specify the file path
pathProteins = 'ProteinsHuman.fasta.txt'
printSeqs = True

# Load sequences
data = loadFasta(path=pathProteins, printData=printSeqs)
