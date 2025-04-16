import random


# This script can be used to generate a randomized set of protein sequences

# User Inputs
numSubs = 10
subLen = 8
fixAA = 'W'
fixPos = 5

# Amino Acids
AA = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I',
      'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
numAA = len(AA) - 1

# Generate random substrate
substrates = []
for countSubs in range(numSubs):
    substrate = ''
    for index in range(subLen):
        if index == fixPos - 1:
            substrate += fixAA
        else:
            indexAA = random.randint(0, numAA)
            substrate += AA[indexAA]
    substrates.append(substrate)

# Print the data
for countSubs, substrate in enumerate(substrates):
    print(countSubs, substrate)

# Save the data
with open('exampleData.txt', 'w') as fileSave:
    for substrate in substrates:
        fileSave.write(f'{substrate}\n')
