import numpy as np
import os
import sys
import pickle as pk
import pandas as pd
import threading
import time
from functions import filePaths, NGS

import matplotlib.pyplot as plt
import networkx as nx
import pandas as pd
import sys
from Trie import Trie



inNumberOfMotifs = 50
inOffset = 2000
inNodeSize = 300
inFontSize = 10
inScaleX = 2
inScaleY = 1



# ===================================== User Inputs ======================================
# Input 1: Select Dataset
inEnzymeName = 'Mpro2'
inBasePath = f'/Users/ca34522/Documents/Research/NGS/{inEnzymeName}'
inFilePath = os.path.join(inBasePath, 'Extracted Data')
inSavePathFigures = os.path.join(inBasePath, 'Figures')
inFileNamesInitial, inFileNamesFinal, inAAPositions = filePaths(enzyme=inEnzymeName)
inSaveFigures = True

# Input 2: Processing The Data
inPlotentropy = False
inPlotEnrichmentMap = True
inPlotEnrichmentMotif = True
inPlotWeblogoMotif = False
inPlotWordCloud = False
inPlotAADistribution = False
inPlotPositionalRF = False # For understanding shannon entropy
inPrintLoadedSubs = True
inPrintSampleSize = True
inPrintCounts = True
inPlotCounts = False
inCountsColorMap = ['white','white','lightcoral','red','firebrick','darkred']
inStDevColorMap = ['white','white','#FF76FA','#FF50F9','#FF00F2','#CA00DF','#BD16FF']
inPrintRF = True
inPrintES = True
inPrintEntopy = True
inPrintMotifData = True
inPrintNumber = 10
inCodonSequence = 'NNS' # Base probabilities of degenerate codons (can be N, S, or K)
inUseCodonProb = False # If True: use "inCodonSequence" for baseline probabilities
                       # If False: use "inFileNamesInitial" for baseline probabilities

# Input 3: Computational Parameters
inFilterSubstrates = True
inFixedResidue = ['Q']
inFixedPosition = [4]
inExcludeResidues = False # Do you want to remove any AAs from your collection of substrate
inExcludedResidue = ['A','A']
inExcludedPosition = [9,10]
inMinimumSubstrateCount = 10
inPrintFixedSubs = True
inMinDeltaS = 0.6
inManualEntropy = False # Manually select important residues in the substrate
inManualFrame = ['R4', 'R5', 'R6', 'R2']
inFixEntireSubstrateFrame = False
inExcludedResidue = ['']
inExcludedPosition = []

# Input 4: PCA
inRunPCA = False
inBinSubsPCA = False
inIndexNTerminus = 2 # Define bounds for binned substrate
inBinnedSubstrateLength = 5 # Define the length of you substrate
inFramePositions = [inIndexNTerminus - 1,
                    inIndexNTerminus + inBinnedSubstrateLength - 1]
inAAPositionsBinned = inAAPositions[inFramePositions[0]:inFramePositions[-1]]

# Input 7: Plot Heatmap
inTitleEnrichmentMap = inEnzymeName
inYLabelEnrichmentMap = 2 # 0 for full Residue name, 1 for 3-letter code, 2 for 1 letter
inPrintSelectedSubstrates = 1 # Set = 1, to print substrates with fixed residue
inFigureSize = (9.5, 8) # (width, height)
inFigureBorders = [0.882, 0.075, 0.117, 0.998]

inEnrichmentColorMap = ['navy','royalblue','dodgerblue','lightskyblue','white','white',
                        'lightcoral','red','firebrick','darkred']

# Input 8: Plot Sequence Motif
inNormLetters = True # Normalize fixed letter heights
inShowWeblogoYTicks = True
inAddHorizontalLines = False
inTitleMotif = inTitleEnrichmentMap
inBigLettersOnTop = False
inFigureSizeMotif = inFigureSize
inFigureBordersMotifYTicks = [0.882, 0.075, 0.07, 0.98] # [top, bottom, left, right]
inFigureBordersMotifMaxYTick = [0.882, 0.075, 0.102, 0.98]
inFigureBordersEnrichmentMotif = [0.882, 0.075, 0.138, 0.98]
inLetterColors = ['darkgreen','firebrick','deepskyblue','pink','navy','black','gold']
                  # Aliphatic, Acidic, Basic, Hydroxyl, Amide, Aromatic, Sulfur

# Input _____: Figure Parameters
inFigureTitleSize = 18
inFigureLabelSize = 16
inFigureTickSize = 13
inShowSampleSize = True # Include the sample size in your figures



# =================================== Setup Parameters ===================================
# Figure parameters
figSizeEM = inFigureSize
figBordersEM = inFigureBorders
colors = {
    'A': inLetterColors[0],
    'R': inLetterColors[2],
    'N': inLetterColors[4],
    'D': inLetterColors[1],
    'C': inLetterColors[6],
    'E': inLetterColors[1],
    'Q': inLetterColors[4],
    'G': inLetterColors[0],
    'H': inLetterColors[2],
    'I': inLetterColors[0],
    'L': inLetterColors[0],
    'K': inLetterColors[2],
    'M': inLetterColors[6],
    'F': inLetterColors[5],
    'P': inLetterColors[0],
    'S': inLetterColors[3],
    'T': inLetterColors[3],
    'W': inLetterColors[5],
    'Y': inLetterColors[5],
    'V': inLetterColors[0]
}

# Colors: Console
white = '\033[38;2;255;255;255m'
silver = '\033[38;2;204;204;204m'
purple = '\033[38;2;189;22;255m'
magenta = '\033[38;2;255;0;128m'
pink = '\033[38;2;255;0;242m'
cyan = '\033[38;2;22;255;212m'
green = '\033[38;2;5;232;49m'
greenLight = '\033[38;2;204;255;188m'
greenDark = '\033[38;2;30;121;13m'
yellow = '\033[38;2;255;217;24m'
orange = '\033[38;2;247;151;31m'
red = '\033[91m'
resetColor = '\033[0m'



# =================================== Initialize Class ===================================
ngs = NGS(enzymeName=inEnzymeName, substrateLength=len(inAAPositions),
          fixedAA=inFixedResidue, fixedPosition=inFixedPosition,
          minCounts=inMinimumSubstrateCount, colorsCounts=inCountsColorMap,
          colorStDev=inStDevColorMap, colorsEM=inEnrichmentColorMap,
          colorsMotif=inLetterColors, xAxisLabels=inAAPositions,
          xAxisLabelsBinned=inAAPositionsBinned, residueLabelType=inYLabelEnrichmentMap,
          titleLabelSize=inFigureTitleSize, axisLabelSize=inFigureLabelSize,
          tickLabelSize=inFigureTickSize, printNumber=inPrintNumber,
          showNValues=inShowSampleSize, saveFigures=inSaveFigures, savePath=inFilePath,
          savePathFigs=inSavePathFigures, setFigureTimer=None)



# ====================================== Load Data =======================================
startTimeLoadData = time.time()
if inFilterSubstrates:
    fixedSubSeq = ngs.fixSubstrateSequence(exclusion=inExcludeResidues,
                                           excludedAA=inExcludedResidue,
                                           excludedPosition=inExcludedPosition)
    datasetTag = f'Fixed {fixedSubSeq}'
else:
    fixedSubSeq = None
    datasetTag = None


loadUnfixedSubstrates = True
if inFilterSubstrates:
    filePathFixedCountsFinal = os.path.join(inFilePath,
                                            f'counts_{inEnzymeName}_'
                                            f'FinalSort_{fixedSubSeq}_'
                                            f'MinCounts_{inMinimumSubstrateCount}')
    filePathFixedSubsFinal = os.path.join(inFilePath,
                                          f'fixedSubs_{inEnzymeName}_'
                                          f'FinalSort_{fixedSubSeq}_'
                                          f'MinCounts_{inMinimumSubstrateCount}')


    # Verify that the files exists
    if (os.path.exists(filePathFixedSubsFinal) and
            os.path.exists(filePathFixedCountsFinal)):
        loadUnfixedSubstrates = False

        # # Load the data: Initial sort
        filePathsInitial = []
        for fileName in inFileNamesInitial:
            filePathsInitial.append(os.path.join(inFilePath, f'counts_{fileName}'))


        # Verify that all files exist
        missingFile = False
        indexMissingFile = []
        for indexFile, path in enumerate(filePathsInitial):
            if os.path.exists(path):
                continue
            else:
                missingFile = True
                indexMissingFile.append(indexFile)
        if missingFile:
            print('\033[91mERROR: File not found at path:')
            for indexMissing in indexMissingFile:
                print(f'     {filePathsInitial[indexMissing]}')
            print(
                f'\nMake sure your path is correctly named, and that you have already '
                f'extracted and counted your NGS data\n')
            sys.exit()
        else:
            countsInitial, countsInitialTotal = ngs.loadCounts(
                filePath=filePathsInitial, files=inFileNamesInitial,
                printLoadedData=inPrintCounts, fileType='Initial Sort', fixedSeq=None)
            # Calculate RF
            initialRF = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                                        fileType='Initial Sort', printRF=inPrintRF)


        # # Load the data: Final sort
        print('================================= Load: Substrate Files '
              '=================================')
        print(f'Loading file at path:\n'
              f'     {greenDark}{filePathFixedSubsFinal}{resetColor}\n\n')
        with open(filePathFixedSubsFinal, 'rb') as file:
            substratesFinal = pk.load(file)
        print(f'Loaded Substrates:{purple} {inEnzymeName} Fixed {fixedSubSeq}')
        iteration = 0
        for substrate, count in substratesFinal.items():
            print(f'     {silver}{substrate}{resetColor}, Count:{red} {count:,}'
                  f'{resetColor}')
            iteration += 1
            if iteration >= 10:
                print('\n')
                break

        # Load: Fixed Counts
        print('================================== Load: Counted Files '
              '==================================')
        print(f'Loading file at path:\n'
              f'     {greenDark}{filePathFixedCountsFinal}{resetColor}\n\n')
        countsFinal = pd.read_csv(filePathFixedCountsFinal, index_col=0)
        countsFinalTotal = sum(countsFinal.iloc[:, 0])
        print(f'Loaded Counts:{purple} {inEnzymeName} Fixed {fixedSubSeq}'
              f'\n{red}{countsFinal}{resetColor}\n\n'
              f'Total substrates:{white} {countsFinalTotal:,}{resetColor}\n\n')
    else:
        print(f'File not found:\n'
              f'     {filePathFixedSubsFinal}'
              f'\n\nLoading substrates and fixing the residue(s):'
              f'{purple} {fixedSubSeq}{resetColor}\n\n')


if loadUnfixedSubstrates:
    def loadSubstrates(filePath, fileNames, fileType, printLoadedData, result):
        subsLoaded, totalSubs = ngs.loadSubstrates(filePath=filePath,
                                                   fileNames=fileNames,
                                                   fileType=fileType,
                                                   printLoadedData=printLoadedData)
        result[fileType] = (subsLoaded, totalSubs)


    # Initialize result dictionary
    loadedResults = {}

    # Create threads for loading initial and final substrates
    threadInitial = threading.Thread(target=loadSubstrates,
                                     args=(inFilePath, inFileNamesInitial,
                                           'Initial Sort', inPrintCounts, loadedResults))
    threadFinal = threading.Thread(target=loadSubstrates,
                                   args=(inFilePath, inFileNamesFinal,
                                         'Final Sort', inPrintCounts, loadedResults))

    # Start the threads
    threadInitial.start()
    threadFinal.start()

    # Wait for the threads to complete
    threadInitial.join()
    threadFinal.join()
    time.sleep(0.5)

    # Retrieve the loaded substrates
    substratesInitial, totalSubsInitial = loadedResults['Initial Sort']
    substratesFinal, totalSubsFinal = loadedResults['Final Sort']

    # Load the data: Initial sort
    filePathsInitial = []
    for fileName in inFileNamesInitial:
        filePathsInitial.append(os.path.join(inFilePath, f'counts_{fileName}'))
    # if '/' in inFilePath:
    #     for fileName in inFileNamesInitial: # _MinCounts_{inMinimumSubstrateCount}
    #         filePathsInitial.append(f'{inFilePath}/counts_{fileName}')
    # else:
    #     for fileName in inFileNamesInitial:
    #         filePathsInitial.append(f'{inFilePath}\\counts_{fileName}')

    # Verify that all files exist
    missingFile = False
    indexMissingFile = []
    for indexFile, path in enumerate(filePathsInitial):
        if os.path.exists(path):
            continue
        else:
            missingFile = True
            indexMissingFile.append(indexFile)
    if missingFile:
        print('\033[91mERROR: File not found at path:')
        for indexMissing in indexMissingFile:
            print(f'     {filePathsInitial[indexMissing]}')
        print(f'\nMake sure your path is correctly named, and that you '
              f'have already extracted and counted your NGS data\n')
        sys.exit()
    else:
        countsInitial, countsInitialTotal = ngs.loadCounts(filePath=filePathsInitial,
                                                           files=inFileNamesInitial,
                                                           printLoadedData=inPrintCounts,
                                                           fileType='Initial Sort',
                                                           fixedSeq=None)
        # Calculate RF
        initialRF = ngs.calculateRF(counts=countsInitial, N=countsInitialTotal,
                                    fileType='Initial Sort', printRF=inPrintRF)

    # Load the data: final sort
    filePathsFinal = []
    if '/' in inFilePath:
        for fileName in inFileNamesFinal:
            filePathsFinal.append(f'{inFilePath}/counts_{fileName}')
    else:
        for fileName in inFileNamesFinal:
            filePathsFinal.append(f'{inFilePath}\\counts_{fileName}')

    # Verify that all files exist
    missingFile = False
    for indexFile, path in enumerate(filePathsFinal):
        if os.path.exists(path):
            continue
        else:
            missingFile = True
            indexMissingFile.append(indexFile)
    if missingFile:
        print('\033[91mERROR: File not found at path:')
        for indexMissing in indexMissingFile:
            print(filePathsFinal[indexMissing])
        print(f'\nMake sure your path is correctly named, and that you '
              f'have already extracted and counted your NGS data\n')
        sys.exit()
    else:
        countsFinal, countsFinalTotal = ngs.loadCounts(filePath=filePathsFinal,
                                                       files=inFileNamesFinal,
                                                       printLoadedData=inPrintCounts,
                                                       fileType='Final Sort',
                                                       fixedSeq=None)

# Calculate RF
finalRF = ngs.calculateRF(counts=countsFinal, N=countsFinalTotal,
                          fileType='Final Sort', printRF=inPrintRF)

# Time
endTime = time.time()
runTime = endTime - startTimeLoadData
if inFilterSubstrates:
    print(f'Runtime:{white} Loading {fixedSubSeq} substrates\n'
          f'     {red} {np.round(runTime, 3)}s{resetColor}\n\n')
else:
    print(f'Runtime:{white} Loading unfixed substrates\n'
          f'     {red} {np.round(runTime, 3)}s{resetColor}\n\n')



# =================================== Define Functions ===================================
def pressKey(event):
    if event.key == 'escape':
        plt.close()



def extractMotif(substrates, indexSubFrame, datasetTag):
    print('================================= Extract Motif '
          '=================================')
    print(f'Dataset: {datasetTag}')
    motifs = {}
    indexStart = min(indexSubFrame)
    indexEnd = max(indexSubFrame)

    # Print substrates
    iteration = 0
    for substrate, count in substrates.items():
        iteration += 1
        print(f'Substrate:{pink} {substrate}{resetColor}\n'
              f'     Count:{red} {count:,}{resetColor}')
        if iteration >= inPrintNumber:
            break
    print('\n')


    # Extract the motifs
    iteration = 0
    for substrate, count in substrates.items():
        motif = substrate[indexStart:indexEnd+1]
        if motif in motifs:
            motifs[motif] += count
        else:
            motifs[motif] = count
        # if iteration >= inNumberOfMotifs:
        #     break
        # iteration += 1
    motifs = dict(sorted(motifs.items(), key=lambda item: item[1], reverse=True))

    # Print motifs
    iteration = 0
    print('Extracted Motifs')
    for  motif, counts in motifs.items():
        iteration += 1
        print(f'Motif:{yellow} {motif}{resetColor}\n'
              f'     Count:{red} {counts:,}{resetColor}')
        if iteration >= inPrintNumber:
            break
    print('\n')

    return motifs



def printTrie(node, level=0, path=""):
    if level == 0:
        print("Trie Structure:")

    # Recursively print the Trie structure
    if node is None:
        return

    # Print the current node's path and level
    print("  " * level + f"Level {level}: {path}")

    # Recursively print all children of the current node
    for char, nodeChild in node.children.items():
        printTrie(nodeChild, level + 1, path + char)



def suffixTrie(subs, entropy, title):
    trie = Trie()  # Initialize Trie

    # Find motif positions based on entropy threshold
    indexPos = []
    for index in entropy.index:
        posEntropy = entropy.loc[index, 'ΔEntropy']
        if posEntropy >= inMinDeltaS:
            indexPos.append(int(index.replace('R', '')) - 1)

    # Make the trie
    motifs = set() # Use a set to avoid duplicate motifs
    countsMotif = 0
    for seq in subs.keys():
        # Extract important AAs from the motif
        motif = ''.join(seq[index] for index in indexPos)
        countsMotif += 1
        motifs.add(motif)
        trie.insert(motif)

    # Plot the Trie
    plotTrie(trie, title, countsMotif)



def addNodesToGraph(node, graph, scaleX, scaleY, offset=inOffset, nodeSize=inNodeSize):
    pos = {}  # Stores node positions
    nodeCountLevel = {}  # Track all nodes per level
    queue = [(node, None, '', 0)]  # (node, parent, char, level)

    # Pass 1: Collect nodes for each level before positioning
    while queue:
        nodeCurrent, parent, char, level = queue.pop(0)
        nodeID = f"{char}-{level}-{id(nodeCurrent)}"

        if level not in nodeCountLevel:
            nodeCountLevel[level] = []
        nodeCountLevel[level].append((nodeCurrent, parent, char, nodeID))

        # Add children to the queue
        for child_char, nodeChild in nodeCurrent.children.items():
            queue.append((nodeChild, nodeID, child_char, level + 1))

    # Pass 2: Assign positions level by level and group children under their parents
    for level, nodes in nodeCountLevel.items():
        nodeNumber = len(nodes)
        for i, (nodeCurrent, parent, char, nodeID) in enumerate(nodes):
            if parent is None:
                pos[nodeID] = (0, 0)  # Root node position
            else:
                parentX, parentY = pos[parent]

                # Calculate the x-position for the children based on their order (index)
                if nodeNumber % 2 == 1:
                    # Odd number of nodes: Center one on parentX
                    posX = parentX + (i - nodeNumber // 2) * offset
                else:
                    # Even number of nodes: Spread symmetrically
                    posX = parentX + (i - (nodeNumber / 2 - 0.5)) * offset

                # Adjust for node size (ensure spacing is enough to prevent overlap)
                # We calculate the "radius" of the node (half of its size)
                # and make sure nodes are spaced at least this far apart
                nodeRadius = nodeSize / 2
                while any(abs(posX - otherX) < nodeSize for otherID, (otherX, otherY) in pos.items()):
                    # If there's overlap, increase the horizontal distance
                    posX += nodeRadius * 2  # Increase spacing

                # Group children under the parent node
                pos[nodeID] = (posX, parentY - scaleY)

            # Print debugging info in a single statement
            print(f'Level: {level}, Total Nodes at Level: {nodeNumber}\n'
                  f'Parent: {parent}\n'
                  f'    X: {parentX if parent else "N/A"}, '
                  f'Y: {parentY if parent else "N/A"}\n'
                  f'Node: {nodeID}\n'
                  f'    Pos: {pos[nodeID]}\n')

            # Add node and edge to graph
            graph.add_node(nodeID, label=char)
            if parent is not None:
                graph.add_edge(parent, nodeID)

    return pos



def plotTrie(trie, title, countsMotif):
    printTrie(trie.root)
    graph = nx.DiGraph()
    pos = addNodesToGraph(trie.root, graph,
                          scaleX=inScaleX,
                          scaleY=inScaleY,
                          offset=inOffset)

    # Get node labels
    labels = {node: data['label'] for node, data in graph.nodes(data=True)}

    # Print: Dataset tag
    if datasetTag is None:
        print(f'\nDataset: '
              f'Suffix Tree - {inEnzymeName} - {countsMotif} - Unfixed\n\n')
    else:
        print(f'\nDataset: '
              f'Suffix Tree - {inEnzymeName} - {countsMotif} - {datasetTag}\n\n')

    # Create figure and connect key event
    fig, ax = plt.subplots(figsize=(12, 8))  # Increased figure size
    fig.canvas.mpl_connect('key_press_event', pressKey)

    # Draw graph with custom colors
    nx.draw(graph, pos, with_labels=True, labels=labels, node_size=inNodeSize,
            node_color="#F18837", font_size=inFontSize, font_weight="bold",
            edge_color="#101010", ax=ax)

    plt.title(title, fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.show()



# ==================================== Evaluate Data =====================================
# Calculate: Average initial RF
RFsum = np.sum(initialRF, axis=1)
initialRFAvg = RFsum / len(initialRF.columns)
initialRFAvg = pd.DataFrame(initialRFAvg, index=initialRFAvg.index,
                            columns=['Average RF'])

if inFilterSubstrates:
    titleSuffixTree = f'{inEnzymeName}: {fixedSubSeq}'
    datasetTag = f'{fixedSubSeq}'

    # Fix AA
    substratesFinal, counFinalTotal = ngs.fixResidue(
        substrates=substratesFinal, minimumCount=inMinimumSubstrateCount,
        exclusion=inExcludeResidues, excludedAA=inExcludedResidue,
        excludePositon=inExcludedPosition, fixedString=fixedSubSeq,
        printRankedSubs=inPrintFixedSubs, sortType='Final Sort')

    # Count fixed substrates
    countsFinal, _ = ngs.countResidues(substrates=substratesFinal,
                                       datasetType='Final Sort',
                                       printCounts=inPrintCounts)
else:
    titleSuffixTree = f'{inEnzymeName}: Fixed {fixedSubSeq}'
    datasetTag = None


# Calculate: Positional entropy
entropy = ngs.calculateEntropy(RF=finalRF, printEntropy=inPrintEntopy)
entropySubFrame = ngs.findSubstrateFrame(entropy=entropy,
                                         minEntropy=0.6,
                                         fixFullFrame=False)
entropyIndexList = list(entropy.index)
entropySubFrameIndex = [entropyIndexList.index(idx) for idx in entropySubFrame.index]

# Extract Motif
motifs = extractMotif(substrates=substratesFinal,
                      indexSubFrame=entropySubFrameIndex,
                      datasetTag=datasetTag)

# Make suffix tree
suffixTrie(subs=motifs, entropy=entropySubFrame, title=titleSuffixTree)
