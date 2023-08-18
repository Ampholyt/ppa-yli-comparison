### This file automates the overview analysis of a model
# @param path to input model
# @param path to output directory
# @output: .csv files

# 0. import, read model:
import sys
import os 
import cobra
from cobra.io import read_sbml_model, write_sbml_model, validate_sbml_model
import pandas as pd

# Check arguments
# expecting: python3 ModelOverview.py <ModelPath> <OutputPath> <verbose>
# example: python3 ModelOverview.py '/Users/ampholyt/Coding/BSEP22/Data/Models/iYli21_v1.xml' '/Users/ampholyt/Coding/BSEP22/Data/ModelAnalysis/ModelOverview/' True
n = len(sys.argv)

_verbose = True
if (n != 3):
    if (n != 4):
        message = 'Wrong number of arguments. Expected 2 or 3, got {n}. \nUsage: python ModelOverview.py <ModelPath> <OutputPath> <verbose>'.format(n=n)
        raise Exception(message)
    _verbose = False

if (_verbose):
    print('This Script will create a bunch of .csv files to give an overview over iYli21. If the path to the model and an output path is given. It will also gernerate a table for all Metabolies.\nYou can turn of the prints by setting the third argument to false')


# paths
_modelPath = sys.argv[1]
_outputDir = sys.argv[2]
_modelPrefix = _modelPath.split('/')[-1].split('.')[0]
_createOutput = True
_extremlyInterested = False
_sep = ';'
if _outputDir[-1] != '/':
    _outputDir += '/'
# check if output directory exists (create if not)
if not os.path.exists(_outputDir):
    os.makedirs(_outputDir)

def checkReadModel():
    '''reads model and checks for errors and print them if verbose'''
    (model, errors) = validate_sbml_model(_modelPath)
    if _verbose:
        print('Model will be checked for errors. If any are found, they will be printed below.')
        print(errors)
    return model

def numberDf(model):
    '''Check numbers of genes, metabolites and reactions and safe them to dataframe'''
    # investigate number and print numbers
    modelName = model.name
    reacNum = len(model.reactions)
    metaboNum = len(model.metabolites)
    geneNum = len(model.genes)
    if _verbose:
        print('The {name} has {reacNum} reactions, {metaboNum} metabolites and {geneNum} genes'.format(name=modelName, reacNum=reacNum, metaboNum=metaboNum, geneNum=geneNum))
    # generate dataframe
    tupleList = [(_modelPrefix,geneNum,metaboNum,reacNum)]
    modelNumbers = pd.DataFrame(tupleList, columns=['modelName','#Genes','#Metabolites','#Reactions'])
    if _createOutput:
        if _verbose:
            print('The model numbers will be saved to a .csv file.', "Path: ", _outputDir+_modelPrefix+'_modelNumbers.csv')
        modelNumbers.to_csv(_outputDir + _modelPrefix + '_modelNumbers.csv', index=False, sep=_sep)
    return modelNumbers

def checkInterestingPathways(model, 
    interestingPathwayList= ["Alkane metabolism", "Amino acid metabolism", "Citrate cycle (TCA cycle)", "Fatty acid metabolism", "Glycerolipid metabolism", "Glycolysis/Gluconeogenesis", "Glyoxylate and dicarboxylate metabolism", "Methane metabolism", "Pentose phosphate pathway", "Transport", "Butanoate metabolism"], 
    interestingPathwaySubset = ["Amino acid metabolism", "Citrate cycle (TCA cycle)", "Exchange", "Fatty acid metabolism", "Glycerolipid metabolism", "Glycolysis/Gluconeogenesis", "Pentose phosphate pathway", "Transport"]):
    '''Check if interesting pathways are found'''
    counter = 0
    for member in model.groups:
        for pathway in interestingPathwaySubset:
            if(member.name == pathway):
                counter += 1

    if (counter == len(interestingPathwaySubset)):
        print(f"All {counter} pathways have been found")

def pathwayExchangeDefault(model, interestingPathwayList):
    '''Creates 3 dataframes: pathway, exchange reactions and default flux values.'''
    # build dict and than dataframe
    pathwayDict = {}
    for member in model.groups:
        for pathway in interestingPathwayList:
            if(member.name == pathway):
                # get the reactions and store them in the dict
                reactionList = list(member.members)
                pathwayDict[pathway] = reactionList

    # build pathway dataframe
    pathways = []
    reactionIds = []
    reactionNames = []
    for pathwayName, reactionList in pathwayDict.items():
        for reaction in reactionList:
            pathways.append(pathwayName.lower())
            reactionIds.append(reaction.id)
            reactionNames.append(reaction.name)

    zipped = zip(pathways, reactionIds, reactionNames)
    pathwayDf = pd.DataFrame(zipped, columns=["pathwayName", "reactionId", "reactionName"])

    # 3a. Generate Exchange reactions table (get all reactions with pathwayName "exchange")
    exchangeReactions = pathwayDf[pathwayDf['pathwayName'].str.contains("exchange")]

    # 4. Generate DefaultExchangeReactionBounds table 
    # Iterate over the reactions and save bounds
    reacBounds = [] 
    for i in range(len(exchangeReactions)):
        reaction1 = exchangeReactions['reactionId'].loc[exchangeReactions.index[i]]
        reac1 = model.reactions.get_by_id(reaction1)
        reacBounds.append(reac1.bounds)

    exchangeReactions["bounds"] = reacBounds

    if _createOutput:
        if _verbose:
            print('The default exchange reaction bounds will be saved to a .csv file.', "\nPath: ", _outputDir+_modelPrefix+'_exchangeReactionBounds.csv')
            print('The pathways dataframe will be saved to a .csv file.', "\nPath: ", _outputDir+_modelPrefix+'_pathways.csv')
        pathwayDf.to_csv(_outputDir + _modelPrefix + '_pathwayTable.csv', sep=_sep, index=False)
        exchangeReactions.to_csv(_outputDir + _modelPrefix + '_exchangeReactionBounds.csv', sep=_sep, index=False)

    return exchangeReactions, pathwayDf
# useful functions
# Interacting with metabolite, gene or reaction lists
def getIds(elementList):
    '''Returns the ids of an element list (list of genes, reactions or metabolites'''
    return [elem.id for elem in elementList]

def getNames(elementList):
    '''Returns the names of the elemts in the given list (list of genes, reactions or metabolites'''
    return [elem.name for elem in elementList]

# Introduce metabolite names into reaction representation
def formulaWithNames(reaction):
    '''Returns the reaction formula with metabolite names instead of ids'''
    formula = reaction.reaction
    for metabolite in reaction.metabolites:
        formula = formula.replace(metabolite.id, metabolite.name)
    return formula

def overviewDF(model, pathwayDf):
    '''Returns an overview dataframe with the following columns: reactionId, reactionName, reactionFormula, reactionBounds, reactionObjective, reactionGenes, reactionPathways'''
    allReactions = model.reactions
    reactionIds = []
    reactionStoichimetries = []
    metaboliteIds = []
    geneIds = []
    reactionNames = []
    metabolieNames = []
    geneNames = []
    reactionWithNames = []

    # following code will use the function getIds() it is defined as usefunction at the beginning of this section
    for reaction in allReactions:
        reactionIds.append(reaction.id)
        reactionStoichimetries.append(reaction.reaction)
        reactionWithNames.append(formulaWithNames(reaction))
        metaboliteIds.append(getIds(reaction.metabolites))
        geneIds.append(getIds(reaction.genes))
        reactionNames.append(reaction.name)
        metabolieNames.append(getNames(reaction.metabolites))
        geneNames.append(getNames(reaction.genes))

    # gathering information into a dataframe    
    zipped = zip(reactionWithNames, reactionIds, reactionNames, metaboliteIds, metabolieNames, geneIds, geneNames, reactionStoichimetries)
    overviewDataframe = pd.DataFrame(zipped, columns=["reactionStoichiometry(names)", "reactionId", "reactionName", "metaboliteIds", "metaboliteNames", "geneIds", "geneNames", "reactionStoichimetry(ids)"])
    # merge pathwayDF and overview dataframe
    mergedDf = pd.merge(pathwayDf[['pathwayName', 'reactionId']], overviewDataframe, on="reactionId")
    if _createOutput:
        if _verbose:
            print('The overview dataframe will be saved to a .csv file.', "\nPath: ", _outputDir+_modelPrefix+'_overviewReactionTable.csv')
        mergedDf.to_csv(_outputDir + _modelPrefix + '_overviewReactionTable.csv', sep=_sep, index=False)
    return mergedDf

def uniqueReactions(model):
    '''Checks if there are reactions in multiple groups'''
    # iterate over each group and check if reactions are unique for each group/pathway 
    groupDict = {}
    for member in model.groups:
        groupDict[member.name] = list(member.members)

    allReactionsInUniqueGroups = True
    for groupName, reactionList in groupDict.items():
        for reaction in reactionList:
            for ggroupName, rreactionList in groupDict.items():
                if (ggroupName == groupName):
                    continue
                for rreaction in rreactionList:
                    if (reaction == rreaction):
                        allReactionsInUniqueGroups = False
                        if _extremlyInterested:
                            print("Found a reaction in 2 groups:")
                            print(f"Reaction from group {groupName} was also found in {ggroupName}")
                            print(f"The reaction found is called {reaction} and is equal to {rreaction}")
    
    if (allReactionsInUniqueGroups):
        print("All reactions are in unique groups")
    else:
        print("Not all reactions are in unique groups")
    return allReactionsInUniqueGroups

def geneReacPath(model, pathwayDf):
    '''Gene dataframe stating which reactions are associated with which genes + Pathway information'''
    # 7. Get genes and reactions per Group (At the end of iYli21_playground toy)
    # generate dataframe of genes and reactions pathways (gene -> reaction -> pathway)
    geneReacData = [(gene, reaction.name, reaction.id) for gene in model.genes for reaction in gene.reactions]
    geneReactionPathwayDf = pd.DataFrame(geneReacData, columns=["gene", "reactionName", "reactionId"])
    # merge with pathwayDf
    geneReactionPathwayDf = pd.merge(geneReactionPathwayDf, pathwayDf.drop('reactionName', axis=1), on="reactionId")
    if _createOutput:
        if _verbose:
            print('The geneReactionPathway dataframe will be saved to a .csv file.', "\nPath: ", _outputDir+_modelPrefix+'_geneReactionPathwayDf.csv')
        geneReactionPathwayDf.to_csv(_outputDir + _modelPrefix + '_geneReactionPathwayTable.csv', sep=_sep, index=False)
    return geneReactionPathwayDf

def main():
    # From Guo et al 2022 Fig. 4
    _interestingPathwayList = ["Alkane metabolism", "Amino acid metabolism", "Citrate cycle (TCA cycle)", "Fatty acid metabolism", "Glycerolipid metabolism", "Glycolysis/Gluconeogenesis", "Glyoxylate and dicarboxylate metabolism", "Methane metabolism", "Pentose phosphate pathway", "Transport", "Butanoate metabolism"]
    # interesting subset
    _interestingPathwaySubset = ["Amino acid metabolism", "Citrate cycle (TCA cycle)", "Exchange", "Fatty acid metabolism", "Glycerolipid metabolism", "Glycolysis/Gluconeogenesis", "Pentose phosphate pathway", "Transport"]
 
    # check and load model
    model = checkReadModel()
    # 1. Genes Metabolites and reactions (Numbers) + Dataframe
    numbersDf = numberDf(model)
    # 2. Check if interesting pathways are found
    checkInterestingPathways(model, _interestingPathwayList, _interestingPathwaySubset)
    # 3. Build 3 dataframes: exchange bounds, exchange reactions and pathway dataframe
    exchangeReactions, pathwayDf = pathwayExchangeDefault(model, _interestingPathwaySubset)
    # 5. Gene, Metabolite, Reaction dataframe
    overviewDf = overviewDF(model, pathwayDf)
    # 6. Check if reactions are only in one Pathway/Group
    uniqueReactions(model)
    # 7. For each gene which reactions are regulated and in which pathways they are (At the end of iYli21_playground toy)
    geneReacPathDf = geneReacPath(model, pathwayDf)

if __name__ == "__main__":
    main()