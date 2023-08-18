from cobra.io import read_sbml_model
import re
import pandas as pd



# generate tabe with name (carbon source), id(carbon source), if growth is possible, name (reaction), 
def fromRow2ReactionInfo(row, model):
    '''gets a row and returns the exchange reaction information (name, id, metabolites, reaction)'''
    #### If we want to add transport reactions we could use split and search 2 times for the id

    ex = re.search("R[0-9]+", row) # will return only the first occurence
    if (not ex):
        print('In fromRow2Reaction: No reaction for the following row found')
        print(row)
        return -1
    
    startPos, endPos = ex.span()
    reactionId = row[startPos:endPos]
    reactionName = model.reactions.get_by_id(reactionId)
    metabolicReaction = reactionName.reaction
    metabolite = list(reactionName.metabolites.keys())
    reversibility = reactionName.reversibility
    print("Exchange reaction info: ")
    print(reactionId, reactionName, metabolicReaction, metabolite, reversibility)
    return metabolicReaction, metabolite, reactionId, reactionName, reversibility

def fromRow2CarbonNames(row, modelGrowth):
    # paperName: split at = and take the first element + 
    equalSplit = row.split('=')
    paperName = equalSplit[0]
    print('\npaperName: ',paperName)
    if (not modelGrowth):
        return '', paperName
    # D-Mannose = EX:R1071,D-mannose exchange,0.0,0.0,TR:R843,D-mannose transport,False,"['H+_p+1', 'D-mannose_C6H12O6', 'H+_p+1', 'D-mannose_C6H12O6']","['m342[C_ex]', 'm1044[C_ex]', 'm10[C_cy]', 'm612[C_cy]']"
    # modelName: split at = (take the second argument) + 
    # split at ',' take the second element + split at ' ' and take the first argument
    # L-Arabinose = EX:R1201,L-arabinose exchange,0.0,0.0,TR:R1200,L-arabinose transport,True,"['L-arabinose_C5H10O5', 'L-arabinose_C5H10O5']","[<Metabolite m1283[C_ex] at 0x7fb93831d5b0>, <Metabolite m308[C_cy] at 0x7fb94b5d7be0>]"

    result = re.search(',((\S)+) ', equalSplit[1]) # take the name without exchange
    modelName = result.group(1)
    print("Model Name !!! ", modelName)
    print("Results !!! ", result)
    return modelName, paperName

def prepareRow(row, modelGrowth = False, model = {}):
    ''' Return the reaction and metabolitc information of a row '''
    # 1 get the carbon source name paper + model name
    modelName, paperName = fromRow2CarbonNames(row, modelGrowth)
    if (not modelGrowth):
        return modelName, paperName, '', '', '', '', ''
    # 2. get reaction id and reaction name
    metabolicReaction, metabolite, reactionId, reactionName, reversibility = fromRow2ReactionInfo(row, model)
    # 3. get id of carbon source
    return modelName, paperName, metabolicReaction, metabolite, reactionId, reactionName, reversibility

def extractInformation(modelPath, filename):
    ''' opens file and check for usable carbon sources '''
    # read model
    model = read_sbml_model(modelPath)
    bioGrowthList = []
    modelGrowthList = []
    # features for the table (only usable carbon sources)
    # carbonSourceNameModel, reactionID, reactionName, metabolites, reaction
    carbonNames = []
    reactionIDs = []
    reactionNames = []
    metabolites = []
    reactions = []
    reversibilitys = []
    # Counter will be meaningful in the next version
    growthBioNoModelCounter = 0
    noGrowthCounter = 0
    growthbioAndModelCounter = 0
    questionMarks = 0
    rowCounter = 0
    # Open the file as f.
    with open(filename) as f:
        for row in f:
            # growth is biologically possible but not for this model (iYli21)
            if (' -/' in row):
                # print(row)
                bioGrowthList.append(1)
                modelGrowthList.append(0)
                growthBioNoModelCounter += 1
            elif ('?' in row):
                # print(row)
                questionMarks += 1
                if ('|' in row):
                    print('Add exchange reaction for this metabolite')
            # growth is biologically not possible
            elif (' -' in row):
                # print(row)
                bioGrowthList.append(0)
                modelGrowthList.append(0)
                noGrowthCounter += 1
            ## growth is possible biologically and in the model
            elif ('=' in row):
                modelGrowth = True
                bioGrowthList.append(1)
                modelGrowthList.append(1)
                modelName, paperName, metabolicReaction, metabolite, reactionId, reactionName, reversibility = prepareRow(row, modelGrowth, model)
                carbonNames.append(modelName)
                reactionIDs.append(reactionId)
                reactionNames.append(reactionName)
                metabolites.append(metabolite)
                reactions.append(metabolicReaction)
                reversibilitys.append(reversibility)            
                growthbioAndModelCounter += 1
            rowCounter += 1   
    return carbonNames, reactionIDs, reactionNames, metabolites, reactions, reversibilitys

def tableForPossibleCarbonSources(modelPath, filename):
    # extract infomation from file and model
    carbonNames, reactionIDs, reactionNames, metabolites, reactions, reversibilitys = extractInformation(modelPath=modelPath,filename=filename)
    # prepare and write table
    zipped = list(zip(carbonNames,reactionIDs, reactionNames, metabolites, reactions, reversibilitys))
    df = pd.DataFrame(zipped, columns=['carbonSourceNameModel', 'reactionID', 'reactionName', 'metabolites', 'reaction', 'reversibilitys'])
    print("DF length", len(df))
    # write file
    df.to_csv('UsableCarbonSources.csv', index=False)
    print(df)
    return df

def main():
    # path to model, path to file
    modelPath = "./iYli21_v1.xml"
    filename = "carbonSourceNames_textFromImage"
    # table for all possible carbon sources (exchange reactions)
    possibleSourcesDf = tableForPossibleCarbonSources(modelPath=modelPath,filename=filename)
    
    # comming soon: table for all carbon sources and if possible to grow on (biology and model compare)
    # --> add further carbon sources you add to the model during your research
    # print('growth not in model but biologically: ', growthBioNoModelCounter)
    # print('questionmarks: ', questionMarks)
    # print('no growth: ', noGrowthCounter)
    # print('growth in model possible: ', growthbioAndModelCounter)
    # print('Number of rows: ', rowCounter)
    # print('create Main!')

if __name__=="__main__":
    main()