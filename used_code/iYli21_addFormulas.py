
# Get iYli21 pass memots (start by adding formula)
# Adding 1261 metabolite formulas

# Imports
import cobra
from cobra.io import read_sbml_model, write_sbml_model # reading models
import re # regex
import pandas as pd # dataframes
import sys

# Check arguments
# expecting <ModelPath> <OutputPath>
n = len(sys.argv)
verbose = True
if (n != 3):
    if (n != 4):
        message = 'Wrong number of arguments. Expecting <ModelPath> <OutputPath> [verbose]'
        raise Exception(message)
    verbose = False

if (verbose):
    print('This Programm will add metabolites to the iYli21 model if the path to the model and an output path is given. And gernerate a table for all Metabolies')

# Setup:
basemodel = sys.argv[1]
outdir = sys.argv[2]
outModelName = ''.join(basemodel.split('/')[-1].split('.')[0]) + "-formulas.xml"

# Loading
model = read_sbml_model(basemodel)

# We start by generating a table of all metabolites. We will use their name and their formula. The latter one is expected to be empty for every metabolite.

# Table of metabolies (get overview)
allMetabolites = model.metabolites
metaNames = []
metaFormulas = []
metaIds = []

for metabolite in allMetabolites:
    metaNames.append(metabolite.name)
    metaFormulas.append(metabolite.formula)
    metaIds.append(metabolite.id)

zipped = zip(metaIds, metaNames, metaFormulas)
metaDf = pd.DataFrame(zipped, columns=["metaboliteId", "metaboliteName", "metaboliteFormula"])
metaDf.head()
metaDf.to_csv(outdir + '/AllMetabolites.csv')

# Search table of metabolites with regex (I looked in the created table and build a regex)
allMetabolites = model.metabolites
newFormulasCounter = 0
for metabolite in allMetabolites:
    findFormula = re.search('_(?P<formula>([A-Z](\d{1,3})?)+$)', metabolite.name)
    if findFormula:
        newFormulasCounter+=1
        metabolite.formula = findFormula.group('formula')

if (verbose):
    print(f"You added {newFormulasCounter} Formulas. While skimming over {len(allMetabolites)} metabolites. This means for {len(allMetabolites) - newFormulasCounter} metabolites we do not have a formula.")

# Table of metabolies (with formula)
allMetabolites = model.metabolites
metaNames = []
metaFormulas = []
metaIds = []

for metabolite in allMetabolites:
    metaNames.append(metabolite.name)
    metaFormulas.append(metabolite.formula)
    metaIds.append(metabolite.id)

zipped = zip(metaIds, metaNames, metaFormulas)
metaDf = pd.DataFrame(zipped, columns=["metaboliteId", "metaboliteName", "metaboliteFormula"])
metaDf.head()
metaDf.to_csv(outdir + '/FormulaAllMetabolites.csv')

# write new model:
write_sbml_model(model, outdir + outModelName)

# reading the new model
modelFormulas = read_sbml_model(outdir+outModelName) 

# skimming over the metabolites
formMetabolites = modelFormulas.metabolites
formCounter = 0
for metabolite in formMetabolites:
    if (metabolite.formula != ""):
        formCounter += 1

if (verbose):
    print(f"We skimmed over {formCounter} formulas.")
