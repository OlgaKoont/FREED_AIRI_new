from functools import partial
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os, sys, pickle
from rdkit.Chem import RDConfig
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


class Reward:
    def __init__(self, property, reward, weight=1.0, preprocess=None):
        self.property = property
        self.reward = reward
        self.weight = weight
        self.preprocess = preprocess

class Reward:
    def __init__(self, property, reward, weight=1.0, preprocess=None):
        self.property = property
        self.reward = reward
        self.weight = weight
        self.preprocess = preprocess

    def __call__(self, input):
        property_values = self.property(input)
        if isinstance(property_values, pd.DataFrame) or isinstance(property_values, pd.Series):
            list_of_rewards_and_properties = []
            for i in range(len(self.reward(property_values))):
                list_of_rewards_and_properties.append((self.weight * self.reward(property_values)[i], 
                                                  property_values.iloc[i]['affinity'] if isinstance(property_values,pd.DataFrame) else None))
        
        else:
            raise ValueError("False datatype")
    
        return list_of_rewards_and_properties


def identity(x):
    return x


def ReLU(x):
    return max(x, 0)


def HSF(x):
    return float(x > 0)


class OutOfRange:
    def __init__(self, lower=None, upper=None, hard=True):
        self.lower = lower
        self.upper = upper
        self.func = HSF if hard else ReLU

    def __call__(self, x):
        y, u, l, f = 0, self.upper, self.lower, self.func
        if u is not None:
            y += f(x - u)
        if l is not None:
            y += f(l - x)
        return y


class PatternFilter:
    def __init__(self, patterns):
        self.structures = list(filter(None, map(Chem.MolFromSmarts, patterns)))

    def __call__(self, molecule):
        return int(any(molecule.HasSubstructMatch(struct) for struct in self.structures))


def MolLogP(m):
    return rdMolDescriptors.CalcCrippenDescriptors(m)[0]


def BRSASCORE(m):
    smi = Chem.MolToSmiles(m)
    scorer = SAScorer()
    score, _ = scorer.calculateScore(smi)
    return score

def check_sascore(m):
  """ Calculates SAScore """
  return sascorer.calculateScore(m)


