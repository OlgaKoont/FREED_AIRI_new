from functools import partial
import pandas as pd
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
    
    def __call__(self, input):
        if self.preprocess:
            input = self.preprocess(input)
        property = self.property(input)
        reward = self.weight * self.reward(property)
        return reward, property

class Reward_DB:
    def __init__(self, property_df, reward, weight=1.0, preprocess=None):
        self.rewards_dict = {}  
        self.property = property_df
        self.reward = reward
        self.weight = weight
        self.preprocess = preprocess

    def add_reward(self, name, reward):
        self.rewards_dict[name] = reward
    
    def get_reward(self, name):
        return self.rewards_dict.get(name)
    
    def __call__(self, input):
       property_values = self.property(input)
       if isinstance(property_values, pd.DataFrame) or isinstance(property_values, pd.Series):
           list_of_rewards_and_properties = []
           for i in range(len(self.reward(property_values))):
               list_of_rewards_and_properties.append((self.weight * self.reward(property_values)[i], 
                                                 property_values.iloc[i]['affinity'] if isinstance(property_values,pd.DataFrame) else None))
       
       else:
           raise ValueError("False datatype")
       print("list_rewards", list_of_rewards_and_properties)
       return list_of_rewards_and_properties

'''
    def __init__(self, df_properties: pd.DataFrame,
                 reward_funcs: list,
                 weights=None,
                 preprocess=None):
        
        rewards_list = []
        
        for prop_func, rew_func in zip(property[:len(reward)], reward[:len(reward)]):
            rewards_list.append(Reward(prop_func=prop_func,
                                       reward_func=rew_func,
                                       weight=weights[len(rewards_list)],
                                       preprocess=preprocess))
        return rewards_list


   class Reward_DB:
    def __init__(self, property_list, reward, weight=1.0, preprocess=None):
        self.property_list = property_list
        self.reward = reward
        self.weight = weight
        self.preprocess = preprocess
    
    def __call__(self, input):
        if self.preprocess:
            input = self.preprocess(input)
        property_list = self.property(input)
        reward = self.weight * self.reward(property)
      
        return reward, property
    def __call__(self, input):
       property_values = self.property(input)
       if isinstance(property_values, pd.DataFrame) or isinstance(property_values, pd.Series):
           list_of_rewards_and_properties = []
           for i in range(len(self.reward(property_values))):
               list_of_rewards_and_properties.append((self.weight * self.reward(property_values)[i], 
                                                 property_values.iloc[i]['affinity'] if isinstance(property_values,pd.DataFrame) else None))
       
       else:
           raise ValueError("False datatype")
       print("list_rewards", list_of_rewards_and_properties)
       return list_of_rewards_and_properties
'''

    #def __call__(self, input):
    #    if self.preprocess:
    #       input = self.preprocess(input)
    #    property = self.property(input)
    #    if isinstance(property,  pd.DataFrame):
    #        # == True:
    #        for i in range(len(property)):
    #            reward = self.weight * self.reward(property['affinity'][i]) 
    #            print(reward, property['affinity'][i])
    #            return reward, property['affinity'][i]
    #    else:
    #        reward = self.weight * self.reward(property)
    #        print(reward, property)
    #        return reward, property


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


