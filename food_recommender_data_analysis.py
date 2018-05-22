#!/usr/bin/env python
"""Food group recommender based on the nutrient content of different food groups.

Use USDA Food Database data given in https://github.com/wesm/pydata-book. The data initialization techniques are directly based on the example codes given in the book. These parts are marked below in the code. 
"""
import json
import pandas as pd
from pandas import Series, DataFrame
import matplotlib.pyplot as plt
import sys


def initialize_nutrient_data():
    """Load USDA data and initialize nutrient data.

    This function is directly based on the data processing examples in https://github.com/wesm/pydata-book.
    """
    # load data
    db = json.load(open('datasets/usda_food/database.json'))
    # select the needed columns
    info_keys = ['description','group','id']
    info = pd.DataFrame(db, columns=info_keys)
    # read nutrients data into a list of DataFrames and add column 'id' to enable merging with info
    nutrients = []
    for rec in db:
        fnuts = pd.DataFrame(rec['nutrients'])
        fnuts['id'] = rec['id']
        nutrients.append(fnuts)

    # concatenate the list into a DataFrame object
    nutrients = pd.concat(nutrients, ignore_index=True)
    # drop duplicates
    nutrients = nutrients.drop_duplicates()
    # rename duplicate columns in info and nutrients
    col_mapping = {'description' : 'food',
                   'group' : 'fgroup'}
    info = info.rename(columns=col_mapping, copy=False)
    col_mapping = {'description' : 'nutrient',
                   'group' : 'nutgroup'}
    nutrients = nutrients.rename(columns=col_mapping, copy=False)
    # merge info and nutrients
    ndata = pd.merge(nutrients, info, on='id', how='outer')
    return ndata


def food_group_nutrient_contents(ndata,nutrient,imagebase):
    """Calculate nutrient content in each food group.

    Make a plot of the median values of the chosen nutrient in each food group. Calculate how much of the chosen nutrient each food group contains on average. The median value plot is an example from the book https://github.com/wesm/pydata-book.
    Parameters:
    ndata: [DataFrame object] nutrient database
    nutrient: name of the nutrient which should be included
    imagebase: base name for images
    Return:
    sorted_d: [Series object] food group and how much of the chosen nutrient it contains in descending order
    """
    gr_nutrient_fgroup_qu05 = ndata.groupby(['nutrient', 'fgroup'])['value'].quantile(0.5)  # median value by food group
    plt.clf()
    # use the default ascending=True for plotting 
    sorted = gr_nutrient_fgroup_qu05['%s' % nutrient].sort_values()
    sorted.plot(kind='barh',figsize=(30,15))
    plt.title('Median contents of %s' % nutrient)
    plt.savefig('%s%s.png' % (imagebase,nutrient))
    # return values sorted in descending order
    sorted_d = sorted.sort_values(ascending=False)
    return sorted_d


def recommend_food_groups(ndata,imagebase,does_have=[],does_not_have=[],strictness_good=10,strictness_bad=10):
    """Recommend food groups based on which nutrients it should or should not have.

    Parameters:
    ndata: [DataFrame object] nutrient database
    imagebase: base name for images
    does_have: list of nutrients which should be included
    does_not_have: list of nutrients which should not be included
    strictness_good: take into account this number of good products (for each nutrient)
    strictness_bad: take into account this number of bad products (for each nutrient)
    """
    good_ones = []
    bad_ones = []

    for nutr in does_have:
        print(nutr)
        good = food_group_nutrient_contents(ndata,nutr,imagebase)
        good_lim = good[:strictness_good]
        print('Food groups having relatively much %s:' % nutr)
        print(good_lim)
        good_ones.append(good_lim)

    for nutr in does_not_have:
        print(nutr)
        bad = food_group_nutrient_contents(ndata,nutr,imagebase)
        bad_lim = bad[:strictness_bad]
        print('Food groups having relatively much %s:' % nutr)
        print(bad_lim)
        bad_ones.append(bad_lim)

    goods = len(good_ones)
    bads = len(bad_ones)

    # make a union of good ones and bad ones, respectively
    good_union = pd.Index([])
    if(goods > 0):
        for i in range(goods):
            good_union = good_union.union(good_ones[i].index)
            
    bad_union = pd.Index([])
    if(bads > 0):
        for i in range(bads):
            bad_union = bad_union.union(bad_ones[i].index)

    # remove the bad ones from the good ones
    if(goods > 0 and bads > 0):
        final_goods = good_union.difference(bad_union)
        recoms = len(final_goods)
        print('Recommended food groups (%d):' % recoms)
        if(recoms > 0):
            prettyprint = pd.DataFrame(final_goods)
            print(prettyprint)
    elif(goods > 0):
        final_goods = good_union
        recoms = len(final_goods)
        print('Recommended food groups (%d):' % recoms)
        if(recoms > 0):
            prettyprint = pd.DataFrame(final_goods)
            print(prettyprint)
    elif(bads > 0):  # only list avoidable food groups
        avoids = len(bad_union)
        print('Avoid food groups (%d):' % avoids)
        if(avoids > 0):
            prettyprint = pd.DataFrame(bad_union)
            print(prettyprint)
    else:
        print('Give lists of good and bad nutrients')



# test example
ndata = initialize_nutrient_data()

does_have = ['Protein','Fiber, total dietary']
does_not_have = ['Starch','Energy']

imagebase = 'food_groups_'

strictness_good = 12  # take into account this number of good products (for each nutrient)
strictness_bad = 5    # take into account this number of bad products (for each nutrient)

recommend_food_groups(ndata,imagebase,does_have=does_have,does_not_have=does_not_have,strictness_good=strictness_good,strictness_bad=strictness_bad)

