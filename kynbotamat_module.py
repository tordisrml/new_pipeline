import pandas as pd
import numpy as np
import datetime
import os
import shutil
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.formula.api as sm

#---------------------------------------------------------------------------
#Function to read tab/space sperated files with files
#---------------------------------------------------------------------------
def readfilecsv(file, columns, sepr):
    print( f'-------------------------------------------' )
    print( f'Reading {file} file .....' )
    print( f'-------------------------------------------' )
    df = pd.read_csv(
        file,
        header=None,
        sep=sepr,
        names=columns)
    return df


#---------------------------------------------------------------------------
#Function to read fixed with files
#---------------------------------------------------------------------------
def readingfilefwf(file,columns,filewidths):
    print( f'----------------------------------------' )
    print( f'Reading the {file} file.....' )
    print( f'----------------------------------------' )
    widths = filewidths
    df = pd.read_fwf(
        file,
        header=None,
        widths=widths,
        names=columns)
    print( f'-------------------------------------------' )
    print( f'Reading the the {file} file finished.....' )
    print( f'-------------------------------------------' )
    return df

#---------------------------------------------------------------------------
#This is a function to create herd x year group fixed effect
#If cows are fewer than 3 in a group the function tries to combine with group
#above or below WITHIN SAME HERD
#---------------------------------------------------------------------------
def hy_grouping(element, df, col, group, df1):
    tdf = []
    for herd, data in df.groupby( element ):
        # get counts and assign initial groups
        counts = data[ col ].value_counts().sort_index().to_frame()
        counts[ 'group' ] = range( counts.shape[ 0 ] )

        while True:
            gcounts = counts.groupby( 'group' ).sum()[ col ]  # group counts
            change = gcounts[ gcounts.values < 4 ]  # groups with too few

            if change.shape[ 0 ] == 0:
                # no changes, exit
                break

            # check how to merge groups
            cgroup = change.index.min()
            groups = gcounts.index.values
            g_ind = list( groups ).index( cgroup )
            if ( g_ind + 1 ) < groups.shape[ 0 ]:
                # merge forward
                ngroup = groups[ g_ind + 1 ]

            elif g_ind > 0:
                # merge backward
                ngroup = groups[ g_ind - 1 ]

            else:
                # no groups to merge
                print( f'Can not merge herd {herd}' )
                break

            counts.loc[ counts[ 'group' ] == cgroup, 'group' ] = ngroup

        # assign groups
        for ind, gdata in counts.iterrows():
            data.loc[ data[ col ] == ind, 'group' ] = gdata[ 'group' ]

        tdf.append( data )

    tdf = pd.concat( tdf )
    #Creation of fixed effect herd + yeargroup
    tdf[ group ] = tdf[ element ].astype( 'str' ) + tdf[ 'group' ].astype( int ).astype( str )
    #Merged into main dataframe
    df1 = pd.merge(left=df1, right=tdf[['id', group]], on='id',
        how='outer').fillna(0, downcast='infer')

    return df1
#---------------------------------------------------------------------------

#Function to count how many daughters have observations for each lactation
def countingoff (df0,trait,traitcount, df1):
    print( f'-------------------------------------------' )
    print( f'Counting offpspring for {trait} .....' )
    print( f'-------------------------------------------' )
    df = df0.loc[(df0[trait] >= 0)]
    #Counting offspring for this trait
    df.loc[:,traitcount] = df.groupby('sire')['sire'].transform('count') #Counting offspring for this trait
    df = df[['sire', traitcount]]
    df = df.drop_duplicates(subset=['sire'])
    df.columns = ['id', traitcount]
    df1 = pd.merge(left=df1, right=df[['id',traitcount]], on='id', how='outer').fillna(0, downcast='infer')
    return df1

xticks = [2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020]
yticks = [75,80,85,90,95,100,105,110]

def plottingmean(axes,df,trait,ylabel,name):
    axes.set_xticks(xticks)
    axes.set_xticklabels(xticks, rotation=90)
    axes.set_yticks(yticks)
    axes.plot(df.groupby('BY')[trait].mean(), label=name, marker='.')
    axes.legend()
    axes.set_ylabel(ylabel)
    axes.grid(axis='y')
    axes.set_xlim(2000,2020)

def plottingmeansns(axes,df,x,y,ylabel,name):
    sns.regplot(x=x, y=y, data=df, ax=axes, label=name)
    axes.set_xticks(xticks)
    axes.set_xticklabels(xticks, rotation=90)
    axes.set_yticks(yticks)
    axes.legend()
    axes.set_ylabel(ylabel)
    axes.grid(axis='y')
    axes.set_xlim(2000,2020)
    result = sm.ols(formula= f'{x} ~ {y}', data=df).fit()
    params = pd.DataFrame(result.params)
    slope = params.loc[params.index == y]
    # axes.set_title(f'{(params.loc[params.index == y].astype(str))},rsquared = {result.rsquared:.2f} ', fontsize=8)
    axes.annotate(f"{params.loc[params.index == y]}",xy = (2006,82))
    axes.annotate(f"rsquared = {result.rsquared:.2f}",xy = (2006,76))
