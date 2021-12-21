import pandas as pd
import numpy as np
import datetime
import os
import sys
import datetime
import time
import math
import shutil
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.formula.api as sm

yticks = [75,80,85,90,95,100,105,110]
xticks = [2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020]

#---------------------------------------------------------------------------
#File with information for program!
control = pd.read_csv(
    'control.txt',
    header=None,
    names=['control']
    )
#---------------------------------------------------------------------------
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
    df = df.loc[(df['id'] > 0)]
    df1 = pd.merge(left=df1, right=df[['id',traitcount]], on='id', how='outer').fillna(0, downcast='infer')
    return df1

def plottingmean(axes,df,trait,ylabel,name):
    axes.set_xticks(xticks)
    axes.set_xticklabels(xticks, rotation=90)
#    axes.set_yticks(yticks)
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
    axes.set_ylim(75,110)
    result = sm.ols(formula= f'{x} ~ {y}', data=df).fit()
    params = pd.DataFrame(result.params)
    slope = params.loc[params.index == y]
    axes.annotate(f"{params.loc[params.index == y]}",xy = (2006,82))
    axes.annotate(f"rsquared = {result.rsquared:.2f}",xy = (2006,76))


def yieldcollect(filename,t1,t2,t3,ID, per1, per2, per3, pe1, pe2, pe3):
    file=filename
    LP0, LP1, LP2, LP3, Wil=[], [], [], [], []
    # Note that the data only includes milk records from dim=5 to dim=305.
    # Therefore, (I think) I should only use these 301 days (rather than 305)
    # Note that these computations are according to appendix G in Mrode.
    for i in range(1,301):
        stdim = -1+2*(i-1)/(305-5)
        LP0.append(1)
        LP1.append(1.2247*stdim)       #L2 in solsaman
        LP2.append(-0.7906 + 2.3717*stdim**2)   #L3 in solsaman
        LP3.append(-2.8067*stdim + 4.6771*stdim**3)   #L4 in solsaman
        Wil.append(math.exp(-0.05*(i+5)))
    Lyield = np.array([LP0,LP1,Wil])
    Lper1 = Lyield[:,45]
    Lper2 = Lyield[:,46:246]
    Lyield = Lyield[:,5:296]
    Lscs = np.array([LP0,LP1,LP2])
    Lscs = Lscs[:,5:296]
    Lfppp = np.array([LP0,LP1,LP2,LP3])
    Lfppp = Lfppp[:,5:296]

    a1,a2,a3,a4,b1,b2,b3,b4,c1,c2,c3,c4,ShortID = [],[],[],[],[],[],[],[],[],[],[],[],[]
    d1,d2,d3,d4,e1,e2,e3,e4,f1,f2,f3,f4, ShortIDPE = [],[],[],[],[],[],[],[],[],[],[],[],[]

    if ((t1 == 'my1')
        | (t1 == 'fy1')
        | (t1 == 'py1')
        | (t1 == 'scs1')
        ):
        print(f'Reading {filename}')
        with open(file, 'r') as fp:
             for cnt, line in enumerate(fp):
                 if line[0]=='4':
                    i = line.split()
                    if i[2]=='1':
                        ShortID.append(i[4])
                        a1.append(float(i[7]))
                    elif i[2]=='2':
                        a2.append(float(i[7]))
                    elif i[2]=='3':
                        a3.append(float(i[7]))
                    elif i[2]=='4':
                        b1.append(float(i[7]))
                    elif i[2]=='5':
                        b2.append(float(i[7]))
                    elif i[2]=='6':
                        b3.append(float(i[7]))
                    elif i[2]=='7':
                        c1.append(float(i[7]))
                    elif i[2]=='8':
                        c2.append(float(i[7]))
                    elif i[2]=='9':
                        c3.append(float(i[7]))
                 elif line[0]=='3' and line.split()[3]=='3':
                    i = line.split()
                    if i[2]=='1':
                        ShortIDPE.append(i[4])
                        d1.append(float(i[7]))
                    elif i[2]=='2':
                        d2.append(float(i[7]))
                    elif i[2]=='3':
                        d3.append(float(i[7]))
                    elif i[2]=='4':
                        e1.append(float(i[7]))
                    elif i[2]=='5':
                        e2.append(float(i[7]))
                    elif i[2]=='6':
                        e3.append(float(i[7]))
                    elif i[2]=='7':
                        f1.append(float(i[7]))
                    elif i[2]=='8':
                        f2.append(float(i[7]))
                    elif i[2]=='9':
                        f3.append(float(i[7]))

        aList = np.array([a1, a2, a3],dtype='float')
        bList = np.array([b1, b2, b3],dtype='float')
        cList = np.array([c1, c2, c3],dtype='float')
        dList = np.array([d1, d2, d3],dtype='float')
        eList = np.array([e1, e2, e3],dtype='float')
        fList = np.array([f1, f2, f3],dtype='float')

        print(file, 'read.')
        print('Solutions written to numpy arrays a, b and c.')

        if ((t1 == 'my1')
            | (t1 == 'fy1')
            | (t1 == 'py1')
            ):

            #Multiply matrices to obtain solutions.
            vec1 = np.dot(aList.transpose(), Lyield)
            vec2 = np.dot(bList.transpose(), Lyield)
            vec3 = np.dot(cList.transpose(), Lyield)
            vec4 = np.dot(dList.transpose(), Lyield)
            vec5 = np.dot(eList.transpose(), Lyield)
            vec6 = np.dot(fList.transpose(), Lyield)

            #Multiply matrices to obtain solutions.
            vec1per = (-(np.dot(aList.transpose(), Lper1))*200) + (np.dot(aList.transpose(), Lper2)).sum(1)
            vec2per = (-(np.dot(bList.transpose(), Lper1))*200) + (np.dot(bList.transpose(), Lper2)).sum(1)
            vec3per = (-(np.dot(cList.transpose(), Lper1))*200) + (np.dot(cList.transpose(), Lper2)).sum(1)

            df = list(zip(ShortID, vec1.sum(1), vec2.sum(1), vec3.sum(1), vec1per, vec2per, vec3per))
            df = sorted(df, key = lambda id: id[0])
            dfPE = list(zip(ShortIDPE, vec4.sum(1), vec5.sum(1), vec6.sum(1)))
            dfPE = sorted(dfPE, key = lambda id: id[0])

            df = pd.DataFrame(df, columns = [ID, t1, t2, t3, per1, per2, per3])
            df[ID] = df[ID].astype(int)

            dfPE = pd.DataFrame(dfPE, columns = [ID, pe1, pe2, pe3])
            dfPE[ID] = dfPE[ID].astype(int)

            df = pd.merge(df,dfPE, how="left", on = ID)

        if (t1 == 'scs1'):

            #Multiply matrices to obtain solutions.
            vec1 = np.dot(aList.transpose(), Lscs)
            vec2 = np.dot(bList.transpose(), Lscs)
            vec3 = np.dot(cList.transpose(), Lscs)
            vec4 = np.dot(dList.transpose(), Lscs)
            vec5 = np.dot(eList.transpose(), Lscs)
            vec6 = np.dot(fList.transpose(), Lscs)

            df = list(zip(ShortID, vec1.sum(1), vec2.sum(1), vec3.sum(1)))
            df = sorted(df, key = lambda id: id[0])
            dfPE = list(zip(ShortIDPE, vec4.sum(1), vec5.sum(1), vec6.sum(1)))
            dfPE = sorted(dfPE, key = lambda id: id[0])

            df = pd.DataFrame(df, columns = [ID, t1, t2, t3])
            df[ID] = df[ID].astype(int)

            dfPE = pd.DataFrame(dfPE, columns = [ID, pe1, pe2, pe3])
            dfPE[ID] = dfPE[ID].astype(int)

            df = pd.merge(df,dfPE, how="left", on = ID)

    if ((t1 == 'pp1')
        | (t1 == 'fp1')
        ):
        print(f'Reading {filename}')
        with open(file, 'r') as fp:
             for cnt, line in enumerate(fp):
                 if line[0]=='4':
                    i = line.split()
                    if i[2]=='1':
                        ShortID.append(i[4])
                        a1.append(float(i[7]))
                    elif i[2]=='2':
                        a2.append(float(i[7]))
                    elif i[2]=='3':
                        a3.append(float(i[7]))
                    elif i[2]=='4':
                        a4.append(float(i[7]))
                    elif i[2]=='5':
                        b1.append(float(i[7]))
                    elif i[2]=='6':
                        b2.append(float(i[7]))
                    elif i[2]=='7':
                        b3.append(float(i[7]))
                    elif i[2]=='8':
                        b4.append(float(i[7]))
                    elif i[2]=='9':
                        c1.append(float(i[7]))
                    elif i[2]=='10':
                        c2.append(float(i[7]))
                    elif i[2]=='11':
                        c3.append(float(i[7]))
                    elif i[2]=='12':
                        c4.append(float(i[7]))
                 elif line[0]=='3' and line.split()[3]=='3':
                    i = line.split()
                    if i[2]=='1':
                        ShortIDPE.append(i[4])
                        d1.append(float(i[7]))
                    elif i[2]=='2':
                        d2.append(float(i[7]))
                    elif i[2]=='3':
                        d3.append(float(i[7]))
                    elif i[2]=='4':
                        d4.append(float(i[7]))
                    elif i[2]=='5':
                        e1.append(float(i[7]))
                    elif i[2]=='6':
                        e2.append(float(i[7]))
                    elif i[2]=='7':
                        e3.append(float(i[7]))
                    elif i[2]=='8':
                        e4.append(float(i[7]))
                    elif i[2]=='9':
                        f1.append(float(i[7]))
                    elif i[2]=='10':
                        f2.append(float(i[7]))
                    elif i[2]=='11':
                        f3.append(float(i[7]))
                    elif i[2]=='12':
                        f4.append(float(i[7]))

        aList = np.array([a1, a2, a3, a4],dtype='float')
        bList = np.array([b1, b2, b3, b4],dtype='float')
        cList = np.array([c1, c2, c3, c4],dtype='float')
        dList = np.array([d1, d2, d3, d4],dtype='float')
        eList = np.array([e1, e2, e3, e4],dtype='float')
        fList = np.array([f1, f2, f3, f4],dtype='float')

        print(file, 'read.')
        print('Solutions written to numpy arrays a, b and c.')

        #Multiply matrices to obtain solutions.
        vec1 = np.dot(aList.transpose(), Lfppp)
        vec2 = np.dot(bList.transpose(), Lfppp)
        vec3 = np.dot(cList.transpose(), Lfppp)
        vec4 = np.dot(dList.transpose(), Lfppp)
        vec5 = np.dot(eList.transpose(), Lfppp)
        vec6 = np.dot(fList.transpose(), Lfppp)

        df = list(zip(ShortID, vec1.sum(1), vec2.sum(1), vec3.sum(1)))
        df = sorted(df, key = lambda id: id[0])
        dfPE = list(zip(ShortIDPE, vec4.sum(1), vec5.sum(1), vec6.sum(1)))
        dfPE = sorted(dfPE, key = lambda id: id[0])

        df = pd.DataFrame(df, columns = [ID, t1, t2, t3])
        df[ID] = df[ID].astype(int)

        dfPE = pd.DataFrame(dfPE, columns = [ID, pe1, pe2, pe3])
        dfPE[ID] = dfPE[ID].astype(int)

        df = pd.merge(df,dfPE, how="left", on = ID)

    return df
