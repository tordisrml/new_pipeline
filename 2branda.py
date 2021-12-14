#This is a program to collect EBV results for Icelandic dairy cows and
#Write to a file to be read by Huppa

#Written by Thordis Thorarinsdottir 2021

import pandas as pd
import numpy as np
import os
import shutil
from matplotlib import pyplot as plt
import seaborn as sns
import statsmodels.formula.api as sm
#user defined functions used by program from kynbotamat_module.py
from kynbotamat_module import readfilecsv
from kynbotamat_module import readingfilefwf
from kynbotamat_module import countingoff
from kynbotamat_module import plottingmean
from kynbotamat_module import plottingmeansns
from kynbotamat_module import yieldcollect

#Option to read SOL files for these 3 traits and collect results, 0 means skip
# and 1 means SOL files will be read
yield_results = 1
fertility_results = 1
conf_results = 1
rankorder_results = 1
long_results = 1

#Option to write scaled results for above traits to seperate files to be read
#later, 0 means skip and 1 means seperate result files will be written to disc
seperate_files = 1

#Option to collect results for a large datafile to be read by Huppa.
#--IF ABOVE OPTIONS ARE SET TO 0 --> program will collect results for these three
    #trait from seperate files!
#--IF ABOVE OPTIONS ARE SET TO 1 --> program will collect results for these three
    #trait straight from the program itself.
#Program always collects results for yield, scs, persistancy and longevity from
    #imported datafiles.
collectresults = 1    #0 means don't, 1 means collect

#Option to write large datafile to disc
writebranda = 1

collectdmufiles = 0

phantomcollection = 0

plotting = 0

#---------------------------------------------------------------------------
#File with information for program!
control = pd.read_csv(
    'control.txt',
    header=None,
    names=['control']
    )
#---------------------------------------------------------------------------
#Info for program
#---------------------------------------------------------------------------
yearmonth = control.loc[15,'control']
brandafile = control.loc[23,'control']
#Scaling year, present year - 5 years
scalingyear = pd.to_numeric(control.loc[12,'control'])
#Scaling objects
imean = pd.to_numeric(control.loc[13,'control'])
isd = pd.to_numeric(control.loc[14,'control'])
#---------------------------------------------------------------------------
#Radnrkodifile
radnrkodifile = '../dmu_data/radnrkodi' #Created by prep_tdm.f
#Format of radnrkodi file created by prep_tdm.f
radnrkodi_columns = ['id','code_id','stada','norec','fix1','fix2', 'fix3','sex']
widths_radnrkodi = [15,9,3,3,6,6,6,2]
#---------------------------------------------------------------------------
#Name of pedigree file
pedigreefile = control.loc[11,'control']
#Format of pedigree file from Huppa
ped_columns = ['id','dam','sire','unused','sex','unused2', 'bullno', 'name', 'farm']
ped_widths = [15,15,15,12,1,2,5,20,20]
#---------------------------------------------------------------------------

#DMU sol files
mysol = f'../DMU/{yearmonth}/my/SOL'
fysol = f'../DMU/{yearmonth}/fy/SOL'
pysol = f'../DMU/{yearmonth}/py/SOL'
scssol = f'../DMU/{yearmonth}/scs/SOL'
fpsol = f'../DMU/{yearmonth}/fp/SOL'
ppsol = f'../DMU/{yearmonth}/pp/SOL'

fertilitysolfile = f'../DMU/{yearmonth}/fer/SOL'
rankordersolfile = f'../DMU/{yearmonth}/rank/SOL'
confsolfile = f'../DMU/{yearmonth}/conf/SOL'
longsolfile = f'../DMU/{yearmonth}/long/SOL'
#Format of SOL files
solcolumns = ['1_code_effect', #2 for fixed and 4 for genetic
    '2_trait_no',   # lact 1, 2 or 3
    '3','4',
    'code_id',   #fixed effects and id's
    '6_no_obs',  #No. of observations in this class
    '7',
    '8_BLUP',  #Estimate/prediction
    '9']   #Solution from the second but last DMU5
sol_widths = [1,3,3,4,12,12,12,20,20]
#---------------------------------------------------------------------------

#Observationfiles used for DMU
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
fertilityobs = '../dmu_data/dmu_fertility.txt'
fertility_columns = ['code_id','HBY','HC1','HC2','HC3','IYM0','IYM1','IYM2',
    'IYM3','AGEi_h','AGEc_1','AGEc_2','AGEc_3','tech_h',
    'CR0','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']
#---------------------------------------------------------------------------
confobs = '../dmu_data/dmu_conformation.txt'
conf_columns = ['code_id','HdomsY','lact','AGEc_1',
    'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap']
#---------------------------------------------------------------------------
rankorderobs = '../dmu_data/dmu_rankorder.txt'
rankorder_columns = ['code_id', 'year', 'mjaltarod', 'gaedarod']
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
longobs = '../dmu_data/dmu_long.txt'
longobs_columns = ['code_id','AGEc_1','CYM1','h5y','herdCY1','L1','L2','L3']
#---------------------------------------------------------------------------

# #Seperate EBV result files
yieldebv = '../results/yieldebv.txt'
fertilityebv = '../results/fertilityebv.txt' #written by this program if option above set to 1
confebv = '../results/conformationebv.txt' #written by this program if option above set to 1
rankorderebv = '../results/rankorderebv.txt' #written by this program if option above set to 1
longebv = '../results/longebv.txt'        #written by this program if option above set to 1
#Accuracy files for yield and scs
accyield = '../results/accuracy.sol' #From programs by JHE
accscs = '../results/accuracy_f.sol' #From programs by JHE

#---------------------------------------------------------------------------
#Columns in EBV files
#---------------------------------------------------------------------------
accyield_columns = ['id','no_daughters_yield', 'yield_acc']
widths_accyield = [15,8,8]

accscs_columns =  ['id','no_daughters_SCS', 'SCS_acc']
widths_accscs = [15,8,8]
#---------------------------------------------------------------------------
#Space seperated files/dataframes
#(Created by this program if option above set to 1)
yieldebv_columns = ['id','my1','my2','my3','fy1','fy2','fy3','py1','py2','py3','fp1','fp2','fp3',
    'pp1','pp2','pp3','scs1','scs2','scs3', 'milkper', 'fatper', 'protper' ]

fertilityebv_columns = ['id','fer_lact1','fer_lact2','fer_lact3','CR0','ICF',
    'IFL','fertility','offCR0','offICF1','offICF2','offICF3']

confebv_columns = ['id','boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap','offconf']

rankorderebv_columns = ['id','mjaltarod', 'gaedarod','offrankorder']

longebv_columns = ['id','L1','L2','L3','offlong1','offlong2','offlong3']


#Phantom group parent results
fertilityebvphg = '../results/fertilityebvphg.txt' #written by this program if option above set to 1
confebvphg = '../results/conformationebvphg.txt' #written by this program if option above set to 1
rankorderebvphg = '../results/rankorderebvphg.txt' #written by this program if option above set to 1
longebvphg = '../results/longebvphg.txt' #written by this program if option above set to 1
#---------------------------------------------------------------------------

#Large file with all solutions scaled
brandafile_columns = ['id',                                             #1
        'my1','my2','my3','fy1','fy2','fy3', #6
        'py1','py2','py3','fp1','fp2','fp3',    #6
        'pp1','pp2','pp3','scs1','scs2','scs3',             #6
        'milkper', 'fatper', 'protper',                                 #3
        'fer_lact1','fer_lact2','fer_lact3','CR0','ICF','IFL',          #6
        'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti', #6
        'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta', #4
        'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',  #5
        'mjaltir', 'skap','mjaltarod', 'gaedarod','L1','L2','L3',              #5
        'myt','fyt','pyt','fpt','ppt',                 #5
        'yieldtotal','fertility','scs','jugur','spenar','mjaltir_t','skap2', #7
        'total',                                                             #1
        'no_daughters_yield', 'yield_acc', 'no_daughters_SCS', 'SCS_acc', #4
        'offconf','offrankorder','offCR0','offICF1','offICF2','offICF3', #6
        'offlong1','offlong2','offlong3'                                        #1
        ] #72 columns
widths_branda = [15,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, #31
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] #41
#---------------------------------------------------------------------------


#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Functions used by program!
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#Function to seperate solutions in SOL file and merging with radnrkodi
def solread(df):
    print( f'----------------------------------------' )
    print( f'Merging sol file with radnrkodi .....' )
    print( f'----------------------------------------' )
    #Only keep genetic effectss
    solp = df[(df['1_code_effect'] == 4) & (df['code_id'] < 0)]
    sol = df[(df['1_code_effect'] == 4) & (df['code_id'] > 0)] #collect animal results
    # merge real ids back into results
    # sol= pd.merge(left=sol[
    #     ['code_id','2_trait_no', '6_no_obs','8_BLUP']
    #     ], right=radnrkodi[['id','code_id']], on='code_id', how='left')
    sol= pd.merge(left=radnrkodi[
        ['id','code_id']
        ], right=sol[['code_id','2_trait_no', '6_no_obs','8_BLUP']], on='code_id', how='left')
    return sol, solp
#Function to find solution for each trait in sol dataframes and merge them by ids
#into a new df
def solutions(df, df1, no, trait, id):
    print( f'----------------------------------------' )
    print( f'Collecting results for {trait} .....' )
    print( f'----------------------------------------' )
    dftrait = df.loc[df['2_trait_no'] == no ] #one dataframe per trait
    dftrait.loc[:,trait] = dftrait['8_BLUP']            #name the trait
    df1 = pd.merge(left=df1, right=dftrait[[id,trait]  #merge into large dataframe
        ], on=id)
    df1 = df1.drop_duplicates(subset=[id])
    df1 = df1.sort_values(by=[id])
    return df1
#---------------------------------------------------------------------------

#-------------------------------------------------------------
#Creating STD and mean for traits from average group
#-------------------------------------------------------------
#Reading in the file used in DMU so animals that have their own observations
#can be located
def ownobs(ownobsfile, columns):
    print( f'-------------------------------------------' )
    print( f'Finding own obs for {ownobsfile} file .....' )
    print( f'-------------------------------------------' )
    df = pd.read_csv(
        ownobsfile,
        header=None,
        sep=' ',
        names=columns)
    #Merging ownobs and id file to bring back einstaklingsnumer used in huppa
    df = pd.merge(left=df, right=radnrkodi[['id','code_id']], on='code_id',  how='inner')
    return df
#-------------------------------------------------------------
#Function to create groups for average calculations
def avegroup(df,df1):
    #Make birth year from id number an integer
    print( f'-------------------------------------------' )
    print( f'Creating average groups.....' )
    print( f'-------------------------------------------' )
    df['BY'] = (df.id.astype(str).str[:4]).astype(int)
    ave_group = df.loc[(df['BY'] == scalingyear)]
    ave_group = pd.merge(left=df1['id'], right=ave_group, on='id', how='left')
    return ave_group

#Funtion to scale traits
#-------------------------------------------------------------
def scaling(df,trait,x,ave_group):
    print( f'-------------------------------------------' )
    print( f'Scaling .....' )
    print( f'-------------------------------------------' )
    SD = ave_group[trait].std()
    mean = ave_group[trait].mean()
    print(trait)
    print(mean)
    print(SD)
    df[trait] = (imean+(((df[trait]- mean )*(x) / SD) * isd ))#.astype(float).round().astype(int)
    return df[trait]
#-------------------------------------------------------------
#Function to read imported result files and combine to a big one
#Fixed with files
def combineresultsfwf(df, file, columns,widths_file):
    print( f'-------------------------------------------' )
    print( f'Reading {file} file .....' )
    print( f'-------------------------------------------' )
    df2 = pd.read_fwf(
        file,
        widths=widths_file,
        header=None,
        names=columns)
    #Merging ownobs and id file to bring back einstaklingsnumer used in huppa
    df = pd.merge(left=df, right=df2, on='id',  how='left')
    df = df.drop_duplicates(subset=['id'])
    df = df.sort_values(by=['id'])
    return df

#-------------------------------------------------------------
#Function to read imported result files and combine to a big one
#Space seperated files
def combineresultscsv(df, file, columns):
    print( f'-------------------------------------------' )
    print( f'Reading {file} file .....' )
    print( f'-------------------------------------------' )
    df2 = pd.read_csv(
        file,
        header=None,
        sep=' ',
        names=columns)
    #Merging ownobs and id file to bring back einstaklingsnumer used in huppa
    df = pd.merge(left=df, right=df2, on='id',  how='left')
    df = df.drop_duplicates(subset=['id'])
    df = df.sort_values(by=['id'])
    return df

#Function to combine results if results are read from dataframes
def combineresultsdf(df, df2,columns):
    print( f'-------------------------------------------' )
    print( f'Combining .....' )
    print( f'-------------------------------------------' )
    df = pd.merge(left=df, right=df2[columns], on='id',  how='left')
    df = df.drop_duplicates(subset=['id'])
    df = df.sort_values(by=['id'])
    return df


#-------------------------------------------------------------
#-------------------------------------------------------------
#START OF PROGRAM!!!
#-------------------------------------------------------------
#-------------------------------------------------------------
#Reading in id codes to replace in SOL files
if ((collectresults == 1)
    | (fertility_results == 1)
    | (conf_results == 1)
    | (rankorder_results == 1)
    | (long_results == 1)
    | (yield_results == 1)
    ):
    radnrkodi = readingfilefwf(radnrkodifile,radnrkodi_columns,widths_radnrkodi)
    radnrkodi = radnrkodi.drop(
        ['norec','fix1','fix2', 'fix3','sex'], axis = 1)

    #Reading pedigree
    ped = readingfilefwf(pedigreefile,ped_columns,ped_widths)

#Reading yield SOL files, collecting and scaling results
#Write results to disc if option above set to 1
if yield_results == 1:

    myebv = yieldcollect(mysol,'my1','my2','my3','code_id','myper1','myper2','myper3')
    yielddf = pd.merge(left=radnrkodi, right=myebv, on='code_id', how='left')

    fyebv = yieldcollect(fysol,'fy1','fy2','fy3','code_id','fyper1','fyper2','fyper3')
    yielddf = pd.merge(left=yielddf, right=fyebv, on='code_id', how='left')

    pyebv = yieldcollect(pysol,'py1','py2','py3','code_id','pyper1','pyper2','pyper3')
    yielddf = pd.merge(left=yielddf, right=pyebv, on='code_id', how='left')

    scsebv = yieldcollect(scssol,'scs1','scs2','scs3','code_id','','','')
    yielddf = pd.merge(left=yielddf, right=scsebv, on='code_id', how='left')

    fpebv = yieldcollect(fpsol,'fp1','fp2','fp3','code_id','','','')
    yielddf = pd.merge(left=yielddf, right=fpebv, on='code_id', how='left')

    ppebv = yieldcollect(ppsol,'pp1','pp2','pp3','code_id','','','')
    yielddf = pd.merge(left=yielddf, right=ppebv, on='code_id', how='left')

    yieldscale = yielddf.loc[yielddf['stada'] == 1]

    yielddf['my1'] = scaling(yielddf,'my1',1, yieldscale)
    yielddf['my2'] = scaling(yielddf,'my2',1, yieldscale)
    yielddf['my3'] = scaling(yielddf,'my3',1, yieldscale)

    yielddf['myper1'] = scaling(yielddf,'myper1',1, yieldscale)
    yielddf['myper2'] = scaling(yielddf,'myper2',1, yieldscale)
    yielddf['myper3'] = scaling(yielddf,'myper3',1, yieldscale)

    yielddf['fy1'] = scaling(yielddf,'fy1',1, yieldscale)
    yielddf['fy2'] = scaling(yielddf,'fy2',1, yieldscale)
    yielddf['fy3'] = scaling(yielddf,'fy3',1, yieldscale)

    yielddf['fyper1'] = scaling(yielddf,'fyper1',1, yieldscale)
    yielddf['fyper2'] = scaling(yielddf,'fyper2',1, yieldscale)
    yielddf['fyper3'] = scaling(yielddf,'fyper3',1, yieldscale)

    yielddf['py1'] = scaling(yielddf,'py1',1, yieldscale)
    yielddf['py2'] = scaling(yielddf,'py2',1, yieldscale)
    yielddf['py3'] = scaling(yielddf,'py3',1, yieldscale)

    yielddf['pyper1'] = scaling(yielddf,'pyper1',1, yieldscale)
    yielddf['pyper2'] = scaling(yielddf,'pyper2',1, yieldscale)
    yielddf['pyper3'] = scaling(yielddf,'pyper3',1, yieldscale)

    yielddf['fp1'] = scaling(yielddf,'fp1',1, yieldscale)
    yielddf['fp2'] = scaling(yielddf,'fp2',1, yieldscale)
    yielddf['fp3'] = scaling(yielddf,'fp3',1, yieldscale)
    yielddf['pp1'] = scaling(yielddf,'pp1',1, yieldscale)
    yielddf['pp2'] = scaling(yielddf,'pp2',1, yieldscale)
    yielddf['pp3'] = scaling(yielddf,'pp3',1, yieldscale)
    yielddf['scs1'] = scaling(yielddf,'scs1',-1, yieldscale)
    yielddf['scs2'] = scaling(yielddf,'scs2',-1, yieldscale)
    yielddf['scs3'] = scaling(yielddf,'scs3',-1, yieldscale)

    yielddf['milkper'] = (yielddf['myper1'] * 0.5 +
                            yielddf['myper2'] * 0.3 +
                            yielddf['myper3'] * 0.2 )

    yielddf['fatper'] = (yielddf['fyper1'] * 0.5 +
                            yielddf['fyper2'] * 0.3 +
                            yielddf['fyper3'] * 0.2 )

    yielddf['protper'] = (yielddf['pyper1'] * 0.5 +
                            yielddf['pyper2'] * 0.3 +
                            yielddf['pyper3'] * 0.2 )

    if seperate_files == 1:
        print('Yield results written to seperate file')
        yield_results = yielddf[yieldebv_columns].astype(float).round(2)#.astype(int)
        yield_results.to_csv(yieldebv, index=False, header=False, sep=' ')
        print(yielddf.iloc[50000:50015])
        print(yielddf.info())
        print(f'Yield results written to {yieldebv}')

    print(yielddf.iloc[500000:500015])
    print(yielddf.info())

#Reading fertilty SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if fertility_results == 1:
    #Reading sol files
    fertilitysol = readingfilefwf(fertilitysolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    fertilitysolids, fertilitysolph = solread(fertilitysol)
    #Seperating solutions by traits
    fertilitydf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    fertilitydf = solutions(fertilitysolids, fertilitydf, 1, 'CR0','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 2, 'ICF1','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 3, 'ICF2','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 4, 'ICF3','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 5, 'IFL1','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 6, 'IFL2','id')
    fertilitydf = solutions(fertilitysolids, fertilitydf, 7, 'IFL3','id')
    #Reading own observations files
    ownobs_fertility = ownobs(fertilityobs,fertility_columns)
    #Creating average groups to scale
    fertility_ave = avegroup(fertilitydf,ownobs_fertility)
    #Scaling results
    fertilitydf['CR0'] = scaling(fertilitydf,'CR0',1, fertility_ave)
    fertilitydf['ICF1'] = scaling(fertilitydf,'ICF1',-1, fertility_ave)
    fertilitydf['ICF2'] = scaling(fertilitydf,'ICF2',-1, fertility_ave)
    fertilitydf['ICF3'] = scaling(fertilitydf,'ICF3',-1, fertility_ave)
    fertilitydf['IFL1'] = scaling(fertilitydf,'IFL1',-1, fertility_ave)
    fertilitydf['IFL2'] = scaling(fertilitydf,'IFL2',-1, fertility_ave)
    fertilitydf['IFL3'] = scaling(fertilitydf,'IFL3',-1, fertility_ave)

    #EBVs by lactations
    fertilitydf['fer_lact1'] = (fertilitydf['IFL1'] * 0.6 + fertilitydf['ICF1'] * 0.4)
    fertilitydf['fer_lact2'] = (fertilitydf['IFL2'] * 0.6 + fertilitydf['ICF2'] * 0.4)
    fertilitydf['fer_lact3'] = (fertilitydf['IFL3'] * 0.6 + fertilitydf['ICF3'] * 0.4)
    #Combination of three lactations to make one trait
    fertilitydf['ICF'] = (fertilitydf['ICF1'] * 0.5 +
                            fertilitydf['ICF2'] * 0.3 + fertilitydf['ICF3'] * 0.2 )
    fertilitydf['IFL'] = (fertilitydf['IFL1'] * 0.5 +
                            fertilitydf['IFL1'] * 0.3 + fertilitydf['IFL2'] * 0.2 )
    #New combined trait for fertility
    fertilitydf['fertility'] = (fertilitydf['CR0'] * 0.2 +
        fertilitydf['ICF'] * 0.3 +
        fertilitydf['IFL'] * 0.5 )

    #Counting daughters with own obs
    ownobs_fertility = pd.merge(left=ownobs_fertility, right=ped[['id','sire']], on='id', how='left')

    fertilitydf = countingoff(ownobs_fertility,'CR0','offCR0', fertilitydf)
    fertilitydf = countingoff(ownobs_fertility,'ICF1','offICF1', fertilitydf)
    fertilitydf = countingoff(ownobs_fertility,'ICF2','offICF2', fertilitydf)
    fertilitydf = countingoff(ownobs_fertility,'ICF3','offICF3', fertilitydf)

    if seperate_files == 1:
        print('Fertility results written to seperate file')
        fer_results = fertilitydf[fertilityebv_columns].astype(float).round(2)#.astype(int)
        fer_results.to_csv(fertilityebv, index=False, header=False, sep=' ')
        print(fertilitydf.iloc[50000:50015])
        print(fertilitydf.info())
        print(f'Fertility results written to {fertilityebv}')

    if phantomcollection == 1:

        fertilityphgdf = fertilitysolph['code_id'].copy()
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 1, 'CR0','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 2, 'ICF1','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 3, 'ICF2','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 4, 'ICF3','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 5, 'IFL1','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 6, 'IFL2','code_id')
        fertilityphgdf = solutions(fertilitysolph, fertilityphgdf, 7, 'IFL3','code_id')
        fertilityphgdf['CR0'] = scaling(fertilityphgdf,'CR0',1, fertility_ave)
        fertilityphgdf['ICF1'] = scaling(fertilityphgdf,'ICF1',-1, fertility_ave)
        fertilityphgdf['ICF2'] = scaling(fertilityphgdf,'ICF2',-1, fertility_ave)
        fertilityphgdf['ICF3'] = scaling(fertilityphgdf,'ICF3',-1, fertility_ave)
        fertilityphgdf['IFL1'] = scaling(fertilityphgdf,'IFL1',-1, fertility_ave)
        fertilityphgdf['IFL2'] = scaling(fertilityphgdf,'IFL2',-1, fertility_ave)
        fertilityphgdf['IFL3'] = scaling(fertilityphgdf,'IFL3',-1, fertility_ave)
        print('Phantom group results collected for fertility')
        print(fertilityphgdf.iloc[0:51])
        print(fertilityphgdf.info())
        if seperate_files == 1:
            print('PHG conformation results written to seperate file')
            fertilityphgdf.to_csv(fertilityebvphg, index=False, header=False, sep=' ')
            print(f'PHG conformation results written to {fertilityebvphg}')
    else:
        print('Phantom group results collected for fertility not collected')
else:
    print('Fertility results not collected')

#Reading conformation SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if conf_results == 1:

    #Reading sol files
    confsol = readingfilefwf(confsolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    confsolids, confsolph = solread(confsol)
    #Seperating solutions by traits
    confdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    confdf = solutions(confsolids, confdf, 1, 'boldypt','id')
    confdf = solutions(confsolids, confdf, 2, 'utlogur','id')
    confdf = solutions(confsolids, confdf, 3, 'yfirlina','id')
    confdf = solutions(confsolids, confdf, 4, 'malabreidd','id')
    confdf = solutions(confsolids, confdf, 5, 'malahalli','id')
    confdf = solutions(confsolids, confdf, 6, 'malabratti','id')
    confdf = solutions(confsolids, confdf, 7, 'stada_haekla_hlid','id')
    confdf = solutions(confsolids, confdf, 8, 'stada_haekla_aftan','id')
    confdf = solutions(confsolids, confdf, 9, 'klaufhalli','id')
    confdf = solutions(confsolids, confdf, 10, 'jugurfesta','id')
    confdf = solutions(confsolids, confdf, 11, 'jugurband','id')
    confdf = solutions(confsolids, confdf, 12, 'jugurdypt','id')
    confdf = solutions(confsolids, confdf, 13, 'spenalengd','id')
    confdf = solutions(confsolids, confdf, 14, 'spenathykkt','id')
    confdf = solutions(confsolids, confdf, 15, 'spenastada','id')
    confdf = solutions(confsolids, confdf, 16, 'mjaltir','id')
    confdf = solutions(confsolids, confdf, 17, 'skap','id')
    #Reading own observations files
    ownobs_conf = ownobs(confobs,conf_columns)
    #Creating average groups to scale
    conf_ave = avegroup(confdf,ownobs_conf)
    #Scaling results
    confdf['boldypt'] = scaling(confdf,'boldypt',1, conf_ave)
    confdf['utlogur'] = scaling(confdf,'utlogur',1, conf_ave)
    confdf['yfirlina'] = scaling(confdf,'yfirlina',1, conf_ave)
    confdf['malabreidd'] = scaling(confdf,'malabreidd',1, conf_ave)
    confdf['malahalli'] = scaling(confdf,'malahalli',1, conf_ave)
    confdf['malabratti'] = scaling(confdf,'malabratti',1, conf_ave)
    confdf['stada_haekla_hlid'] = scaling(confdf,'stada_haekla_hlid',1, conf_ave)
    confdf['stada_haekla_aftan'] = scaling(confdf,'stada_haekla_aftan',1, conf_ave)
    confdf['klaufhalli'] = scaling(confdf,'klaufhalli',1, conf_ave)
    confdf['jugurfesta'] = scaling(confdf,'jugurfesta',1, conf_ave)
    confdf['jugurband'] = scaling(confdf,'jugurband',1, conf_ave)
    confdf['jugurdypt'] = scaling(confdf,'jugurdypt',1, conf_ave)
    confdf['spenalengd'] = scaling(confdf,'spenalengd',1, conf_ave)
    confdf['spenathykkt'] = scaling(confdf,'spenathykkt',1, conf_ave)
    confdf['spenastada'] = scaling(confdf,'spenastada',1, conf_ave)
    confdf['mjaltir'] = scaling(confdf,'mjaltir',1, conf_ave)
    confdf['skap'] = scaling(confdf,'skap',1, conf_ave)

    #Counting daughters with own obs
    ownobs_conf = pd.merge(left=ownobs_conf, right=ped[['id','sire']], on='id', how='left')

    confdf = countingoff(ownobs_conf,'boldypt','offconf', confdf)

    if seperate_files == 1:
        print('Conformation results written to seperate file')
        conf_results = confdf[confebv_columns].astype(float).round(2)#.astype(int)
        conf_results.to_csv(confebv, index=False, header=False, sep=' ')
        print(f'Conformation results written to {confebv}')

        print(confdf.iloc[50000:50015])
        print(confdf.info())

    if phantomcollection == 1:

        confphgdf = confsolph['code_id'].copy()
        confphgdf = solutions(confsolph, confphgdf, 1, 'boldypt','code_id')
        confphgdf = solutions(confsolph, confphgdf, 2, 'utlogur','code_id')
        confphgdf = solutions(confsolph, confphgdf, 3, 'yfirlina','code_id')
        confphgdf = solutions(confsolph, confphgdf, 4, 'malabreidd','code_id')
        confphgdf = solutions(confsolph, confphgdf, 5, 'malahalli','code_id')
        confphgdf = solutions(confsolph, confphgdf, 6, 'malabratti','code_id')
        confphgdf = solutions(confsolph, confphgdf, 7, 'stada_haekla_hlid','code_id')
        confphgdf = solutions(confsolph, confphgdf, 8, 'stada_haekla_aftan','code_id')
        confphgdf = solutions(confsolph, confphgdf, 9, 'klaufhalli','code_id')
        confphgdf = solutions(confsolph, confphgdf, 10, 'jugurfesta','code_id')
        confphgdf = solutions(confsolph, confphgdf, 11, 'jugurband','code_id')
        confphgdf = solutions(confsolph, confphgdf, 12, 'jugurdypt','code_id')
        confphgdf = solutions(confsolph, confphgdf, 13, 'spenalengd','code_id')
        confphgdf = solutions(confsolph, confphgdf, 14, 'spenathykkt','code_id')
        confphgdf = solutions(confsolph, confphgdf, 15, 'spenastada','code_id')
        confphgdf = solutions(confsolph, confphgdf, 16, 'mjaltir','code_id')
        confphgdf = solutions(confsolph, confphgdf, 17, 'skap','code_id')
        confphgdf['boldypt'] = scaling(confphgdf,'boldypt',1, conf_ave)
        confphgdf['utlogur'] = scaling(confphgdf,'utlogur',1, conf_ave)
        confphgdf['yfirlina'] = scaling(confphgdf,'yfirlina',1, conf_ave)
        confphgdf['malabreidd'] = scaling(confphgdf,'malabreidd',1, conf_ave)
        confphgdf['malahalli'] = scaling(confphgdf,'malahalli',1, conf_ave)
        confphgdf['malabratti'] = scaling(confphgdf,'malabratti',1, conf_ave)
        confphgdf['stada_haekla_hlid'] = scaling(confphgdf,'stada_haekla_hlid',1, conf_ave)
        confphgdf['stada_haekla_aftan'] = scaling(confphgdf,'stada_haekla_aftan',1, conf_ave)
        confphgdf['klaufhalli'] = scaling(confphgdf,'klaufhalli',1, conf_ave)
        confphgdf['jugurfesta'] = scaling(confphgdf,'jugurfesta',1, conf_ave)
        confphgdf['jugurband'] = scaling(confphgdf,'jugurband',1, conf_ave)
        confphgdf['jugurdypt'] = scaling(confphgdf,'jugurdypt',1, conf_ave)
        confphgdf['spenalengd'] = scaling(confphgdf,'spenalengd',1, conf_ave)
        confphgdf['spenathykkt'] = scaling(confphgdf,'spenathykkt',1, conf_ave)
        confphgdf['spenastada'] = scaling(confphgdf,'spenastada',1, conf_ave)
        confphgdf['mjaltir'] = scaling(confphgdf,'mjaltir',1, conf_ave)
        confphgdf['skap'] = scaling(confphgdf,'skap',1, conf_ave)
        print('Phantom group results collected for conformation ')
        print(confphgdf.iloc[0:51])
        print(confphgdf.info())
        if seperate_files == 1:
            print('PHG conformation results written to seperate file')
            confphgdf.to_csv(confebvphg, index=False, header=False, sep=' ')
            print(f'PHG conformation results written to {confebvphg}')
    else:
        print('Phantom group results collected for conformation not collected')
else:
    print('Conformation results not collected')

#Reading rankorder SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if rankorder_results == 1:

    #Reading sol files
    rankordersol = readingfilefwf(rankordersolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    rankordersolids, rankordersolph = solread(rankordersol)
    #Seperating solutions by traits
    rankorderdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    rankorderdf = solutions(rankordersolids, rankorderdf, 1, 'mjaltarod','id')
    rankorderdf = solutions(rankordersolids, rankorderdf, 2, 'gaedarod','id')
    #Reading own observations files
    ownobs_rankorder = ownobs(rankorderobs,rankorder_columns)
    #Creating average groups to scale
    rankorder_ave = avegroup(rankorderdf,ownobs_rankorder)
    #Scaling results
    rankorderdf['mjaltarod'] = scaling(rankorderdf,'mjaltarod',-1, rankorder_ave)
    rankorderdf['gaedarod'] = scaling(rankorderdf,'gaedarod',-1, rankorder_ave)

    #Counting daughters with own obs
    ownobs_rankorder = pd.merge(left=ownobs_rankorder, right=ped[['id','sire']], on='id', how='left')

    rankorderdf = countingoff(ownobs_rankorder,'mjaltarod','offrankorder', rankorderdf)

    if seperate_files == 1:
        print('Rankorder results written to seperate file')
        rank_results = rankorderdf[rankorderebv_columns].astype(float).round().astype(int)
        rank_results.to_csv(rankorderebv, index=False, header=False, sep=' ')
        print(f'Rankorder results written to {rankorderebv}')

    print(rankorderdf.iloc[50000:50015])
    print(rankorderdf.info())
    print(rankorderdf)

    if phantomcollection == 1:

        rankorderphgdf = rankordersolph['code_id'].copy()
        rankorderphgdf = solutions(rankordersolph, rankorderphgdf, 1, 'mjaltarod','code_id')
        rankorderphgdf = solutions(rankordersolph, rankorderphgdf, 2, 'gaedarod','code_id')
        rankorderphgdf['mjaltarod'] = scaling(rankorderphgdf,'mjaltarod',-1, rankorder_ave)
        rankorderphgdf['gaedarod'] = scaling(rankorderphgdf,'gaedarod',-1, rankorder_ave)
        print('Phantom group results collected for rank order')
        print(rankorderphgdf.iloc[0:51])
        print(rankorderphgdf.info())
        if seperate_files == 1:
            print('PHG Rankorder results written to seperate file')
            rankorderphgdf.to_csv(rankorderebvphg, index=False, header=False, sep=' ')
            print(f'PHG Rankorder results written to {rankorderebvphg}')
    else:
        print('Phantom group results collected for rank order not collected')
else:
    print('Rank order results not collected')


#Reading longevity SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if long_results == 1:

    #Reading sol files
    longsol = readingfilefwf(longsolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    longsolids, longsolph = solread(longsol)
    #Seperating solutions by traits
    longdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    longdf = solutions(longsolids, longdf, 1, 'L1','id')
    longdf = solutions(longsolids, longdf, 2, 'L2','id')
    longdf = solutions(longsolids, longdf, 3, 'L3','id')
    # longdf = longdf.loc[(longdf['id'] < 900000000000000)]
    # unscaledlong = longdf.copy()
    # unscaledlong.columns = ['id', 'L1x', 'L2x', 'L3x']
    #Reading own observations files
    ownobs_long = ownobs(longobs,longobs_columns)
    #Creating average groups to scale
    long_ave = avegroup(longdf,ownobs_long)
    #Scaling results
    longdf['L1'] = scaling(longdf,'L1',1, long_ave)
    longdf['L2'] = scaling(longdf,'L2',1, long_ave)
    longdf['L3'] = scaling(longdf,'L3',1, long_ave)

    #Counting daughters with own obs
    ownobs_long = pd.merge(left=ownobs_long, right=ped[['id','sire']], on='id', how='left')

    longdf = countingoff(ownobs_long,'L1','offlong1', longdf)
    longdf = countingoff(ownobs_long,'L2','offlong2', longdf)
    longdf = countingoff(ownobs_long,'L3','offlong3', longdf)

    if seperate_files == 1:
        print('Longgevity results written to seperate file')
        long_results = longdf[longebv_columns].astype(float).round().astype(int)
        long_results.to_csv(longebv, index=False, header=False, sep=' ')
        print(f'Rankorder results written to {longebv}')

    print(longdf.iloc[50000:50015])
    print(longdf.info())

    if phantomcollection == 1:

        longphgdf = longsolph['code_id'].copy()
        longphgdf = solutions(longsolph, longphgdf, 1, 'L1','code_id')
        longphgdf = solutions(longsolph, longphgdf, 2, 'L2','code_id')
        longphgdf = solutions(longsolph, longphgdf, 3, 'L3','code_id')
        longphgdf['L1'] = scaling(longphgdf,'L1',1, long_ave)
        longphgdf['L2'] = scaling(longphgdf,'L2',1, long_ave)
        longphgdf['L3'] = scaling(longphgdf,'L3',1, long_ave)
        print('Phantom group results collected for longevity')
        print(longphgdf.iloc[0:51])
        print(longphgdf.info())
        if seperate_files == 1:
            print('PHG Rankorder results written to seperate file')
            longphgdf.to_csv(longebvphg, index=False, header=False, sep=' ')
            print(f'PHG Rankorder results written to {longebvphg}')
    else:
        print('Phantom group results collected for longevity not collected')

else:
    print('Longevity results not collected')
#-------------------------------------------------------------
#-------------------------------------------------------------



#-------------------------------------------------------------
#-------------------------------------------------------------
#Program reads imported datafiles with EBVs for yield, scs, persistancy and
    #longevity along with accuracy for yield and scs.
    #animals that don't have EBVs for yield are removed.
#Program then collects results for feritlity, conformation and rank order, either
    #from seperate files or from program itself.
#-------------------------------------------------------------
#-------------------------------------------------------------
if collectresults == 1:
    results = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results

    if yield_results == 0:
        results = combineresultscsv(results,yieldebv,yieldebv_columns)
    elif yield_results == 1:
        results = combineresultsdf(results,yielddf,yieldebv_columns)

    if fertility_results == 0:
        results = combineresultscsv(results,fertilityebv,fertilityebv_columns)
    elif fertility_results == 1:
        results = combineresultsdf(results,fertilitydf,fertilityebv_columns)

    if conf_results == 0:
        results = combineresultscsv(results,confebv,confebv_columns)
    elif conf_results == 1:
        results = combineresultsdf(results,confdf,confebv_columns)

    if rankorder_results == 0:
        results = combineresultscsv(results,rankorderebv,rankorderebv_columns)
    elif rankorder_results == 1:
        results = combineresultsdf(results,rankorderdf,rankorderebv_columns)

    if long_results == 0:
        results = combineresultscsv(results,longebv,longebv_columns)
    elif long_results == 1:
        results = combineresultsdf(results,longdf,longebv_columns)

    results = combineresultsfwf(results,accyield,accyield_columns,widths_accyield) #accuracy yield
    results = combineresultsfwf(results,accscs,accscs_columns,widths_accscs) #accuracy scs

    results = results[(results['my1'].notnull().astype(int) == 1)] #only results for animals in yield file

    #Combination of three lactations to make one trait
    results['myt'] = (results['my1'] * 0.5 +
                            results['my2'] * 0.3 +
                            results['my3'] * 0.2 )

    results['fyt'] = (results['fy1'] * 0.5 +
                            results['fy2'] * 0.3 +
                            results['fy3'] * 0.2 )

    results['pyt'] = (results['py1'] * 0.5 +
                            results['py2'] * 0.3 +
                            results['py3'] * 0.2 )

    results['fpt'] = (results['fp1'] * 0.5 +
                            results['fp2'] * 0.3 +
                            results['fp3'] * 0.2 )

    results['ppt'] = (results['pp1'] * 0.5 +
                            results['pp2'] * 0.3 +
                            results['pp3'] * 0.2 )

    results['yieldtotal'] = (results['fyt'] * 0.47 +
                            results['pyt'] * 0.48 +
                            results['ppt'] * 0.05 )

    results['scst'] = (results['scs1'] * 0.5 +
                            results['scs2'] * 0.3 +
                            results['scs2'] * 0.2 )

#Creation of a total grade for udder
    results['jugur'] = (results['jugurfesta'] * 0.35 +
                            results['jugurband'] * 0.15 +
                            results['jugurdypt'] * 0.5 )

#Creation of a total grade for tits
    results['spenar'] = (results['spenalengd'] * 0.3 +
                            (200 - results['spenathykkt']) * 0.3 +
                            results['spenastada'] * 0.4 )

#Creation of a total grade for milking
    results['mjaltir_t'] = (results['mjaltir'] * 0.6 +
                            results['mjaltarod'] * 0.4 )

#Creation of a total grade for bulls with grade for longevity
    results.loc[(results['offlong1'] >= 20), 'total'] = (
                    results['yieldtotal']*0.36 +
                    results['fertility']*0.10 +
                    results['scst']*0.08 +
                    results['jugur']*0.10 +
                    results['spenar']*0.10 +
                    results['mjaltir_t']*0.08 +
                    results['skap']*0.08 +
                    results['L3']*0.10)

#Creation of a total grade for cows
    results.loc[(results['offlong1'] < 20), 'total'] = (
                    results['yieldtotal']*0.36 +
                    results['fertility']*0.11 +
                    results['scst']*0.09 +
                    results['jugur']*0.11 +
                    results['spenar']*0.13 +
                    results['mjaltir_t']*0.10 +
                    results['skap']*0.10 +
                    results['L3']*0.0)

    print(x.iloc[500:515])
    print(x.info())

    print(results.iloc[500000:500015])
    print(results.info())
#Scaled EBVs written to disc
    if writebranda == 1:
        results.loc[:, 'skap2'] = results['skap']
        branda = results[brandafile_columns].astype(int)  #72 columns

        np.savetxt(brandafile, branda,
        fmt='%15s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s')


        print(branda.iloc[500000:500015])
        print(branda.info())

if collectdmufiles == 1:
    #---------------------------------------------------------------------------
    #This function that collects dmu1 and dmu5 files into one directory
    #---------------------------------------------------------------------------
    def collectfiles(trait,year):
        path = f'../DMU/{year}/dmufiles' #creation of dir dmufiles
        isExist = os.path.exists(path) #first program checks if dir already exits
        if not isExist:
            os.makedirs(path);
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')

        path = f'../DMU/{year}/dmufiles/dmu1_{trait}' #copies dmu1 from trait file
        isExist = os.path.exists(path)                #and renames with trait
        if not isExist:                               #checks if file exists already
            shutil.copy(f'../DMU/{year}/{trait}/dmu1', path)
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')

        path = f'../DMU/{year}/dmufiles/dmu5_{trait}' #same as above, dmu5 files copied
        isExist = os.path.exists(path)
        if not isExist:
            shutil.copy(f'../DMU/{year}/{trait}/dmu5', path)
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')


    collectfiles('fertility',yearmonth)
    collectfiles('conformation',yearmonth)
    collectfiles('rankorder',yearmonth)


if plotting == 1:
    #---------------------------------------------------------------------------
    #If plotting == 1 then program will create plots with genetic trends
    #Functions that create plots are in kynbotamat_module
    #---------------------------------------------------------------------------
    branda['BY'] = (branda.id.astype(str).str[:4]).astype(int)
    branda2000 = branda.loc[(branda['BY'] >= 2000 )]

    #Mean EBV per birth year is found for all traits
    argangar = (branda2000.groupby('BY')['my1','my2','my3','fy1','fy2','fy3',
    'py1','py2','py3','fp1','fp2','fp3',
    'pp1','pp2','pp3','scs1','scs2','scs3',
    'milkper', 'fatper', 'protper',
    'fer_lact1','fer_lact2','fer_lact3','CR0','ICF','IFL',
    'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap','mjaltarod', 'gaedarod','longevity',
    'myt','fyt','pyt','fpt','ppt',
    'yieldtotal','fertility','scs','jugur','spenar','mjaltir_t','skap2','total'
        ].mean()).reset_index()
    argangar.to_csv('../results/meanargangar.txt' , index=False, header=True, sep=' ')

    xticks = [2000,2002,2004,2006,2008,2010,2012,2014,2016,2018,2020]
    yticks = [75,80,85,90,95,100,105,110]

    #Creating figure and 16 subplots
    fig, ((ax1, ax2, ax3, ax4),
        (ax5, ax6, ax7, ax8),
        ( ax9, ax10 ,ax11, ax12),
        (ax13, ax14, ax15, ax16))  = plt.subplots(4,4, sharey=True, sharex=True)

    #This figure shows regression
    plottingmeansns(ax1,argangar,'BY','yieldtotal','Skalaðar einkunnir','Afurðir heildareinkun')
    plottingmeansns(ax2,argangar,'BY','myt','','Mjólk kg')
    plottingmeansns(ax3,argangar,'BY','fyt','','Fita kg')
    plottingmeansns(ax4,argangar,'BY','pyt','','Prótein kg')

    plottingmeansns(ax5,argangar,'BY','fpt','Skalaðar einkunnir','Fitu %')
    plottingmeansns(ax6,argangar,'BY','ppt','','Prótein %')
    plottingmeansns(ax7,argangar,'BY','milkper','','Mjólkurúthald')
    plottingmeansns(ax8,argangar,'BY','scs','','Frumutala')

    plottingmeansns(ax9,argangar,'BY','CR0','Skalaðar einkunnir','Kvígur, fanghlutfall við fyrstu sæðingu')
    plottingmeansns(ax10,argangar,'BY','ICF','','Bil milli burðar og fyrstu sæðingar')
    plottingmeansns(ax11,argangar,'BY','IFL','','Bil milli fyrstu og seinustu sæðingar')
    plottingmeansns(ax12,argangar,'BY','boldypt','','Boldýpt')

    plottingmeansns(ax13,argangar,'BY','utlogur','Skalaðar einkunnir','Útlögur')
    plottingmeansns(ax14,argangar,'BY','jugur','','Júgur')
    plottingmeansns(ax15,argangar,'BY','spenar','','Spenar')
    plottingmeansns(ax16,argangar,'BY','total','','Heildareinkunn')

    fig.suptitle('Aðhvarf á meðalkynbótamat árganga', fontsize=26, fontweight ="bold")
    plt.subplots_adjust(left=0.07, bottom=0.08, right=0.96, top=None, wspace=0.05, hspace=0.11)
    fig.set_size_inches([20, 10])
    plt.savefig('../figures/regressionbirthyear20211127.png')

    #Creating figure and 16 subplots
    fig, ((ax1, ax2, ax3, ax4),
        (ax5, ax6, ax7, ax8),
        ( ax9, ax10 ,ax11, ax12),
        (ax13, ax14, ax15, ax16))  = plt.subplots(4,4, sharey=True, sharex=True)

    plottingmean(ax1,branda2000,'my1','Skalaðar einkunnir','Mjólk kg 1. mjalt')
    plottingmean(ax1,branda2000,'my2','Skalaðar einkunnir','Mjólk kg 2. mjalt')
    plottingmean(ax1,branda2000,'my3','Skalaðar einkunnir','Mjólk kg 3. mjalt')
    plottingmean(ax2,branda2000,'fy1','','Fita kg 1. mjalt')
    plottingmean(ax2,branda2000,'fy2','','Fita kg 2. mjalt')
    plottingmean(ax2,branda2000,'fy3','','Fita kg 3. mjalt')
    plottingmean(ax3,branda2000,'py1','','Prótein kg 1. mjalt')
    plottingmean(ax3,branda2000,'py2','','Prótein kg 2. mjalt')
    plottingmean(ax3,branda2000,'py3','','Prótein kg 3. mjalt')
    plottingmean(ax4,branda2000,'yieldtotal','','Heildarafurðaeinkun')

    plottingmean(ax5,branda2000,'fp1','Skalaðar einkunnir','Fitu % 1. mjalt')
    plottingmean(ax5,branda2000,'fp2','Skalaðar einkunnir','Fitu % 2. mjalt')
    plottingmean(ax5,branda2000,'fp3','Skalaðar einkunnir','Fitu % 3. mjalt')
    plottingmean(ax6,branda2000,'pp1','','Prótein % 1. mjalt')
    plottingmean(ax6,branda2000,'pp2','','Prótein % 2. mjalt')
    plottingmean(ax6,branda2000,'pp3','','Prótein % 3. mjalt')
    plottingmean(ax7,branda2000,'scs1','','Frumutala 1. mjalt')
    plottingmean(ax7,branda2000,'scs2','','Frumutala 2. mjalt')
    plottingmean(ax7,branda2000,'scs3','','Frumutala 3. mjalt')
    plottingmean(ax8,branda2000,'milkper','','Mjólkurúthald')
    plottingmean(ax8,branda2000,'fatper','','Fituúthald')
    plottingmean(ax8,branda2000,'protper','','Próteinúthald')

    plottingmean(ax9,branda2000,'CR0','Skalaðar einkunnir','Kvígur, fanghlutfall við fyrstu sæðingu')
    plottingmean(ax9,branda2000,'ICF','Skalaðar einkunnir','Bil milli burðar og fyrstu sæðingar')
    plottingmean(ax9,branda2000,'IFL','Skalaðar einkunnir','Bil milli fyrstu og seinustu sæðingar')
    plottingmean(ax10,branda2000,'boldypt','','Boldýpt')
    plottingmean(ax10,branda2000,'utlogur','','Útlögur')
    plottingmean(ax10,branda2000,'yfirlina','','Yfirlína')
    plottingmean(ax11,branda2000,'malabreidd','','Malabreidd')
    plottingmean(ax11,branda2000,'malahalli','','Malahalli')
    plottingmean(ax11,branda2000,'malabratti','','Malabratti')
    plottingmean(ax12,branda2000,'stada_haekla_hlid','','Staða hækla - hlið')
    plottingmean(ax12,branda2000,'stada_haekla_aftan','','Staða hækla - aftan')
    plottingmean(ax12,branda2000,'klaufhalli','','Klaufhalli')

    plottingmean(ax13,branda2000,'jugurfesta','Skalaðar einkunnir','Júgurfesta')
    plottingmean(ax13,branda2000,'jugurband','Skalaðar einkunnir','Júgurband')
    plottingmean(ax13,branda2000,'jugurdypt','Skalaðar einkunnir','Júgurdýpt')
    plottingmean(ax14,branda2000,'spenalengd','','Spenalengd')
    plottingmean(ax14,branda2000,'spenathykkt','','Spenaþykkt')
    plottingmean(ax14,branda2000,'spenastada','','Spenastaða')
    plottingmean(ax15,branda2000,'mjaltir','','Mjaltir')
    plottingmean(ax15,branda2000,'skap','','Skap')
    plottingmean(ax16,branda2000,'total','','Heildareinkunn')

    fig.suptitle('Meðalkynbótamat árganga', fontsize=26, fontweight ="bold")
    plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

    fig.set_size_inches([20, 10])
    plt.savefig('../figures/meanbybirthyear20211127.png')

    # plt.show()




#
