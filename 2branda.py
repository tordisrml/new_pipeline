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

#Option to read SOL files and collect results, 0 means skip
# and 1 means SOL files will be read
yield_collect = 1
plotyieldPE = 1

fertility_collect = 1
conf_collect = 1
rankorder_collect = 1
long_collect = 1

#An option to collect phantom parent group results to seperate files.
#Does not collect for yield only other traits. Should be 0 by defult
phantomcollection = 0

#An option to collect and plot unscaled EBVs for the traits. Should be 1 by default.
collectunscaled = 0
plotunscaled = 0

#Option to write scaled results for above traits to seperate files to be read
#later, 0 means skip and 1 means seperate result files will be written to disc
seperate_files = 0

#Option to collect results for a large datafile to be read by Huppa.
#--IF ABOVE OPTIONS ARE SET TO 0 --> program will collect results from seperate files!
#--IF ABOVE OPTIONS ARE SET TO 1 --> program will collect results straight from the program itself.
collectresults = 1    #0 means don't, 1 means collect

#Option to write large datafile to disc. 1 by defult
writebranda = 1 #(only when collectresults = 1 )
#Option to write seperate excel file only with nautastöðvar bulls (only when writebranda = 1 )
nautastod = 0
#Option to plot branda file
plottingbranda = 0

#Option to show plots on computer screen. User must have launched x-ming
# Should be 0 by default.
pltshow = 0
#---------------------------------------------------------------------------
# PART 1 - Defining variables and file names
#---------------------------------------------------------------------------
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
yearmonth = control.loc[0,'control']
brandafile = f'../{yearmonth}/results/branda{yearmonth}'
#Scaling year, present year - 5 years
scalingyear = pd.to_numeric(control.loc[13,'control'])
#Scaling objects
imean = pd.to_numeric(control.loc[14,'control'])
isd = pd.to_numeric(control.loc[15,'control'])
#---------------------------------------------------------------------------
#Radnrkodifile
radnrkodifile = f'../{yearmonth}/dmu_data/radnrkodi' #Created by prep_tdm.f
#Format of radnrkodi file created by prep_tdm.f
radnrkodi_columns = ['id','code_id','stada','norec','fix1','fix2', 'fix3','sex']
widths_radnrkodi = [15,9,3,3,6,6,6,2]
#---------------------------------------------------------------------------
#Name of pedigree file
pedigreefile = control.loc[12,'control']
#Format of pedigree file from Huppa
ped_columns = ['id','dam','sire','unused','sex','unused2', 'bullno', 'name', 'farm']
ped_widths = [15,15,15,12,1,2,5,20,20]
#---------------------------------------------------------------------------

#DMU sol files
mysol = f'../{yearmonth}/DMU/my/SOL'
fysol = f'../{yearmonth}/DMU/fy/SOL'
pysol = f'../{yearmonth}/DMU/py/SOL'
scssol = f'../{yearmonth}/DMU/scs/SOL'
fpsol = f'../{yearmonth}/DMU/fp/SOL'
ppsol = f'../{yearmonth}/DMU/pp/SOL'

fertilitysolfile = f'../{yearmonth}/DMU/fer/SOL'
rankordersolfile = f'../{yearmonth}/DMU/rank/SOL'
confsolfile = f'../{yearmonth}/DMU/conf/SOL'
longsolfile = f'../{yearmonth}/DMU/long/SOL'
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
fertilityobs = f'../{yearmonth}/dmu_data/dmu_fertility.txt'
fertility_columns = ['code_id','HBY','HC1','HC2','HC3','IYM0','IYM1','IYM2',
    'IYM3','AGEi_h','AGEc_1','AGEc_2','AGEc_3','tech_h',
    'CR0','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']
#---------------------------------------------------------------------------
confobs = f'../{yearmonth}/dmu_data/dmu_conformation.txt'
conf_columns = ['code_id','HdomsY','lact','AGEc_1',
    'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap']
#---------------------------------------------------------------------------
rankorderobs = f'../{yearmonth}/dmu_data/dmu_rankorder.txt'
rankorder_columns = ['code_id', 'year', 'mjaltarod', 'gaedarod']
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
longobs = f'../{yearmonth}/dmu_data/dmu_long.txt'
longobs_columns = ['code_id','AGEc_1','CYM1','h5y','herdCY1','L1','L2','L3']
#---------------------------------------------------------------------------

# #Seperate EBV result files
yieldebv = f'../{yearmonth}/results/yieldebv.txt'
fertilityebv = f'../{yearmonth}/results/fertilityebv.txt' #written by this program if option above set to 1
confebv = f'../{yearmonth}/results/conformationebv.txt' #written by this program if option above set to 1
rankorderebv = f'../{yearmonth}/results/rankorderebv.txt' #written by this program if option above set to 1
longebv = f'../{yearmonth}/results/longebv.txt'        #written by this program if option above set to 1
#Accuracy files for yield and scs
accyield = f'../{yearmonth}/results/accuracy.sol' #From programs by JHE
accscs = f'../{yearmonth}/results/accuracy_f.sol' #From programs by JHE

#---------------------------------------------------------------------------
#Columns in EBV files
#---------------------------------------------------------------------------
accyield_columns = ['code_id','ownrec','nopar','offt','offn','offyield','yield_acc','SE']
widths_accyield = [6,4,5,7,6,6,7,7]

accscs_columns =  ['code_id','ownrec','nopar','offt','offn','offscs','scs_acc','SE']
widths_accscs = [6,4,5,7,6,6,7,7]
#---------------------------------------------------------------------------
#Space seperated files/dataframes
#(Created by this program if option above set to 1)
yieldebv_columns = ['id','my1','my2','my3','fy1','fy2','fy3','py1','py2','py3','fp1','fp2','fp3',
    'pp1','pp2','pp3','scs1','scs2','scs3', 'milkper', 'fatper', 'protper' ]

fertilityebv_columns = ['id','fer_lact1','fer_lact2','fer_lact3','CR0','ICF',
    'IFL','fertility','offCR0','offICF1','offICF2','offICF3']

confebv_columns = ['id','boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap','offconf']

rankorderebv_columns = ['id','mjaltarod', 'gaedarod','offrankorder']

longebv_columns = ['id','L1','L2','L3','offlong1','offlong2','offlong3']


#Phantom group parent results
fertilityebvphg = f'../{yearmonth}/results/fertilityebvphg.txt' #written by this program if option above set to 1
confebvphg = f'../{yearmonth}/results/conformationebvphg.txt' #written by this program if option above set to 1
rankorderebvphg = f'../{yearmonth}/results/rankorderebvphg.txt' #written by this program if option above set to 1
longebvphg = f'../{yearmonth}/results/longebvphg.txt' #written by this program if option above set to 1
#---------------------------------------------------------------------------

#Large file with all solutions scaled
brandafile_columns = ['id',                                             #1
        'my1','my2','my3','fy1','fy2','fy3', #6
        'py1','py2','py3','fp1','fp2','fp3',    #6
        'pp1','pp2','pp3','scs1','scs2','scs3',             #6
        'milkper', 'fatper', 'protper',                                 #3
        'fer_lact1','fer_lact2','fer_lact3','CR0','ICF','IFL',          #6
        'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti', #6
        'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta', #4
        'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',  #5
        'mjaltir', 'skap','mjaltarod', 'gaedarod','L1','L2','L3',              #7
        'myt','fyt','pyt','fpt','ppt',                 #5
        'yieldtotal','fertility','scst','skrokkur','jugur','spenar','mjaltir_t', #7
        'total',                                                             #1
        'offyield', 'yield_acc', 'offscs', 'scs_acc', #4
        'offconf','offrankorder','offCR0','offICF1','offICF2','offICF3', #6
        'offlong1','offlong2','offlong3'                                        #3
        ] #76 columns
widths_branda = [15,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4, #
    4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4] #
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

#Funtion to write out seperate result files
#-------------------------------------------------------------
def seperateresults(df,columns,ebvfile,trait):
        print(f'{trait} results written to seperate file')
        dfresults = df[columns].astype(float).round(2)#.astype(int)
        dfresults[['id']] = dfresults[['id']].astype(int)
        dfresults.to_csv(ebvfile, index=False, header=False, sep=' ')
        print(f'{trait} results written to {confebv}')

        print(df.iloc[50000:50015])
        print(df.info())
#-------------------------------------------------------------
#Function to read imported result files and combine to a big one
#Fixed with files
def combineresultsfwf(df, file, columns,widths_file,t1,t2):
    print( f'-------------------------------------------' )
    print( f'Reading {file} file .....' )
    print( f'-------------------------------------------' )
    df2 = pd.read_fwf(
        file,
        widths=widths_file,
        header=None,
        names=columns)
    #Merging ownobs and id file to bring back einstaklingsnumer used in huppa
    df2 = pd.merge(left=df2, right=radnrkodi[['id','code_id']], on='code_id',  how='left').fillna(0, downcast='infer')
    df = pd.merge(left=df, right=df2[['id',t1,t2]], on='id',  how='left')
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

#---------------------------------------------------------------------------
# PART 2 - Collection of SOL results
#---------------------------------------------------------------------------
#-------------------------------------------------------------
#-------------------------------------------------------------
#START OF PROGRAM!!!
#-------------------------------------------------------------
#-------------------------------------------------------------
#Reading in id codes to replace in SOL files
if ((collectresults == 1)
    | (fertility_collect == 1)
    | (conf_collect == 1)
    | (rankorder_collect == 1)
    | (long_collect == 1)
    | (yield_collect == 1)
    ):
    radnrkodi = readingfilefwf(radnrkodifile,radnrkodi_columns,widths_radnrkodi)
    radnrkodi = radnrkodi.drop(
        ['norec','fix1','fix2', 'fix3','sex'], axis = 1)

    #Reading pedigree
    ped = readingfilefwf(pedigreefile,ped_columns,ped_widths)

#Reading yield SOL files, collecting and scaling results
#Write results to disc if option above set to 1
if (yield_collect == 1):

    myebv = yieldcollect(mysol,'my1','my2','my3','code_id','myper1','myper2','myper3','my1PE','my2PE','my3PE')
    yielddf = pd.merge(left=radnrkodi, right=myebv, on='code_id', how='left')

    fyebv = yieldcollect(fysol,'fy1','fy2','fy3','code_id','fyper1','fyper2','fyper3','fy1PE','fy2PE','fy3PE')
    yielddf = pd.merge(left=yielddf, right=fyebv, on='code_id', how='left')

    pyebv = yieldcollect(pysol,'py1','py2','py3','code_id','pyper1','pyper2','pyper3','py1PE','py2PE','py3PE')
    yielddf = pd.merge(left=yielddf, right=pyebv, on='code_id', how='left')

    scsebv = yieldcollect(scssol,'scs1','scs2','scs3','code_id','','','','scs1PE','scs2PE','scs3PE')
    yielddf = pd.merge(left=yielddf, right=scsebv, on='code_id', how='left')

    fpebv = yieldcollect(fpsol,'fp1','fp2','fp3','code_id','','','','fp1PE','fp2PE','fp3PE')
    yielddf = pd.merge(left=yielddf, right=fpebv, on='code_id', how='left')

    ppebv = yieldcollect(ppsol,'pp1','pp2','pp3','code_id','','','','pp1PE','pp2PE','pp3PE')
    yielddf = pd.merge(left=yielddf, right=ppebv, on='code_id', how='left')

    if collectunscaled == 1:
        unscaledyield = yielddf[['id','my1','my2','my3','fy1','fy2','fy3',
            'py1','py2','py3','scs1','scs2','scs3','fp1','fp2','fp3','pp1','pp2','pp3']]

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

    yieldPE = yielddf[['id','my1PE','my2PE','my3PE','fy1PE',
        'fy2PE','fy3PE','py1PE','py2PE','py3PE','scs1PE','scs2PE','scs3PE',
        'fp1PE','fp2PE','fp3PE','pp1PE','pp2PE','pp3PE']]

    if plotyieldPE == 1:
        yieldPE['BY'] = (yieldPE.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6))  = plt.subplots(3,2, sharex=True)

        plottingmean(ax1,yieldPE,'my1PE','PE','Mjólk kg 1. mjalt')
        plottingmean(ax1,yieldPE,'my2PE','PE','Mjólk kg 2. mjalt')
        plottingmean(ax1,yieldPE,'my3PE','PE','Mjólk kg 3. mjalt')

        plottingmean(ax2,yieldPE,'fy1PE','','Fita kg 1. mjalt')
        plottingmean(ax2,yieldPE,'fy2PE','','Fita kg 2. mjalt')
        plottingmean(ax2,yieldPE,'fy3PE','','Fita kg 3. mjalt')

        plottingmean(ax3,yieldPE,'py1PE','PE','Prótein kg 1. mjalt')
        plottingmean(ax3,yieldPE,'py2PE','PE','Prótein kg 2. mjalt')
        plottingmean(ax3,yieldPE,'py3PE','PE','Prótein kg 3. mjalt')

        plottingmean(ax4,yieldPE,'scs1PE','','Frumutala 1. mjalt')
        plottingmean(ax4,yieldPE,'scs2PE','','Frumutala 2. mjalt')
        plottingmean(ax4,yieldPE,'scs3PE','','Frumutala 3. mjalt')

        plottingmean(ax5,yieldPE,'fp1PE','PE','Fita % 1. mjalt')
        plottingmean(ax5,yieldPE,'fp2PE','PE','Fita % 2. mjalt')
        plottingmean(ax5,yieldPE,'fp3PE','PE','Fita % 3. mjalt')

        plottingmean(ax6,yieldPE,'pp1PE','','Prótein % 1. mjalt')
        plottingmean(ax6,yieldPE,'pp2PE','','Prótein % 2. mjalt')
        plottingmean(ax6,yieldPE,'pp3PE','','Prótein % 3. mjalt')

        fig.suptitle(f'Meðal PE árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([18, 9])
        plt.savefig(f'../{yearmonth}/figures/PEyieldmeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if plotunscaled == 1:
        unscaledyield['BY'] = (unscaledyield.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6))  = plt.subplots(3,2, sharex=True)

        plottingmean(ax1,unscaledyield,'my1','Óskalað DMU kynbótamat','Mjólk kg 1. mjalt')
        plottingmean(ax1,unscaledyield,'my2','Óskalað DMU kynbótamat','Mjólk kg 2. mjalt')
        plottingmean(ax1,unscaledyield,'my3','Óskalað DMU kynbótamat','Mjólk kg 3. mjalt')

        plottingmean(ax2,unscaledyield,'fy1','','Fita kg 1. mjalt')
        plottingmean(ax2,unscaledyield,'fy2','','Fita kg 2. mjalt')
        plottingmean(ax2,unscaledyield,'fy3','','Fita kg 3. mjalt')

        plottingmean(ax3,unscaledyield,'py1','Óskalað DMU kynbótamat','Prótein kg 1. mjalt')
        plottingmean(ax3,unscaledyield,'py2','Óskalað DMU kynbótamat','Prótein kg 2. mjalt')
        plottingmean(ax3,unscaledyield,'py3','Óskalað DMU kynbótamat','Prótein kg 3. mjalt')

        plottingmean(ax4,unscaledyield,'scs1','','Frumutala 1. mjalt')
        plottingmean(ax4,unscaledyield,'scs2','','Frumutala 2. mjalt')
        plottingmean(ax4,unscaledyield,'scs3','','Frumutala 3. mjalt')

        plottingmean(ax5,unscaledyield,'fp1','Óskalað DMU kynbótamat','Fita % 1. mjalt')
        plottingmean(ax5,unscaledyield,'fp2','Óskalað DMU kynbótamat','Fita % 2. mjalt')
        plottingmean(ax5,unscaledyield,'fp3','Óskalað DMU kynbótamat','Fita % 3. mjalt')

        plottingmean(ax6,unscaledyield,'pp1','','Prótein % 1. mjalt')
        plottingmean(ax6,unscaledyield,'pp2','','Prótein % 2. mjalt')
        plottingmean(ax6,unscaledyield,'pp3','','Prótein % 3. mjalt')

        fig.suptitle(f'Meðal óskalað kynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([18, 9])
        plt.savefig(f'../{yearmonth}/figures/unscaledyieldmeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if seperate_files == 1:
        seperateresults(yielddf,yieldebv_columns,yieldebv,'Yield')

else:
    print('Rank order results not collected')

#Reading fertilty SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if (fertility_collect == 1):
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

    if collectunscaled == 1:
        unscaledfer = fertilitydf[['id','CR0','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']]

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

    if plotunscaled == 1:
        unscaledfer['BY'] = (unscaledfer.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1, ax2, ax3))  = plt.subplots(3, sharex=True)

        plottingmean(ax1,unscaledfer,'CR0','Óskalað DMU kynbótamat','Fanghlutfall 1. sæðingar, kvígur')
        plottingmean(ax2,unscaledfer,'ICF1','Óskalað DMU kynbótamat','Bil burður 1. sæð. - 1. mjalt')
        plottingmean(ax2,unscaledfer,'ICF2','Óskalað DMU kynbótamat','Bil burður 1. sæð. - 2. mjalt')
        plottingmean(ax2,unscaledfer,'ICF3','Óskalað DMU kynbótamat','Bil burður 1. sæð. - 3. mjalt')

        plottingmean(ax3,unscaledfer,'IFL1','Óskalað DMU kynbótamat','Bil 1. sæð. til síðasta sæð. - 1. mjalt')
        plottingmean(ax3,unscaledfer,'IFL2','Óskalað DMU kynbótamat','Bil 1. sæð. til síðasta sæð. - 2. mjalt')
        plottingmean(ax3,unscaledfer,'IFL3','Óskalað DMU kynbótamat','Bil 1. sæð. til síðasta sæð. - 3. mjalt')

        fig.suptitle(f'Meðal óskalað kynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([18, 9])
        plt.savefig(f'../{yearmonth}/figures/unscaledfertilitymeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if seperate_files == 1:
        seperateresults(fertilitydf,fertilityebv_columns,fertilityebv,'Fertility')

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
if (conf_collect == 1):

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
    confdf = solutions(confsolids, confdf, 7, 'stada_hh','id')
    confdf = solutions(confsolids, confdf, 8, 'stada_ha','id')
    confdf = solutions(confsolids, confdf, 9, 'klaufhalli','id')
    confdf = solutions(confsolids, confdf, 10, 'jugurfesta','id')
    confdf = solutions(confsolids, confdf, 11, 'jugurband','id')
    confdf = solutions(confsolids, confdf, 12, 'jugurdypt','id')
    confdf = solutions(confsolids, confdf, 13, 'spenalengd','id')
    confdf = solutions(confsolids, confdf, 14, 'spenathykkt','id')
    confdf = solutions(confsolids, confdf, 15, 'spenastada','id')
    confdf = solutions(confsolids, confdf, 16, 'mjaltir','id')
    confdf = solutions(confsolids, confdf, 17, 'skap','id')

    if collectunscaled == 1:
        unscaledconf = confdf[['id','boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
            'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
            'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
            'mjaltir', 'skap']]

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
    confdf['stada_hh'] = scaling(confdf,'stada_hh',1, conf_ave)
    confdf['stada_ha'] = scaling(confdf,'stada_ha',1, conf_ave)
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

    if plotunscaled == 1:
        unscaledconf['BY'] = (unscaledconf.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1, ax2),(ax3, ax4),(ax5, ax6))  = plt.subplots(3,2, sharex=True)

        plottingmean(ax1,unscaledconf,'boldypt','Óskalað DMU kynbótamat','Boldýpt')
        plottingmean(ax1,unscaledconf,'utlogur','Óskalað DMU kynbótamat','Útlögur')
        plottingmean(ax1,unscaledconf,'yfirlina','Óskalað DMU kynbótamat','Yfirlína')
        plottingmean(ax2,unscaledconf,'malabreidd','','Malabreidd')
        plottingmean(ax2,unscaledconf,'malahalli','','Malahalli')
        plottingmean(ax2,unscaledconf,'malabratti','','Malabratti')
        plottingmean(ax3,unscaledconf,'stada_hh','Óskalað DMU kynbótamat','Staða hækla - hlið')
        plottingmean(ax3,unscaledconf,'stada_ha','Óskalað DMU kynbótamat','Staða hækla - aftan')
        plottingmean(ax3,unscaledconf,'klaufhalli','Óskalað DMU kynbótamat','Klaufhalli')
        plottingmean(ax4,unscaledconf,'jugurfesta','','Júgurfesta')
        plottingmean(ax4,unscaledconf,'jugurband','','Júgurband')
        plottingmean(ax4,unscaledconf,'jugurdypt','','Júgurdýpt')
        plottingmean(ax5,unscaledconf,'spenalengd','Óskalað DMU kynbótamat','Spenalengd')
        plottingmean(ax5,unscaledconf,'spenathykkt','Óskalað DMU kynbótamat','Spenaþykkt')
        plottingmean(ax5,unscaledconf,'spenastada','Óskalað DMU kynbótamat','Spenastaða')
        plottingmean(ax6,unscaledconf,'mjaltir','','Mjaltir')
        plottingmean(ax6,unscaledconf,'skap','','Skap')

        fig.suptitle(f'Meðal óskalað kynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([18, 9])
        plt.savefig(f'../{yearmonth}/figures/unscaledconformationmeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if seperate_files == 1:
        seperateresults(confdf,confebv_columns,confebv,'Conformation')

    if phantomcollection == 1:

        confphgdf = confsolph['code_id'].copy()
        confphgdf = solutions(confsolph, confphgdf, 1, 'boldypt','code_id')
        confphgdf = solutions(confsolph, confphgdf, 2, 'utlogur','code_id')
        confphgdf = solutions(confsolph, confphgdf, 3, 'yfirlina','code_id')
        confphgdf = solutions(confsolph, confphgdf, 4, 'malabreidd','code_id')
        confphgdf = solutions(confsolph, confphgdf, 5, 'malahalli','code_id')
        confphgdf = solutions(confsolph, confphgdf, 6, 'malabratti','code_id')
        confphgdf = solutions(confsolph, confphgdf, 7, 'stada_hh','code_id')
        confphgdf = solutions(confsolph, confphgdf, 8, 'stada_ha','code_id')
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
        confphgdf['stada_hh'] = scaling(confphgdf,'stada_hh',1, conf_ave)
        confphgdf['stada_ha'] = scaling(confphgdf,'stada_ha',1, conf_ave)
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
if (rankorder_collect == 1):

    #Reading sol files
    rankordersol = readingfilefwf(rankordersolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    rankordersolids, rankordersolph = solread(rankordersol)
    #Seperating solutions by traits
    rankorderdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    rankorderdf = solutions(rankordersolids, rankorderdf, 1, 'mjaltarod','id')
    rankorderdf = solutions(rankordersolids, rankorderdf, 2, 'gaedarod','id')

    if collectunscaled == 1:
        unscaledrank = rankorderdf[['id','mjaltarod','gaedarod']]

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

    if plotunscaled == 1:
        unscaledrank['BY'] = (unscaledrank.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1, ax2))  = plt.subplots(2, sharex=True)

        plottingmean(ax1,unscaledrank,'mjaltarod','Óskalað DMU kynbótamat','Mjaltaröð')
        plottingmean(ax2,unscaledrank,'gaedarod','Óskalað DMU kynbótamat','Gæðaröð')

        fig.suptitle(f'Meðal óskalað kynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([18, 9])
        plt.savefig(f'../{yearmonth}/figures/unscaledrankorderymeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if seperate_files == 1:
        seperateresults(rankorderdf,rankorderebv_columns,rankorderebv,'Rankorder')

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
if (long_collect == 1):

    #Reading sol files
    longsol = readingfilefwf(longsolfile,solcolumns,sol_widths)
    #Sepererating phantom groups from known animals in SOL file and merging real id's
    longsolids, longsolph = solread(longsol)
    #Seperating solutions by traits
    longdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    longdf = solutions(longsolids, longdf, 1, 'L1','id')
    longdf = solutions(longsolids, longdf, 2, 'L2','id')
    longdf = solutions(longsolids, longdf, 3, 'L3','id')

    if collectunscaled == 1:
        unscaledlong = longdf[['id','L1','L2','L3']]

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

    if plotunscaled == 1:
        unscaledlong['BY'] = (unscaledlong.id.astype(str).str[:4]).astype(int)
        #Creating figure and 16 subplots
        fig, ((ax1))  = plt.subplots(1, sharex=True)

        plottingmean(ax1,unscaledlong,'L1','Óskalað DMU kynbótamat','Ending, dagar 1. burður - lok 1. mjalt.')
        plottingmean(ax1,unscaledlong,'L2','Óskalað DMU kynbótamat','Ending, dagar 1. burður - lok 2. mjalt.')
        plottingmean(ax1,unscaledlong,'L3','Óskalað DMU kynbótamat','Ending, dagar 1. burður - lok 3. mjalt.')

        fig.suptitle(f'Meðal óskalað kynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
        plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

        fig.set_size_inches([12, 10])
        plt.savefig(f'../{yearmonth}/figures/unscaledlongevitymeanbybirthyear{yearmonth}.png')

        if pltshow == 1:
            plt.show()

    if seperate_files == 1:
        seperateresults(longdf,longebv_columns,longebv,'Longevity')

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


#---------------------------------------------------------------------------
# PART 3 - Collection of all results and write to file
#---------------------------------------------------------------------------
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

    if (yield_collect == 0):
        results = combineresultscsv(results,yieldebv,yieldebv_columns)
    elif (yield_collect == 1):
        results = combineresultsdf(results,yielddf,yieldebv_columns)

    if (fertility_collect == 0):
        results = combineresultscsv(results,fertilityebv,fertilityebv_columns)
    elif (fertility_collect == 1):
        results = combineresultsdf(results,fertilitydf,fertilityebv_columns)

    if (conf_collect == 0):
        results = combineresultscsv(results,confebv,confebv_columns)
    elif (conf_collect == 1):
        results = combineresultsdf(results,confdf,confebv_columns)

    if (rankorder_collect == 0):
        results = combineresultscsv(results,rankorderebv,rankorderebv_columns)
    elif (rankorder_collect == 1):
        results = combineresultsdf(results,rankorderdf,rankorderebv_columns)

    if (long_collect == 0):
        results = combineresultscsv(results,longebv,longebv_columns)
    elif (long_collect == 1):
        results = combineresultsdf(results,longdf,longebv_columns)

    results = combineresultsfwf(results,accyield,accyield_columns,widths_accyield,'offyield','yield_acc') #accuracy yield
    results = combineresultsfwf(results,accscs,accscs_columns,widths_accscs,'offscs','scs_acc') #accuracy scs

    results['yield_acc'] = results['yield_acc'] * 100
    results['scs_acc'] = results['scs_acc'] * 100

    results[['offyield', 'yield_acc','offscs', 'scs_acc']
        ] = results[['offyield', 'yield_acc'
        ,'offscs', 'scs_acc']].fillna(0).astype(int)

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

#Creation of a total grade for milking
    results['skrokkur'] = (results['boldypt'] * 0.25 +
                            results['utlogur'] * 0.25 +
                            results['yfirlina'] * 0.2 +
                            results['malabreidd'] * 0.3)

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

    print(results.iloc[500000:500015])
    print(results.info())
#Scaled EBVs written to disc
    if writebranda == 1:
        branda = results[brandafile_columns].fillna(0, downcast='infer').astype(int)

        # branda.to_csv(f'../{yearmonth}/results/branda{yearmonth}newformat', index=False, header=False, sep=';')


        brandaold = branda[brandafile_columns].copy()
        brandaold.loc[:,'bandmal_eldra'] = 0
        brandaold.loc[:,'bolur_eldra'] = 0
        brandaold.loc[:,'malir_eldra'] = 0
        brandaold.loc[:,'fotstada_eldra'] = 0
        brandaold.loc[:,'jugurlag_festa_eldra'] = 0
        brandaold.loc[:,'spenalengd_lag_eldra'] = 0
        brandaold.loc[:,'mjaltir_eldra'] = 0
        brandaold.loc[:,'skap_eldra'] = 0
        brandaold['skap2'] = brandaold['skap']
        brandaold.loc[:,'eigin_afurdir'] = 0
        brandaold.loc[:,'no_daughters_utlit_gamla'] = 0
        brandaold.loc[:,'utlit_gamla_acc'] = 0
        brandaold.loc[:,'utlit_nyja_acc'] = 0
        brandaold.loc[:,'mjaltir_acc'] = 0

        brandaold2 = brandaold[['id',                                                         #1
            'my1','my2','my3','fy1','fy2','fy3','py1','py2','py3','fp1','fp2','fp3',    #12
            'pp1','pp2','pp3','fer_lact1','fer_lact2','fer_lact3','scs1','scs2','scs3',  #9

            'bandmal_eldra', 'bolur_eldra', 'malir_eldra','fotstada_eldra',
            'jugurlag_festa_eldra', 'spenalengd_lag_eldra', 'mjaltir_eldra', 'skap_eldra', #8

            'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
            'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
            'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
            'mjaltir', 'skap','mjaltarod', 'gaedarod','L3','myt','fyt','pyt','fpt','ppt', #25
            'eigin_afurdir','yieldtotal','fertility','scst','skrokkur','jugur','spenar','mjaltir_t',
            'skap2','total', 'offyield', 'yield_acc', 'offscs', 'scs_acc',
            'no_daughters_utlit_gamla','utlit_gamla_acc','offconf','utlit_nyja_acc',
            'offrankorder','mjaltir_acc','offlong3','milkper', 'fatper', 'protper',
            'CR0','ICF','IFL']]                                                                   #27    = 82 dálkur!!!



        np.savetxt(f'../{yearmonth}/results/branda{yearmonth}oldformat', brandaold2,
        fmt='%15s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
%3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s')





#nýtt format en fixed width!
#         np.savetxt(brandafile, branda,
#         fmt='%15s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
# %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
# %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
# %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
# %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s \
# %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s %3s')

        print(branda.iloc[500000:500015])
        print(branda.info())

    #---------------------------------------------------------------------------
    # PART 4 - Extra steps
    #---------------------------------------------------------------------------
        if nautastod == 1:
            # ----------------------------------------------------
            #Counting the number of daughters bulls have in the pedigree file
            # ----------------------------------------------------
            dams = ped.loc[ped['sex'] == 2]
            dams.loc[:,'sire_count'] = dams.groupby('sire')['sire'].transform('count')
            sires = dams[['sire', 'sire_count']]
            sires.columns = ['id', 'daughters']
            sires = sires.drop_duplicates(subset=['id'])

            #Merging with saman dataframe so sires have their daughter count
            brandaped = pd.merge(left=branda, right=sires, on='id', how='outer')
            #Merging with pedigree to see bullno and farm id
            brandaped = pd.merge(left=ped[['id','sex','bullno', 'name', 'farm']], right=brandaped, on='id')

            # ----------------------------------------------------
            #This creates an excel file with bulls from Nautstöðin á Hesti
            # ----------------------------------------------------
            nautastod = brandaped.loc[(brandaped['bullno'] > 0) & (brandaped['farm'] == 'Nautastöðin Hesti')]
            #Sorting the dataframe by birthyear and nautastöðvarnúmer
            nautastod['BY'] = (nautastod.id.astype(str).str[:4]).astype(int)
            nautastod = nautastod.sort_values(by=['BY','bullno'],ascending=False)
            #Correct order of columns for excel file
            # Rename the traits for excel file
            # nautastod.columns = [] ---> possible to rename columns
            # Creating the excel file
            nautastod.to_excel(f'../{yearmonth}/results/nautastodblup{yearmonth}.xlsx', index=False, header=True)

        if plottingbranda == 1:                                       #---------------------------------------------------------------------------
            #---------------------------------------------------------------------------
            #If plotting == 1 then program will create plots with genetic trends
            #Functions that create plots are in kynbotamat_module
            #---------------------------------------------------------------------------
            branda['BY'] = (branda.id.astype(str).str[:4]).astype(int)
            branda2000 = branda.loc[(branda['BY'] >= 2000 ) & (branda['BY'] <= 3000 )]

            #Mean EBV per birth year is found for all traits
            argangar = ((branda2000.groupby('BY')['my1','my2','my3','fy1','fy2','fy3',
            'py1','py2','py3','fp1','fp2','fp3',
            'pp1','pp2','pp3','scs1','scs2','scs3',
            'milkper', 'fatper', 'protper',
            'fer_lact1','fer_lact2','fer_lact3','CR0','ICF','IFL',
            'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
            'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
            'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
            'mjaltir', 'skap','mjaltarod', 'gaedarod','L1','L2','L3',
            'myt','fyt','pyt','fpt','ppt',
            'yieldtotal','fertility','scst','jugur','spenar','mjaltir_t','total'
                ].mean()).reset_index()).astype(float).round(2).astype(int)
            argangar.to_csv(f'../{yearmonth}/results/meanargangar.txt' , index=False, header=True, sep=' ')

            #Creating figure and 16 subplots
            fig, ((ax1, ax2, ax3, ax4),
                (ax5, ax6, ax7, ax8),
                ( ax9, ax10 ,ax11, ax12),
                (ax13, ax14, ax15, ax16),
                (ax17, ax18, ax19, ax20))  = plt.subplots(5,4, sharey=True, sharex=True)

            #This figure shows regression
            plottingmeansns(ax1,argangar,'BY','yieldtotal','Skalaðar einkunnir','Afurðir heildareinkun')
            plottingmeansns(ax2,argangar,'BY','myt','','Mjólk kg')
            plottingmeansns(ax3,argangar,'BY','fyt','','Fita kg')
            plottingmeansns(ax4,argangar,'BY','pyt','','Prótein kg')

            plottingmeansns(ax5,argangar,'BY','fpt','Skalaðar einkunnir','Fitu %')
            plottingmeansns(ax6,argangar,'BY','ppt','','Prótein %')
            plottingmeansns(ax7,argangar,'BY','milkper','','Mjólkurúthald')
            plottingmeansns(ax8,argangar,'BY','scst','','Frumutala')

            plottingmeansns(ax9,argangar,'BY','CR0','Skalaðar einkunnir','Kvígur, fanghlutfall við fyrstu sæðingu')
            plottingmeansns(ax10,argangar,'BY','ICF','','Bil milli burðar og fyrstu sæðingar')
            plottingmeansns(ax11,argangar,'BY','IFL','','Bil milli fyrstu og seinustu sæðingar')
            plottingmeansns(ax12,argangar,'BY','L3','','Ending, 1. burður til loka 3. mjaltaskeiðs')

            plottingmeansns(ax13,argangar,'BY','malabreidd','Skalaðar einkunnir','Malarbreidd')
            plottingmeansns(ax14,argangar,'BY','boldypt','','Boldýpt')
            plottingmeansns(ax15,argangar,'BY','utlogur','Skalaðar einkunnir','Útlögur')
            plottingmeansns(ax16,argangar,'BY','skap','','Skap')

            plottingmeansns(ax17,argangar,'BY','mjaltir_t','Skalaðar einkunnir','Mjaltir, heildareinkunn')
            plottingmeansns(ax18,argangar,'BY','jugur','','Júgur')
            plottingmeansns(ax19,argangar,'BY','spenar','','Spenar')
            plottingmeansns(ax20,argangar,'BY','total','','Heildareinkunn')

            fig.suptitle(f'Aðhvarf á meðalkynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
            plt.subplots_adjust(left=0.07, bottom=0.08, right=0.96, top=None, wspace=0.05, hspace=0.11)
            fig.set_size_inches([18, 9])
            plt.savefig(f'../{yearmonth}/figures/regressionbirthyear{yearmonth}.png')

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
            plottingmean(ax12,branda2000,'stada_hh','','Staða hækla - hlið')
            plottingmean(ax12,branda2000,'stada_ha','','Staða hækla - aftan')
            plottingmean(ax12,branda2000,'klaufhalli','','Klaufhalli')

            plottingmean(ax13,branda2000,'jugurfesta','Skalaðar einkunnir','Júgurfesta')
            plottingmean(ax13,branda2000,'jugurband','Skalaðar einkunnir','Júgurband')
            plottingmean(ax13,branda2000,'jugurdypt','Skalaðar einkunnir','Júgurdýpt')
            plottingmean(ax14,branda2000,'spenalengd','','Spenalengd')
            plottingmean(ax14,branda2000,'spenathykkt','','Spenaþykkt')
            plottingmean(ax14,branda2000,'spenastada','','Spenastaða')
            plottingmean(ax15,branda2000,'mjaltir','','Mjaltir')
            plottingmean(ax15,branda2000,'skap','','Skap')
            plottingmean(ax15,branda2000,'L3','','Ending, 1. burður til loka 3. mjaltaskeiðs')
            plottingmean(ax16,branda2000,'total','','Heildareinkunn')

            fig.suptitle(f'Meðalkynbótamat árganga {yearmonth}', fontsize=26, fontweight ="bold")
            plt.subplots_adjust(left=0.05, bottom=0.07, right=0.97, top=0.94, wspace=0.05, hspace=0.09)

            fig.set_size_inches([18, 9])
            plt.savefig(f'../{yearmonth}/figures/meanbybirthyear{yearmonth}.png')

            #Correlation heat map for traits in dataframe!
            sns.set(font_scale=0.6)
            plt.figure(figsize=(8,8))
            sns.heatmap(branda2000[['myt','fyt','pyt','fpt','ppt','milkper', 'fatper', 'protper',
            'yieldtotal','scst','CR0','ICF','IFL', 'fertility',
            'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
            'stada_hh', 'stada_ha', 'klaufhalli', 'jugurfesta',
            'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
            'mjaltir', 'skap','mjaltarod', 'gaedarod','L1','L2','L3', 'total'
            ]].corr(), annot=True, cmap='coolwarm')
            plt.title('Fylgni eiginleika (skalað kynbótamat)', fontsize = 20)

            plt.subplots_adjust(left=0.06, bottom=0.09, right=1, top=0.95, wspace=0.20, hspace=0.20)

            fig.set_size_inches([18, 9])
            plt.savefig(f'../{yearmonth}/figures/correlationheatmap{yearmonth}.png')

            if pltshow == 1:
                plt.show()

#
