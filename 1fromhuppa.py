#This program prepares DMU runs for all traits
#---------------------------------------------------------------------------
#Thordis Thorarinsdottir 2021

import pandas as pd
import numpy as np
import datetime
import os
import shutil
import subprocess
from matplotlib import pyplot as plt

#user defined functions used by program from kynbotamat_module.py
from kynbotamat_module import readfilecsv
from kynbotamat_module import readingfilefwf
from kynbotamat_module import hy_grouping
from kynbotamat_module import plottingmean

#Options in program!

#Creating DMU trait files
prepfordmu = 0

preptdm = 0 #run preptdm.f program

#Create DMU datafiles for these traits
fertilitydmu = 1
confdmu = 1
rankorderdmu = 1
longdmu = 1

#Plotting tdm1 datafile to check if OK
plottdm = 1

#Compile and run accuracy fortran programs for yield
accuracytdm = 1
#---------------------------------------------------------------------------
#File with information for program!
control = pd.read_csv(
    'control.txt',
    header=None,
    names=['control']
    )
#---------------------------------------------------------------------------
yearmonth = control.loc[0,'control']  #year and month of breeding value estamation
#---------------------------------------------------------------------------
# PART 1 - Prepping for DMU runs
#---------------------------------------------------------------------------

#Creating a new Year-month folder for breeding evaluation. From now on all
#files will be gathered there
def createpath(x):
    path = x #creation of a dir for year.month of BV estimation
    isExist = os.path.exists(path) #checks first if dir already exists
    if not isExist:
        os.makedirs(path);
        print(f'{path} is created!')
    else:
        print(f'{path} exists!')

createpath(f'../{yearmonth}')
createpath(f'../{yearmonth}/dmu_data')
createpath(f'../{yearmonth}/results')
createpath(f'../{yearmonth}/figures')

if prepfordmu == 1:
    #---------------------------------------------------------------------------
    #This function creates a directory for current DMU run
    #---------------------------------------------------------------------------
    def prep(trait,year):
        path = f'../{year}/DMU/{trait}' #creation of a dir for year.month of BV estimation
        isExist = os.path.exists(path) #checks first if dir already exists
        if not isExist:
            os.makedirs(path);
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')

        dirpath = f'../{year}/DMU/{trait}/dir' #copies dir files from dir folder to trait folder
        isExist = os.path.exists(dirpath)
        if not isExist:
            shutil.copy(f'dir/{trait}.dir', dirpath)
            print(f'{dirpath} is created!')
        else:
            print(f'{dirpath} exists!')

        parpath = f'../{year}/DMU/{trait}/{trait}.par' #copies dir files from dir folder to trait folder
        isExist = os.path.exists(parpath)
        if not isExist:
            shutil.copy(f'dir/{trait}.par', parpath)
            print(f'{parpath} is created!')
        else:
            print(f'{parpath} exists!')

        path = f'../{year}/DMU/{trait}/dmu1.sh' #copies dmu1.sh trait folder
        isExist = os.path.exists(path)
        if not isExist:
            shutil.copy(f'dmu1.sh', f'../{year}/DMU/{trait}')
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')

        path = f'../{year}/DMU/{trait}/dmu5.sh' #copies dmu5.sh trait folder
        isExist = os.path.exists(path)
        if not isExist:
            shutil.copy(f'dmu5.sh', f'../{year}/DMU/{trait}')
            print(f'{path} is created!')
        else:
            print(f'{path} exists!')

    #yearmonth variable is read from info file
    prep('my',yearmonth)
    prep('fy',yearmonth)
    prep('py',yearmonth)
    prep('fp',yearmonth)
    prep('pp',yearmonth)
    prep('scs',yearmonth)
    prep('fer',yearmonth)
    prep('conf',yearmonth)
    prep('rank',yearmonth)
    prep('long',yearmonth)


#---------------------------------------------------------------------------
# PART 1 - SORT PEDIGREE, MERGE TDM FILES AND RUN PREP_TDM IN SHELL!
#---------------------------------------------------------------------------
if preptdm == 1 :
    print('Part 1 of program')
    #Name of pedigree file
    pedigreefile = control.loc[12,'control']
    sortedped = control.loc[11,'control']
    oldtdm = control.loc[16,'control']
    newtdm = control.loc[17,'control']
    tdmfile = f'../{yearmonth}/dmu_data/tdm.csv'

    sortingped = f'sort +0.0 -0.15 {pedigreefile} -o {sortedped}'
    combinetdm = f"cat {oldtdm} {newtdm} | awk -F',' '!a[$1$5$6$7$8$9$10$11]++' > {tdmfile}"
    subprocess.call(sortingped, shell=True)
    print('Sorting pedigree in shell done')

    subprocess.call(combinetdm, shell=True)
    shutil.copy(tdmfile, f'../tdmfile/tdm{yearmonth}.csv')
    print('Combining TDM files in shell done')

    subprocess.check_output('gfortran -o preptdm preptdm_6.f', shell=True)

    print('Compiling preptdm.f done')

    subprocess.call('./preptdm', shell=True)

else:
    print( f'No prepping for TDM' )


#---------------------------------------------------------------------------
#PART 2 - INFORMATION FOR PROGRAM
#---------------------------------------------------------------------------

#Radnrkodifile
radnrkodifile = f'../{yearmonth}/dmu_data/radnrkodi' #Created by prep_tdm.f
#Format of radnrkodi file created by prep_tdm.f
radnrkodi_columns = ['id','code_id','stada','norec','fix1','fix2', 'fix3','sex']
widths_radnrkodi = [15,9,3,3,6,6,6,2]

#Date of the fertility-data collection from Huppa!
collectiondate = pd.to_datetime(control.loc[20,'control'], format='%Y%m%d')
#---------------------------------------------------------------------------
#Insemination file from Huppa for fertility
#Datafile 1, info about inseminations. One line for every ins.
insfile = control.loc[18,'control']
insfile_columns = ['id','ins','tech','comment','lact']
#Cow file from Huppa for fertility
cowfile = control.loc[19,'control']
cowfile_columns = ['id','herd','birth','death','calv1','calv2','calv3','calv4']

#Conformation file from Huppa
conformationfile = control.loc[21,'control']
conformationfile_columns=['id','domsdagur','birth','calv1','calv2','unused',
'staerd','boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurjafnvaegi', 'jugurfesta',
'jugurband', 'jugurdypt', 'spenagerd', 'spenalengd', 'spenathykkt', 'spenastada', 'spenaoddur',
'mjaltir', 'skap', 'aukaspenar', 'haed','herd']
widths_conformationfile = [15,8,8,8,8,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,7]

#Rank order file from Huppa
rankorderfile = control.loc[22,'control']
rankorderfile_columns=['instring']
widths_rankorderfile = [38]

#---------------------------------------------------------------------------
#Observationfiles used for DMU created by this program
fertilityobs = f'../{yearmonth}/dmu_data/dmu_fertility.txt'
fertilityobs_columns = ['code_id','HBY','HC1','HC2','HC3','IYM0','IYM1','IYM2',
    'IYM3','AGEi_h','AGEc_1','AGEc_2','AGEc_3','tech_h',
    'CRh','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']

confobs = f'../{yearmonth}/dmu_data/dmu_conformation.txt'
confobs_columns = ['code_id','HdomsY','lact','AGEc_1',
    'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenathykkt', 'spenastada',
    'mjaltir', 'skap']

rankorderobs = f'../{yearmonth}/dmu_data/dmu_rankorder.txt'
rankorderobs_columns = ['code_id', 'year', 'mjaltarod', 'gaedarod']

longobs = f'../{yearmonth}/dmu_data/dmu_long.txt'
longobs_columns = ['code_id','AGEc_1','CYM1','h5y','herdCY1','L1','L2','L3']
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
#Function to check if there are values in certain columns, used in fertility
def check(df,trait1,trait2,trait3,trait4,heifer,first,second,third,newcolumn, outcome):
        df.loc[
        (df[trait1].notnull().astype(int) == heifer) &
        (df[trait2].notnull().astype(int) == first) &
        (df[trait3].notnull().astype(int) == second) &
        (df[trait4].notnull().astype(int) == third)
        ,newcolumn] = outcome
        return df[newcolumn]
#---------------------------------------------------------------------------

#---------------------------------------------------------------------------
# PART 3 - CREATION OF FERTILITY, CONFORMATION AND RANK ORDER DMU FILES
#---------------------------------------------------------------------------

#-------------------------------------------------------------
#Reading in ranrkodi to replace IDs with code ids, radnrkodi is created by preptdm.f
if (fertilitydmu == 1) | (confdmu == 1) | (rankorderdmu == 1) | (longdmu == 1) | (plottdm == 1) :
    radnrkodi = readingfilefwf(radnrkodifile,radnrkodi_columns,widths_radnrkodi)
    radnrkodi = radnrkodi.drop(
        ['stada','norec','fix1','fix2', 'fix3','sex'], axis = 1)

#---------------------------------------------------------------------------
#Start of fertility program
#---------------------------------------------------------------------------
if fertilitydmu == 1:
    #Reading files
    ins_df = readfilecsv(insfile,insfile_columns, '\t')
    cows_df = readfilecsv(cowfile,cowfile_columns, '\t')

    #'ins','birth','death','calv1','calv2','calv3','calv4' formatted into dates
    ins_df['ins'] = pd.to_datetime(ins_df['ins'], format='%Y%m%d')

    cows_df[['birth','death','calv1','calv2','calv3','calv4']] = cows_df[
        ['birth','death','calv1','calv2','calv3','calv4']
        ].apply(lambda x: pd.to_datetime(x, format='%Y%m%d'))

    #---------------------------------------------------------------------------
    #Number of calvings in data counted
    #This is so wrong observations can be cleaned away
    #---------------------------------------------------------------------------
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',0,0,0,0,'no_calv',0)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,0,0,0,'no_calv',1)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,0,0,'no_calv',2)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,1,0,'no_calv',3)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,1,1,'no_calv',4)

    #Faulty obs. collected and marked as 9
    cows_df.loc[
    (cows_df['no_calv'] != 0) & (cows_df['no_calv'] != 1) & (cows_df['no_calv'] != 2)
    & (cows_df['no_calv'] != 3) & (cows_df['no_calv'] != 4),
    'no_calv'] = 9
    #---------------------------------------------------------------------------

    #Datafiles sorted before merge
    cows_df = cows_df.sort_values(by=['id'])
    ins_df = ins_df.sort_values(by=['id','ins'])

    #---------------------------------------------------------------------------
    #Two files merged into one file where info about cows follows each ins.
    df = pd.merge(left=ins_df, right=cows_df,on='id' )
    print( f'Gripalisti and saedingar have been merged' )
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    # In which lactation does the ins occur?
    # --------------------------------------
    # Heifers who are inseminated but do not calve
    df.loc[
    (df['no_calv'] == 0), 'ins_lact'] = 7
    #Heifers that calve inseminations
    df.loc[
    (df['no_calv'] > 0) & (df['ins'] < df['calv1']), 'ins_lact'] = 0
    #Lactations 1
    df.loc[
    (df['no_calv'] >= 1) & (df['ins'] > df['calv1']) & (df['ins'] < df['calv2']),
    'ins_lact'] = 1
    df.loc[
    (df['no_calv'] == 1) & (df['ins'] > df['calv1']), 'ins_lact'] = 1
    #Lactations 2
    df.loc[
    (df['no_calv'] >= 2) & (df['ins'] > df['calv2']) & (df['ins'] < df['calv3']),
    'ins_lact'] = 2
    df.loc[
    (df['no_calv'] == 2) & (df['ins'] > df['calv2']), 'ins_lact'] = 2
    #Lactations 3
    df.loc[
    (df['no_calv'] >= 3) & (df['ins'] > df['calv3']) & (df['ins'] < df['calv4']),
    'ins_lact'] = 3
    df.loc[
    (df['no_calv'] == 3) & (df['ins'] > df['calv3']), 'ins_lact'] = 3
    # --------------------------------------
    #Wrong ins located and marked as 99
    df.loc[
    (df['ins_lact'] != 0.0) & (df['ins_lact'] != 7.0) & (df['ins_lact'] != 1.0)
    & (df['ins_lact'] != 2.0) & (df['ins_lact'] != 3.0),
    'ins_lact'] = 99.0
    #Ins that happen before 2008 are marked as 88
    df.loc[
    (df['ins'] < '2008-01-01'),
    'ins_lact'] = 88.0
    #Ins have comments located are marked as 77
    df.loc[
    (df['comment'].notnull().astype(int) == 1),
    'ins_lact'] = 77.0
    # --------------------------------------
    print( f'Lactation for inseminations have been found' )
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    #Splitting info into lactations and counting ins per lact.
    #---------------------------------------------------------------------------
    heifers_df = df.loc[(df['ins_lact'] == 0.0) | (df['ins_lact'] == 7.0)]
    lact1_df = df.loc[df['ins_lact'] == 1.0]
    lact2_df = df.loc[df['ins_lact'] == 2.0]
    lact3_df = df.loc[df['ins_lact'] == 3.0]
    #---------------------------------------------------------------------------

    #Wrong ins were marked and collected
    wrongins = df.loc[(df['ins_lact'] == 99.0) | (df['ins_lact'] == 88.0) |
    (df['ins_lact'] == 77.0)]
    wrongins.loc[:,'check'] = 1
    wrongins = wrongins[['id','check']]
    wrongins.columns = ['id', 'ins_lact']
    wrongins = wrongins.drop_duplicates(subset=['id'])
    #They are merged into the cows dataframe so cows can be cleaned away later
    cows_df = pd.merge(left=cows_df[
    ['id','herd','birth','death','calv1','calv2','calv3','calv4','no_calv']
    ], right=wrongins[['id','ins_lact']], on='id', how='left')
    #---------------------------------------------------------------------------

    #---------------------------------------------------------------------------
    #This is a function to count ins per lactation and collect the first and last ins
    #Results are merged into main dataframe
    def inscount (df1,firstins, lastins, T, tech, main):
        df1.loc[:,'T'] = df1.groupby('id')['id'].transform('count') #Total ins
        df1.loc[:,'No'] = df1.groupby('id')['id'].cumcount() + 1 #number of each ins
        #The last ins
        df1.loc[(df1['T'] == df1['No']), lastins] = df1['ins']
        #Tota number of heifer ins
        df1.loc[(df1['T'] == df1['No']), T] = df1['T']
        #First ins
        df1.loc[(df1['No'] == 1), firstins] = df1['ins']
        #First heifer tech
        df1.loc[(df1['No'] == 1), tech] = df1['tech']
        #First and last observations locaed
        first = df1[df1[firstins] == df1['ins'] ]
        last = df1[df1[lastins] == df1['ins']]
        #First and last observations put in seperate dataframes
        first = first.loc[:, ['id',firstins,tech]]
        last = last.loc[:, ['id',lastins,T]]
        #First and last heifer ins merged into one dataframe
        df1 = pd.merge(left=first, right=last)
        #Removed obs that happen within 302 days before collection date!
        df1['check'] = (collectiondate - df1[firstins]).dt.days
        df1 = df1.loc[(df1['check'] > 302 )]
        df1 = df1.drop(['check'], axis = 1)
        main = main.merge(df1, on='id', how='outer')
        return main
    #---------------------------------------------------------------------------
    cows_df = inscount(heifers_df, 'first_h', 'last_h', 'T_h', 'tech_h',cows_df )
    cows_df = inscount(lact1_df, 'first_1', 'last_1', 'T_1', 'tech_1',cows_df )
    cows_df = inscount(lact2_df, 'first_2', 'last_2', 'T_2', 'tech_2',cows_df )
    cows_df = inscount(lact3_df, 'first_3', 'last_3', 'T_3', 'tech_3',cows_df )
    #---------------------------------------------------------------------------
    print( f'First and last inseminations have been found' )

    #---------------------------------------------------------------------------
    #Feritlity traits and traits for cleaning and fixed effecsts
    #---------------------------------------------------------------------------
    #Age at first ins as unborn heifers
    cows_df.loc[:,'AGEi_h'] = (cows_df['first_h'] - cows_df['birth']).dt.days
    #Age at calving at each lact
    cows_df.loc[:,'AGEc_1'] = (cows_df['calv1'] - cows_df['birth']).dt.days
    cows_df.loc[:,'AGEc_2'] = (cows_df['calv2'] - cows_df['birth']).dt.days
    cows_df.loc[:,'AGEc_3'] = (cows_df['calv3'] - cows_df['birth']).dt.days
    #Gestation lenght
    cows_df.loc[:,'gest1'] = (cows_df['calv1'] - cows_df['last_h']).dt.days
    cows_df.loc[:,'gest2'] = (cows_df['calv2'] - cows_df['last_1']).dt.days
    cows_df.loc[:,'gest3'] = (cows_df['calv3'] - cows_df['last_2']).dt.days
    #Calving interval
    cows_df.loc[:,'CI12'] = (cows_df['calv2'] - cows_df['calv1']).dt.days
    cows_df.loc[:,'CI23'] = (cows_df['calv3'] - cows_df['calv2']).dt.days
    cows_df.loc[:,'CI34'] = (cows_df['calv4'] - cows_df['calv3']).dt.days
    #---------------------------------------------------------------------------
    #Fertility trait, days between calving and first ins
    cows_df.loc[:,'ICF1'] = (cows_df['first_1'] - cows_df['calv1']).dt.days
    cows_df.loc[:,'ICF2'] = (cows_df['first_2'] - cows_df['calv2']).dt.days
    cows_df.loc[:,'ICF3'] = (cows_df['first_3'] - cows_df['calv3']).dt.days
    #Fertility trait, days between first and last ins
    cows_df.loc[:,'IFLh'] = (cows_df['last_h'] - cows_df['first_h']).dt.days
    cows_df.loc[:,'IFL1'] = (cows_df['last_1'] - cows_df['first_1']).dt.days
    cows_df.loc[:,'IFL2'] = (cows_df['last_2'] - cows_df['first_2']).dt.days
    cows_df.loc[:,'IFL3'] = (cows_df['last_3'] - cows_df['first_3']).dt.days
    #1-4 days count as one ins period so IFL set to zero
    cows_df.loc[(cows_df['IFLh'] <= 4) &(cows_df['IFLh'] >= 1),'IFLh'] = 0
    cows_df.loc[(cows_df['IFL1'] <= 4) &(cows_df['IFL1'] >= 1),'IFL1'] = 0
    cows_df.loc[(cows_df['IFL2'] <= 4) &(cows_df['IFL2'] >= 1),'IFL2'] = 0
    cows_df.loc[(cows_df['IFL3'] <= 4) &(cows_df['IFL3'] >= 1),'IFL3'] = 0
    #---------------------------------------------------------------------------
    #Obs where IFL is registered as 0 but no calving has occured after the ins
    #are then set to IFL = 21 days
    cows_df.loc[(cows_df['IFL1'] == 0) &(cows_df['calv2'].notnull().astype(int) == 0),'IFL1'] = 21
    cows_df.loc[(cows_df['IFL2'] == 0) &(cows_df['calv3'].notnull().astype(int) == 0),'IFL2'] = 21
    cows_df.loc[(cows_df['IFL3'] == 0) &(cows_df['calv4'].notnull().astype(int) == 0),'IFL3'] = 21
    #---------------------------------------------------------------------------
    #Fertility trait, conception rate at first calving for heifers
    #1 means the cow seems to be pregnant after first ins
    #0 means the cow did not get pregnant after first ins
    cows_df.loc[
    (cows_df['IFLh'] == 0) &
    (cows_df['calv1'].notnull().astype(int) == 1), #calv1, there is an obsv.
    'CRh'] = 1
    cows_df.loc[
    (cows_df['IFLh'] >= 5) &
    (cows_df['calv1'].notnull().astype(int) == 1),
    'CRh'] = 0
    cows_df.loc[
    (cows_df['IFLh'] >= 0) &
    (cows_df['calv1'].notnull().astype(int) == 0), #calv1, there is no obsv.
    'CRh'] = 0

    print( f'Fertility traits have been created' )
    #---------------------------------------------------------------------------
    #---------------------------------------------------------------------------
    #Cleaning starts, animals who are sorted out can be collected into a seperate file
    # a) Age at first calving 550-1100 days
    # b) Age at first ins 270-900 days c) Calving interval 208-600 days
    # d) Gestastion lenght 260-302 days e) Number of ins 1-8 per lact
    # f) ICF 20-230 days g) IFL 0-365 days
    #Other cleaning:
    #Herd missing, Insemanation years only after 2008-01-01
    #Records for later lactations were excluded if information about previous lactations
        #were not available. Ins with comments cleaned away
    #--------------------------------------------------------------------------
    # a) Age at first calving 550-1100 days
    cows_df.loc[
    (cows_df['AGEc_1'] < 550) | (cows_df['AGEc_1'] > 1100),
    'wrong'] = 'a'
    # b) Age at first ins 270-900 days
    cows_df.loc[
    (cows_df['AGEi_h'] < 270) | (cows_df['AGEi_h'] > 900),
    'wrong'] = 'b'
    # c) Calving interval 208-600 days
    cows_df.loc[
    (cows_df['CI12'] < 208) | (cows_df['CI12'] > 600),#lact1
    'first_1'] = np.nan
    cows_df.loc[
    (cows_df['CI23'] < 208) | (cows_df['CI23'] > 600), #lact2
    'first_2'] = np.nan
    cows_df.loc[
    (cows_df['CI34'] < 208) | (cows_df['CI34'] > 600), #lact3
    'first_3'] = np.nan
    # c) Calving interval 208-600 days
    # d) Gestastion lenght 260-302 days
    cows_df.loc[
    (cows_df['gest1'] < 260) | (cows_df['gest1'] > 302),
    'first_1'] = np.nan
    cows_df.loc[
    (cows_df['gest2'] < 260) | (cows_df['gest2'] > 302),
    'first_2'] = np.nan
    cows_df.loc[
    (cows_df['gest3'] < 260) | (cows_df['gest3'] > 302),
    'first_3'] = np.nan
    # e) Number of ins 1-8 per lact
    cows_df.loc[(cows_df['T_h'] > 8) , 'first_h'] = np.nan
    cows_df.loc[(cows_df['T_1'] > 8) , 'first_1'] = np.nan
    cows_df.loc[(cows_df['T_2'] > 8) , 'first_2'] = np.nan
    cows_df.loc[(cows_df['T_3'] > 8) , 'first_3'] = np.nan
    # f) ICF 20-230 days
    cows_df.loc[(cows_df['ICF1'] < 20) | (cows_df['ICF1'] > 230), 'first_1'] = np.nan
    cows_df.loc[(cows_df['ICF2'] < 20) | (cows_df['ICF2'] > 230), 'first_2'] = np.nan
    cows_df.loc[(cows_df['ICF3'] < 20) | (cows_df['ICF3'] > 230), 'first_3'] = np.nan
    # g) IFL 0-365 days
    cows_df.loc[(cows_df['IFL1'] > 365), 'first_1'] = np.nan
    cows_df.loc[(cows_df['IFL2'] > 365), 'first_2'] = np.nan
    cows_df.loc[(cows_df['IFL3'] > 365), 'first_3'] = np.nan

    #---------------------------------------------------------------------------
    #Checking if cows have the correct order of ins registered
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',1,0,0,0,'check',1000)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',1,1,0,0,'check',1100)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',1,1,1,0,'check',1110)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',1,1,1,1,'check',1111)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',0,1,0,0,'check',100)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',0,1,1,0,'check',110)
    cows_df['check'] = check(cows_df,'first_h','first_1','first_2','first_3',0,1,1,1,'check',111)


    #Above all wrong obsv. were marked, now observ. are sorted into a correct
    #datafile and a wrong datafile
    #---------------------------------------------------------------------------
    #Observations that fulfill every condition are collected
    data_use = cows_df[(
        cows_df['id'] < 900000000000000) & (  #
        cows_df['herd'].notnull().astype(int) == 1) & (  #herd missing
        cows_df['birth'].notnull().astype(int) == 1) & (  #no birth year
        cows_df['wrong'].notnull().astype(int) == 0) & ( #AFC or AFI wrong
        cows_df['ins_lact'] != 1.0) & ( #ins did not register on a lactation #ins happens before 2008 #ins has a comment
        cows_df['no_calv'] != 9.0) & ( #Order of calving not correct
        cows_df['check'].notnull().astype(int) == 1)] #correct orders of ins.
    print( f'First cleaning is finished' )

    #Observations that DON'T fulfill every condition are collected can be collected
    # data_not_used = cows_df[(
    #     cows_df['herd'].notnull().astype(int) == 0) | (  #herd missing
    #     cows_df['birth'].notnull().astype(int) == 0) | ( #no birth year
    #     cows_df['wrong'].notnull().astype(int) == 1) | ( #AFC or AFI wrong
    #     cows_df['ins_lact'] == 1.0) | (  #ins did not register on a lactation #ins happens before 2008 #ins has a comment
    #     cows_df['no_calv'] == 9.0) & (  #Order of calving not correct
    #     cows_df['check'].notnull().astype(int) == 0)]
    #data_not_used.to_excel("../data/data_not_used.xlsx", index=False)

    #Dropping unused columns
    data_use = data_use.drop(
        ['wrong', 'check', 'ins_lact', 'no_calv',
        'last_h','last_1','last_2','last_3'], axis = 1)
    #---------------------------------------------------------------------------
    #Filling in for missing values for technician in heifer ins
    data_use['tech_h'] = data_use['tech_h'].fillna(0, downcast='infer')
    #Count no of obs for each tech
    data_use['tech_h_c'] = data_use.groupby('tech_h')['tech_h'].transform('count')
    #Delete if less than 3
    data_use = data_use[(data_use['tech_h_c'] > 2) ]
    #---------------------------------------------------------------------------
    #Fixed effects - Age at first ins in MONTHS - Age at calving in MONTHS
    data_use[['AGEi_h','AGEc_1','AGEc_2','AGEc_3']] = data_use[
        ['AGEi_h','AGEc_1','AGEc_2','AGEc_3']
        ].apply(
        lambda s: (s / 30.5).apply(np.ceil).fillna(0).astype(int))
    #Counting how many cows are in each age group
    data_use['age0'] = data_use.groupby('AGEi_h')['AGEi_h'].transform('count')
    data_use['age1'] = data_use.groupby('AGEc_1')['AGEc_1'].transform('count')
    data_use['age2'] = data_use.groupby('AGEc_2')['AGEc_2'].transform('count')
    data_use['age3'] = data_use.groupby('AGEc_3')['AGEc_3'].transform('count')
    #Delete rows that are less then 3 in a group
    data_use = data_use[
        (data_use['age0'] > 2) &
        (data_use['age1'] > 2) &
        (data_use['age2'] > 2) &
        (data_use['age3'] > 2)]

    #---------------------------------------------------------------------------
    #Creating year columns for herd-year fixed effect
    data_use[['BY','C1','C2','C3']] = data_use[
        ['birth','calv1','calv2','calv3']
        ].apply(
        lambda s: s.dt.strftime('%Y'))

    #Create a dataframe to group herd-year groups
    ids = data_use['id'].copy()
    herd_by = data_use[['id','herd','birth','BY','first_h']].copy()
    herd_c1 = data_use[['id','herd','calv1','C1','first_1']].copy()
    herd_c2 = data_use[['id','herd','calv2','C2','first_2']].copy()
    herd_c3 = data_use[['id','herd','calv3','C3','first_3']].copy()

    #Only include relevant rows
    herd_by = herd_by.loc[(herd_by['herd'].notnull().astype(int) == 1) &
        (herd_by['birth'].notnull().astype(int) == 1) &
        (herd_by['first_h'].notnull().astype(int) == 1)]

    herd_c1 = herd_c1.loc[(herd_c1['herd'].notnull().astype(int) == 1) &
        (herd_c1['calv1'].notnull().astype(int) == 1) &
        (herd_c1['first_1'].notnull().astype(int) == 1)]

    herd_c2 = herd_c2.loc[(herd_c2['herd'].notnull().astype(int) == 1) &
        (herd_c2['calv2'].notnull().astype(int) == 1) &
        (herd_c2['first_2'].notnull().astype(int) == 1)]

    herd_c3 = herd_c3.loc[(herd_c3['herd'].notnull().astype(int) == 1) &
        (herd_c3['calv3'].notnull().astype(int) == 1) &
        (herd_c3['first_3'].notnull().astype(int) == 1)]

    #Sort the dataframes by relevant date
    herd_by = herd_by.sort_values(by=['herd','birth'])
    herd_c1 = herd_c1.sort_values(by=['herd','calv1'])
    herd_c2 = herd_c2.sort_values(by=['herd','calv2'])
    herd_c3 = herd_c3.sort_values(by=['herd','calv3'])

    #---------------------------------------------------------------------------
    #Calling function above to create herd-year groups in all lactations
    #---------------------------------------------------------------------------
    print( f'Starting with herd x birth year' )
    ids = hy_grouping('herd', herd_by, 'BY', 'HBY', ids)

    print( f'Starting with herd x calving year 1' )
    ids = hy_grouping('herd', herd_c1, 'C1', 'HC1', ids)

    print( f'Starting with herd x calving year 2' )
    ids = hy_grouping('herd', herd_c2, 'C2', 'HC2', ids)

    print( f'Starting with herd x calving year 3' )
    ids = hy_grouping('herd', herd_c3, 'C3', 'HC3', ids)

    #---------------------------------------------------------------------------
    #Merging herd-year groups with main dataframe
    data_use = pd.merge(left=data_use, right=ids, on='id', how='left')
    #---------------------------------------------------------------------------

    #Counting how many cows are in each herd-year group
    data_use['HBY_c'] = data_use.groupby('HBY')['HBY'].transform('count')
    data_use['HC1_c'] = data_use.groupby('HC1')['HC1'].transform('count')
    data_use['HC2_c'] = data_use.groupby('HC2')['HC2'].transform('count')
    data_use['HC3_c'] = data_use.groupby('HC3')['HC3'].transform('count')
    #Delete rows that are alone in a group
    data_use = data_use[
        (data_use['HBY_c'] > 2) &
        (data_use['HC1_c'] > 2) &
        (data_use['HC2_c'] > 2) &
        (data_use['HC3_c'] > 2)]

    #---------------------------------------------------------------------------
    #Creating year columns for ins-year-month fixed effect
    data_use[['IY0','IY1','IY2','IY3']] = data_use[
        ['first_h','first_1','first_2','first_3']
        ].apply(
        lambda s: s.dt.strftime('%Y').replace('NaT', '0').astype(int))
    data_use[['IM0','IM1','IM2','IM3']] = data_use[
        ['first_h','first_1','first_2','first_3']
        ].apply(
        lambda s: s.dt.strftime('%m').replace('NaT', '0').astype(int))
    #---------------------------------------------------------------------------
    #Create a dataframe for ins-year-month fixed effect
    ids = data_use['id'].copy()
    ym0 = data_use[['id','IY0','IM0']].copy()
    ym1 = data_use[['id','IY1','IM1']].copy()
    ym2 = data_use[['id','IY2','IM2']].copy()
    ym3 = data_use[['id','IY3','IM3']].copy()

    #Only include relevant rows
    ym0 = ym0.loc[(ym0['IY0'] > 0)]
    ym1 = ym1.loc[(ym1['IY1'] > 0)]
    ym2 = ym2.loc[(ym2['IY2'] > 0)]
    ym3 = ym3.loc[(ym3['IY3'] > 0)]
    #Sort the dataframes by relevant date
    ym0 = ym0.sort_values(by=['IY0','IM0'])
    ym1 = ym1.sort_values(by=['IY1','IM1'])
    ym2 = ym2.sort_values(by=['IY2','IM2'])
    ym3 = ym3.sort_values(by=['IY3','IM3'])

    #---------------------------------------------------------------------------
    #Calling function above to create ins-year-month fixed effect
    #---------------------------------------------------------------------------
    print( f'Starting with ins year x month 0' )
    ids = hy_grouping('IY0', ym0, 'IM0', 'IYM0', ids)
    print( f'Starting with ins year x month 1' )
    ids = hy_grouping('IY1', ym1, 'IM1', 'IYM1', ids)
    print( f'Starting with ins year x month 2' )
    ids = hy_grouping('IY2', ym2, 'IM2', 'IYM2', ids)
    print( f'Starting with ins year x month 3' )
    ids = hy_grouping('IY3', ym3, 'IM3', 'IYM3', ids)
    #---------------------------------------------------------------------------
    #Merging herd-year groups with main dataframe
    data_use = pd.merge(left=data_use, right=ids, on='id', how='left')
    #---------------------------------------------------------------------------

    #Counting how many cows are in each  ins-year-month group
    data_use['IYM0_c'] = data_use.groupby('IYM0')['IYM0'].transform('count')
    data_use['IYM1_c'] = data_use.groupby('IYM1')['IYM1'].transform('count')
    data_use['IYM2_c'] = data_use.groupby('IYM2')['IYM2'].transform('count')
    data_use['IYM3_c'] = data_use.groupby('IYM3')['IYM3'].transform('count')
    #Delete rows that are less then 3
    data_use = data_use[(data_use['IYM0_c'] > 2) & (data_use['IYM1_c'] > 2) &
        (data_use['IYM2_c'] > 2) & (data_use['IYM3_c'] > 2)]
    #---------------------------------------------------------------------------

    #Filling in missin values for DMU, -999.0 for reals
    #Real columns
    realc = ['CRh','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']
    data_use[realc] = data_use[realc].fillna(-999.0)
    #---------------------------------------------------------------------------
    data_use = pd.merge(left=data_use, right=radnrkodi[['id','code_id']], on='id', how='left')
    data_use = data_use.sort_values(by=['code_id'])
    data_use = data_use[data_use['code_id'].notnull().astype(int) == 1]
    data_use['code_id'] = data_use['code_id'].astype(int)
    #---------------------------------------------------------------------------
    #Creating the DMU datafile
    dmu_fertility = data_use[fertilityobs_columns].copy()
    # DMU datafile
    dmu_fertility.to_csv(fertilityobs, index=False, header=False, sep=' ')

    print(dmu_fertility.iloc[50000:50015])
    print(dmu_fertility.info())

else:
    print( f'DMU fertility file not created' )

#---------------------------------------------------------------------------
#Start of conformation  program
#---------------------------------------------------------------------------
if confdmu == 1:

    df = readingfilefwf(conformationfile,conformationfile_columns,widths_conformationfile)

    #Only cows allowed with id
    df = df.loc[(df['id'].notnull().astype(int) == 1)]
    #ID changed to integer
    df['id'] = df['id'].astype(int)
    df = df.loc[df['id'] < 900000000000000]

    #Dropping unused columns
    df = df.drop(
        ['unused', 'haed', 'aukaspenar','staerd',
        'spenaoddur','spenagerd','jugurjafnvaegi'], axis = 1)

    #Cleaning away values that are not allowed
    df = (df.loc[(df['id'].notnull().astype(int) == 1) &
                (df['boldypt'] != 0) &
                (df['utlogur'] != 0) &
                (df['yfirlina'] != 0) &
                (df['malabreidd'] != 0) &
                (df['malahalli'] != 0) &
                (df['malabratti'] != 0) &
                (df['stada_haekla_hlid'] != 0) &
                (df['stada_haekla_aftan'] != 0) &
                (df['klaufhalli'] != 0) &
                (df['jugurfesta'] != 0) &
                (df['jugurband'] != 0) &
                (df['jugurdypt'] != 0) &
                (df['spenalengd'] != 0) &
                (df['spenathykkt'] != 0) &
                (df['spenastada'] != 0) &
                (df['mjaltir'] != 0) &
                (df['skap'] != 0) &
                (df['domsdagur'] > 20000101)
                ])

    #Creating datetime from date columns
    df[['domsdagur','birth','calv1','calv2']] = df[
        ['domsdagur','birth','calv1','calv2']
        ].apply(
        lambda x: pd.to_datetime(x, format='%Y%m%d'))

    #Only allowed if birthdate and calving date 1 are registered
    df = df.loc[(df['birth'].notnull().astype(int) == 1) & (df['calv1'].notnull().astype(int) == 1)]

    #Creating a judging year column
    df['domsY'] = df['domsdagur'].apply(lambda s: s.strftime('%Y'))

    #Finding info about which lactation the judging occured
    #Lactations 1
    df.loc[
    (df['domsdagur'] > df['calv1']) & (df['domsdagur'] < df['calv2']),
    'lact'] = 1
    df.loc[
    (df['domsdagur'] > df['calv1']) & (df['calv2'].notnull().astype(int) == 0), 'lact'] = 1
    #Lactations 2
    df.loc[
    (df['domsdagur'] > df['calv2']), 'lact'] = 2

    #Only allowed on lact 1 or 2
    df = df.loc[(df['lact'].notnull().astype(int) == 1)]

    #Creating a age at first calving variable
    df.loc[:,'AGEc_1'] = (df['calv1'] - df['birth']).dt.days
    df = df.loc[(df['AGEc_1'] > 550)] #Only older than 550 days at first calving allowed
    df.loc[:,'AGEc_1'] = (df['AGEc_1']/ 30.5).astype(int) #turning into months
    #Deleting if fewer than 3 cows in an age group!
    df['AGEc_1c'] = df.groupby('AGEc_1')['AGEc_1'].transform('count')
    df = df.loc[(df['AGEc_1c'] > 2)]
    #Only cows with herd allowed
    df = (df.loc[(df['herd'].notnull().astype(int) == 1)])#.astype(int)

    #---------------------------------------------------------------------------
    #Sort the dataframes by relevant date
    df = df.sort_values(by=['herd','domsY'])
    #---------------------------------------------------------------------------
    #Calling function above to create herd-year groups in all lactations
    #---------------------------------------------------------------------------
    df = hy_grouping('herd', df, 'domsY', 'HdomsY', df)

    #Counting how many cows are in each herd-year group
    df['HdomsYc'] = df.groupby('HdomsY')['HdomsY'].transform('count')
    #Delete rows that are alone in a group
    df = df[(df['HdomsYc'] > 2) ]

    df = pd.merge(left=df, right=radnrkodi[['id','code_id']], on='id', how='left')
    df = df.sort_values(by=['code_id'])
    df = df[df['code_id'].notnull().astype(int) == 1]
    df[['code_id','lact']] = df[['code_id','lact']].astype(int)
    #---------------------------------------------------------------------------
    #Creating the DMU datafile
    dmu_conf = df[confobs_columns].copy()

    # DMU datafile
    dmu_conf.to_csv(confobs, index=False, header=False, sep=' ')

    print(dmu_conf.iloc[80000:80015])
    print(dmu_conf.info())

else:
    print( f'DMU conformation file not created' )

#---------------------------------------------------------------------------
#Start of rank order program
#---------------------------------------------------------------------------
if rankorderdmu == 1:

    #Loading in main datafile with milking traits
    #Maybe a new datafile needs to be made!
    df = readingfilefwf(rankorderfile,rankorderfile_columns,widths_rankorderfile)

    #Format of file is strange, no other way than to read in as one whole string
    #and extract the needed variables
    df['id'] = (df.instring.astype(str).str[:15]).astype(int)
    df['year'] = (df.instring.astype(str).str[15:19]).astype(int)
    df['mjaltarod'] = (df.instring.astype(str).str[19:20]).astype(int)
    df['gaedarod'] = (df.instring.astype(str).str[26:27]).astype(int)

    df['BY'] = (df.id.astype(str).str[:4]).astype(int)

    df = df.loc[(df['BY'] < 2100) & (df['BY'] > 1990) ]
    df = df.loc[(df['mjaltarod'] > 0) & (df['gaedarod'] > 0) ]
    df = df.drop(['instring'], axis = 1)

    df = df.loc[df['id'] < 900000000000000]

    df = pd.merge(left=df, right=radnrkodi[['id','code_id']], on='id', how='left')
    df = df.sort_values(by=['code_id'])
    df = df[df['code_id'].notnull().astype(int) == 1]
    df[['code_id']] = df[['code_id']].astype(int)

    dmurankorder = df[rankorderobs_columns]

    # DMU datafile
    dmurankorder.to_csv(rankorderobs, index=False, header=False, sep=' ')

    print(dmurankorder.iloc[80000:80015])
    print(dmurankorder.info())

else:
    print( f'DMU rank order file not created' )


#---------------------------------------------------------------------------
#Start of longevity  program
#---------------------------------------------------------------------------
if longdmu == 1:

    #Cow file from fertility evaluation read
    cows_df = readfilecsv(cowfile,cowfile_columns, '\t')

    #Dates formatted
    cows_df[['birth','death','calv1','calv2','calv3','calv4']] = cows_df[
    ['birth','death','calv1','calv2','calv3','calv4']
    ].apply(lambda x: pd.to_datetime(x, format='%Y%m%d'))

    #Checking if calving are OK
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',0,0,0,0,'no_calv',0)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,0,0,0,'no_calv',1)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,0,0,'no_calv',2)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,1,0,'no_calv',3)
    cows_df['check'] = check(cows_df,'calv1','calv2','calv3','calv4',1,1,1,1,'no_calv',4)

    #Faulty obs. collected and marked as 9
    cows_df.loc[
    (cows_df['no_calv'] != 0) & (cows_df['no_calv'] != 1) & (cows_df['no_calv'] != 2)
    & (cows_df['no_calv'] != 3) & (cows_df['no_calv'] != 4),
    'no_calv'] = 9

    #Cows only allowed if no of calving and check has registered
    cows_df = cows_df.loc[(cows_df['no_calv'].notnull().astype(int) == 1) & (cows_df['check'].notnull().astype(int) == 1)]

    #Cows only allowed if they have dates for birth and first caving
    cows_df = cows_df.loc[(cows_df['birth'].notnull().astype(int) == 1)]
    cows_df = cows_df.loc[(cows_df['calv1'].notnull().astype(int) == 1)]

    #Find age at first calving
    cows_df.loc[:,'AGEc_1'] = (cows_df['calv1'] - cows_df['birth']).dt.days
    #Only allowed within a certain range (from NAV)
    cows_df.loc[(cows_df['AGEc_1'] < 450) | (cows_df['AGEc_1'] > 1280),'wrong'] = 'a'
    cows_df = cows_df.loc[(cows_df['wrong'].notnull().astype(int) == 0)]

    #Finding in which lactation a culling has occured
    cows_df.loc[
    (cows_df['no_calv'] == 0) & (cows_df['death'].notnull().astype(int) == 1), 'culled'] = 0
    cows_df.loc[
    (cows_df['no_calv'] == 1) & (cows_df['death'].notnull().astype(int) == 1), 'culled'] = 1
    cows_df.loc[
    (cows_df['no_calv'] == 2) & (cows_df['death'].notnull().astype(int) == 1), 'culled'] = 2
    cows_df.loc[
    (cows_df['no_calv'] == 3) & (cows_df['death'].notnull().astype(int) == 1), 'culled'] = 3

    #Creation of L1, L2 and L3
    cows_df.loc[(cows_df['calv2'].notnull().astype(int) == 1),'L1'] = (cows_df['calv2'] - cows_df['calv1']).dt.days
    cows_df.loc[(cows_df['L1'].notnull().astype(int) == 1),'L1'] = 365 #If cow calved again, set to max value
    #If cow is culled then days from calv 1 to culling
    cows_df.loc[(cows_df['culled'] == 1),'L1'] = (cows_df['death'] - cows_df['calv1']).dt.days
    #If more then max days, set to max days
    cows_df.loc[(cows_df['L1'] > 365), 'L1'] = 365

    cows_df.loc[(cows_df['calv3'].notnull().astype(int) == 1),'L2'] = (cows_df['calv3'] - cows_df['calv1']).dt.days
    cows_df.loc[(cows_df['L2'].notnull().astype(int) == 1),'L2'] = 730
    cows_df.loc[(cows_df['culled'] == 2),'L2'] = (((cows_df['death'] - cows_df['calv2']).dt.days) + cows_df['L1'])
    cows_df.loc[(cows_df['L2'] > 730), 'L2'] = 730

    cows_df.loc[(cows_df['calv4'].notnull().astype(int) == 1),'L3'] = (cows_df['calv4'] - cows_df['calv1']).dt.days
    cows_df.loc[(cows_df['L3'].notnull().astype(int) == 1),'L3'] = 1095
    cows_df.loc[(cows_df['culled'] == 3),'L3'] = (((cows_df['death'] - cows_df['calv3']).dt.days) + cows_df['L2'])
    cows_df.loc[(cows_df['L3'] > 1095), 'L3'] = 1095

    #Check the collectiondate if cow has been able to collect values for traits
    cows_df.loc[
        (cows_df['calv1'].notnull().astype(int) == 1) &
        (cows_df['death'].notnull().astype(int) == 0) &
        ((collectiondate - cows_df['calv1']).dt.days > 365)
        ,'L1'] = 365

    cows_df.loc[
        (cows_df['calv2'].notnull().astype(int) == 1) &
        (cows_df['death'].notnull().astype(int) == 0) &
        ((collectiondate - cows_df['calv2']).dt.days > 365)
        ,'L2'] =730

    cows_df.loc[
        (cows_df['calv3'].notnull().astype(int) == 1) &
        (cows_df['death'].notnull().astype(int) == 0) &
        ((collectiondate - cows_df['calv3']).dt.days > 365)
        ,'L3'] =1095

    #Only cows allowed if they have L1 obs
    cows_df = cows_df.loc[(cows_df['L1'].notnull().astype(int) == 1)]

    #An extra check for cows and their must have observations
    cows_df = cows_df.loc[(cows_df['L1'] > 0)]
    cows_df = cows_df.loc[(cows_df['L2'] > 0) | (cows_df['L2'].notnull().astype(int) == 0) ]
    cows_df = cows_df.loc[(cows_df['L3'] > 0) | (cows_df['L3'].notnull().astype(int) == 0)]

    #Set age to months
    cows_df['AGEc_1'] = (cows_df['AGEc_1']/30.5).astype(int)

    #Create fixed effect and count in groups
    cows_df['age1'] = cows_df.groupby('AGEc_1')['AGEc_1'].transform('count')
    cows_df = cows_df[(cows_df['age1'] > 2)]

    #Collect years from dates to create fixed effects
    cows_df.loc[:, 'BY'] = (cows_df['birth'].dt.strftime('%Y')).astype(int)
    cows_df.loc[:, 'CY1'] = (cows_df['calv1'].dt.strftime('%Y')).astype(int)
    cows_df.loc[:, 'CM1'] = (cows_df['calv1'].dt.strftime('%m')).astype(int)

    #Creation of fixed effect herd * 5 year period
    def herd5year(df,year):
        df.loc[((df[year] >= 1994) &  (df[year] <= 1998)),'ygroup'] = 1
        df.loc[((df[year] >= 1999) &  (df[year] <= 2003)),'ygroup'] = 2
        df.loc[((df[year] >= 2004) &  (df[year] <= 2008)),'ygroup'] = 3
        df.loc[((df[year] >= 2009) &  (df[year] <= 2013)),'ygroup'] = 4
        df.loc[((df[year] >= 2014)&  (df[year] <= 2018)),'ygroup'] = 5
        df.loc[((df[year] >= 2019)),'ygroup'] = 6
        df.loc[:,'h5y'] = (df['herd'].astype('str') + df['ygroup'].astype('str')).astype(float).astype(int)
        return df['h5y']

    cows_df['h5y'] = herd5year(cows_df,'BY')
    cows_df['h5y_c'] = cows_df.groupby('h5y')['h5y'].transform('count')
    cows_df = cows_df[(cows_df['h5y_c'] > 5)]

    #Create seperate df for herd year grouping
    df = cows_df[['id','herd','CY1', 'CM1']]

    df = hy_grouping('CY1', df, 'CM1', 'CYM1', df)
    df['CYM1c'] = df.groupby('CYM1')['CYM1'].transform('count')
    df = df[(df['CYM1c'] > 2)]

    df = hy_grouping('herd', df, 'CY1', 'herdCY1', df)

    #Merging into main df
    cows_df = pd.merge(left=cows_df, right=df, on='id', how='left')

    #Checking size of groups
    cows_df['HCY1_c'] = cows_df.groupby('herdCY1')['herdCY1'].transform('count')
    cows_df = cows_df[(cows_df['HCY1_c'] > 2)]

    #Filling with -999 for DMU
    realc = ['L1','L2','L3']
    cows_df[realc] = cows_df[realc].fillna(-999.0)

    #Getting code ids
    cows_df = pd.merge(left=cows_df, right=radnrkodi, on='id', how='left')
    cows_df = cows_df.sort_values(by=['code_id'])

    dmu_long = cows_df[longobs_columns]

    # DMU datafile
    dmu_long.to_csv(longobs, index=False, header=False, sep=' ')

    print(dmu_long.iloc[80000:80015])
    print(dmu_long.info())

else:
    print( f'DMU longevity file not created' )

#---------------------------------------------------------------------------
# PART 5 - Plotting TDM phenotypic data
#---------------------------------------------------------------------------

if plottdm == 1:

    tdm1file = f'../{yearmonth}/dmu_data/tdm1.dat'

    tdm1file_columns = ['code_id','hy','htd','agec','m','lp1','lp2','lp3','lp4','wil',
        'my1','my2','my3','dim','herd','calvm']
    widths_tdm1file = [8,8,8,5,5,8,9,9,9,9,8,8,8,9,8,8]

    tdm1 = readingfilefwf(tdm1file,tdm1file_columns,widths_tdm1file)

    tdm1 = pd.merge(left=tdm1, right=radnrkodi[['id','code_id']], on='code_id', how='left')

    tdm1['BY'] = (tdm1.id.astype(str).str[:4]).astype(int)

    tdm11 = tdm1.loc[(tdm1['my1'] > 0)]
    tdm12 = tdm1.loc[(tdm1['my2'] > 0)]
    tdm13 = tdm1.loc[(tdm1['my3'] > 0)]


    tdm11.loc[:, 'totalmy1'] = tdm11.groupby('code_id')['my1'].transform('sum')
    tdm12.loc[:, 'totalmy2'] = tdm12.groupby('code_id')['my2'].transform('sum')
    tdm13.loc[:, 'totalmy3'] = tdm13.groupby('code_id')['my3'].transform('sum')


    print(tdm1.iloc[800000:800015])
    print(tdm1.info())

    fig, (ax1) = plt.subplots(1, sharey=True)

    plottingmean(ax1,tdm11,'totalmy1','Mjólkurmagn','Mjaltaskeið 1')
    plottingmean(ax1,tdm12,'totalmy2','Mjólkurmagn','Mjaltaskeið 2')
    plottingmean(ax1,tdm13,'totalmy3','Mjólkurmagn','Mjaltaskeið 3')
    ax1.set_title(f'Mean sum of TD yield by lact by birth year, {yearmonth}', fontsize=10, fontweight ="bold")

    fig.set_size_inches([10, 5])
    plt.savefig(f'../{yearmonth}/figures/tdmyieldsumcheck.png')

    # plt.show()
else:
    print( f'TDM check plot not created' )

#---------------------------------------------------------------------------
# PART 6 - Accuracy fortran progarms
#---------------------------------------------------------------------------

if accuracytdm == 1:

    subprocess.call(f'gfortran -o acc1 acctdm1.f', shell=True)
    subprocess.call(f'gfortran -o acc2 acctdm2.f', shell=True)
    subprocess.call(f'gfortran -o acc2fr acctdm2_fr.f', shell=True)

    subprocess.call('./acc1', shell=True)
    subprocess.call('./acc2', shell=True)
    subprocess.call('./acc2fr', shell=True)

else:
    print( f'No accuracy calc. done' )
