#This is a program to collect EBV results for Icelandic dairy cows and
#Write to a file to be read by Huppa

#Written by Thordis Thorarinsdottir 2021

import pandas as pd
import numpy as np
import os
import shutil

#Option to read SOL files for these 3 traits and collect results, 0 means skip
# and 1 means SOL files will be read
fertility_results = 0
conf_results = 0
rankorder_results = 0
#Option to write scaled results for above traits to seperate files to be read
#later, 0 means skip and 1 means seperate result files will be written to disc
seperate_files = 0

#Option to collect results for a large datafile to be read by Huppa.
#--IF ABOVE OPTIONS ARE SET TO 0 --> program will collect results for these three
    #trait from seperate files!
#--IF ABOVE OPTIONS ARE SET TO 1 --> program will collect results for these three
    #trait straight from the program itself.
#Program always collects results for yield, scs, persistancy and longevity from
    #imported datafiles.
collectresults = 0    #0 means don't, 1 means collect

#Option to write large datafile to disc
branda = 0

collectdmufiles = 1

#---------------------------------------------------------------------------
#File with information for program!
info = pd.read_csv(
    'info',
    header=None,
    names=['info']
    )
#---------------------------------------------------------------------------
#Info for program
#---------------------------------------------------------------------------
yearmonth = info.loc[21,'info']

brandafile = info.loc[24,'info']
#DMU sol files
fertilitysolfile = '../DMU/fertilityDMU/SOL'
rankordersolfile = '../DMU/rankorderDMU/SOL'
confsolfile = '../DMU/conformationDMU/SOL'
#Observationfile used for DMU
fertilityobs = '../dmu_data/dmu_fertility.txt'
rankorderobs = '../dmu_data/dmu_rankorder.txt'
confobs = '../dmu_data/dmu_conformation.txt'
#Seperate EBV result files
tdmebv = '../results/tdmebv.txt'    #From programs by JHE
scsebv = '../results/scsebv.txt'    #From programs by JHE
perebv = '../results/perebv.txt'    #From programs by JHE
fertilityebv = '../results/fertilityebv.txt' #written by this program if option above set to 1
confebv = '../results/conformationebv.txt' #written by this program if option above set to 1
rankorderebv = '../results/rankorderebv.txt' #written by this program if option above set to 1
longgebv = '../results/ending.txt'          #from BHB
#Accuracy files for yield and scs
accyield = '../results/accuracy.sol' #From programs by JHE
accscs = '../results/accuracy_f.sol' #From programs by JHE
#Radnrkodifile
radnrkodifile = '../dmu_data/radnrkodi' #Created by prep_tdm.f
#Name of pedigree file
pedigreefile = info.loc[1,'info']
#Scaling year, present year - 5 years
scalingyear = pd.to_numeric(info.loc[3,'info'])
#Scaling objects
imean = pd.to_numeric(info.loc[4,'info'])
isd = pd.to_numeric(info.loc[5,'info'])

#---------------------------------------------------------------------------
#Columns in DMU files/own observation files
fertility_columns = ['code_id','HBY','HC1','HC2','HC3','IYM0','IYM1','IYM2',
    'IYM3','AGEi_h','AGEc_1','AGEc_2','AGEc_3','tech_h',
    'CR0','ICF1','ICF2','ICF3','IFL1','IFL2','IFL3']
conf_columns = ['code_id','HdomsY','lact','AGEc_1',
    'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenaþykkt', 'spenastada',
    'mjaltir', 'skap']
rankorder_columns = ['code_id', 'year', 'mjaltarod', 'gaedarod']
#---------------------------------------------------------------------------

#Columns in EBV files
#---------------------------------------------------------------------------
#Fixed with files
tdmebv_columns = ['id','milk_kg1','milk_kg2','milk_kg3',
    'fat_kg1','fat_kg2','fat_kg3',
    'prot_kg1','prot_kg2','prot_kg3',
    'fat_%1','fat_%2','fat_%3',
    'prot_%1','prot_%2','prot_%3','eigin_afurdir', 'yieldtotal']
widths_tdmebv = [15,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4]

scsebv_columns = ['id','scs1','scs2','scs3','scs']
widths_scsebv = [15,4,4,4,4]

perebv_columns = ['id','milkper', 'fatper', 'protper']
widths_perebv = [15,4,4,4]

longgebv_columns = ['id','longevity', 'no_daughters_longevity']
widths_longgebv = [15,4,4]

accyield_columns = ['id','no_daughters_yield', 'yield_acc']
widths_accyield = [15,8,8]

accscs_columns =  ['id','no_daughters_SCS', 'SCS_acc']
widths_accscs = [15,8,8]
#---------------------------------------------------------------------------

#Space seperated files/dataframes
#(Created by this program if option above set to 1)
fertilityebv_columns = ['id','fer_lact1','fer_lact2','fer_lact3','CR0','ICF',
    'IFL','fertility','offCR0','offICF1','offICF2','offICF3']

confebv_columns = ['id','boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti',
    'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta',
    'jugurband', 'jugurdypt', 'spenalengd', 'spenaþykkt', 'spenastada',
    'mjaltir', 'skap','offconf']

rankorderebv_columns = ['id','mjaltarod', 'gaedarod','offrankorder']

#---------------------------------------------------------------------------
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

#Format of radnrkodi file created by prep_tdm.f
radnrkodi_columns = ['id','code_id','stada','norec','fix1','fix2', 'fix3','sex']
widths_radnrkodi = [15,9,3,3,6,6,6,2]

#Format of pedigree file from Huppa
ped_columns = ['id','dam','sire','unused','sex','unused2', 'bullno', 'name', 'farm']
ped_widths = [15,15,15,12,1,2,5,20,20]

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Functions used by program!
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#Function to read sol files
#---------------------------------------------------------------------------
def readingfile(file,columns,filewidths):
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
#Function to collect phantom group results
def phantomsol(sol):
    solp = sol[(sol['1_code_effect'] == 4) & (sol['code_id'] < 0)]
    return solp

#---------------------------------------------------------------------------
#Function to find solution for each trait in sol dataframes and merge them by ids
#into a new df
def solutions(df, df1, no, trait, id):
    print( f'----------------------------------------' )
    print( f'Collecting results for {trait} .....' )
    print( f'----------------------------------------' )
    #Only keep genetic effectss
    sol = df[(df['1_code_effect'] == 4) & (df['code_id'] > 0)] #collect animal results
    # merge real ids back into results
    sol= pd.merge(left=sol[
        ['code_id','2_trait_no', '6_no_obs','8_BLUP']
        ], right=radnrkodi[['id','code_id']], on='code_id', how='left')

    dftrait = sol.loc[sol['2_trait_no'] == no ] #one dataframe per trait
    dftrait.loc[:,trait] = dftrait['8_BLUP']            #name the trait
    df1 = pd.merge(left=df1, right=dftrait[[id,trait]  #merge into large dataframe
        ], on=id)
    df1 = df1.drop_duplicates(subset=['id'])
    df1 = df1.sort_values(by=['id'])
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
#-------------------------------------------------------------
#-------------------------------------------------------------
#START OF PROGRAM!!!
#-------------------------------------------------------------
#-------------------------------------------------------------
#Reading in id codes to replace in SOL files
if fertility_results == 1 | conf_results == 1 | rankorder_results == 1 :
    radnrkodi = readingfile(radnrkodifile,radnrkodi_columns,widths_radnrkodi)
    radnrkodi = radnrkodi.drop(
        ['stada','norec','fix1','fix2', 'fix3','sex'], axis = 1)

#Reading fertilty SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if fertility_results == 1:
    #Reading sol files
    fertilitysol = readingfile(fertilitysolfile,solcolumns,sol_widths)
    #Phantom group results
    fertilityph = phantomsol(fertilitysol)
    #Seperating solutions by traits
    fertilitydf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'CR0','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'ICF1','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'ICF2','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'ICF3','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'IFL1','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'IFL2','id')
    fertilitydf = solutions(fertilitysol, fertilitydf, 1, 'IFL3','id')
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
    #Reading pedigree
    ped = readingfile(pedigreefile,ped_columns,ped_widths)
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

else:
    print('Fertility results not collected')

#Reading conformation SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if conf_results == 1:
    #Reading sol files
    confsol = readingfile(confsolfile,solcolumns,sol_widths)
    #Phantom group results
    confph = phantomsol(confsol)
    #Seperating solutions by traits
    confdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results
    confdf = solutions(confsol, confdf, 1, 'boldypt','id')
    confdf = solutions(confsol, confdf, 2, 'utlogur','id')
    confdf = solutions(confsol, confdf, 3, 'yfirlina','id')
    confdf = solutions(confsol, confdf, 4, 'malabreidd','id')
    confdf = solutions(confsol, confdf, 5, 'malahalli','id')
    confdf = solutions(confsol, confdf, 6, 'malabratti','id')
    confdf = solutions(confsol, confdf, 7, 'stada_haekla_hlid','id')
    confdf = solutions(confsol, confdf, 8, 'stada_haekla_aftan','id')
    confdf = solutions(confsol, confdf, 9, 'klaufhalli','id')
    confdf = solutions(confsol, confdf, 10, 'jugurfesta','id')
    confdf = solutions(confsol, confdf, 11, 'jugurband','id')
    confdf = solutions(confsol, confdf, 12, 'jugurdypt','id')
    confdf = solutions(confsol, confdf, 13, 'spenalengd','id')
    confdf = solutions(confsol, confdf, 14, 'spenaþykkt','id')
    confdf = solutions(confsol, confdf, 15, 'spenastada','id')
    confdf = solutions(confsol, confdf, 16, 'mjaltir','id')
    confdf = solutions(confsol, confdf, 17, 'skap','id')
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
    confdf['spenaþykkt'] = scaling(confdf,'spenaþykkt',1, conf_ave)
    confdf['spenastada'] = scaling(confdf,'spenastada',1, conf_ave)
    confdf['mjaltir'] = scaling(confdf,'mjaltir',1, conf_ave)
    confdf['skap'] = scaling(confdf,'skap',1, conf_ave)

    #Counting daughters with own obs
    #Reading pedigree
    ped = readingfile(pedigreefile,ped_columns,ped_widths)
    ownobs_conf = pd.merge(left=ownobs_conf, right=ped[['id','sire']], on='id', how='left')

    confdf = countingoff(ownobs_conf,'boldypt','offconf', confdf)

    if seperate_files == 1:
        print('Conformation results written to seperate file')
        conf_results = confdf[confebv_columns].astype(float).round(2)#.astype(int)
        conf_results.to_csv(confebv, index=False, header=False, sep=' ')
        print(f'Conformation results written to {confebv}')

        print(confdf.iloc[50000:50015])
        print(confdf.info())
else:
    print('Conformation results not collected')

#Reading rankorder SOL file, collecting and scaling results and counting daughters
#Write results to disc if option above set to 1
if rankorder_results == 1:
    #Reading sol files
    rankordersol = readingfile(rankordersolfile,solcolumns,sol_widths)
    #Phantom group results
    rankorderph = phantomsol(rankordersol)
    #Seperating solutions by traits
    rankorderdf = radnrkodi['id'].copy()  #Creating a dataframe to merge trait results

    rankorderdf = solutions(rankordersol, rankorderdf, 1, 'mjaltarod','id')
    rankorderdf = solutions(rankordersol, rankorderdf, 1, 'gaedarod','id')
    #Reading own observations files
    ownobs_rankorder = ownobs(rankorderobs,rankorder_columns)
    #Creating average groups to scale
    rankorder_ave = avegroup(rankorderdf,ownobs_rankorder)
    #Scaling results
    rankorderdf['mjaltarod'] = scaling(rankorderdf,'mjaltarod',-1, rankorder_ave)
    rankorderdf['gaedarod'] = scaling(rankorderdf,'gaedarod',-1, rankorder_ave)

    #Counting daughters with own obs
    #Reading pedigree
    ped = readingfile(pedigreefile,ped_columns,ped_widths)
    ownobs_rankorder = pd.merge(left=ownobs_rankorder, right=ped[['id','sire']], on='id', how='left')

    rankorderdf = countingoff(ownobs_rankorder,'mjaltarod','offrankorder', rankorderdf)

    if seperate_files == 1:
        print('Rankorder results written to seperate file')
        rank_results = rankorderdf[rankorderebv_columns].astype(float).round().astype(int)
        rank_results.to_csv(rankorderebv, index=False, header=False, sep=' ')
        print(f'Rankorder results written to {rankorderebv}')

    print(rankorderdf.iloc[50000:50015])
    print(rankorderdf.info())
else:
    print('Rank order results not collected')


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

    results = combineresultsfwf(results,tdmebv,tdmebv_columns,widths_tdmebv) #test day
    results = combineresultsfwf(results,scsebv,scsebv_columns,widths_scsebv) #somatic cell score
    results = combineresultsfwf(results,perebv,perebv_columns,widths_perebv) #persistancy
    results = combineresultsfwf(results,longgebv,longgebv_columns,widths_longgebv)    #longevity
    results = combineresultsfwf(results,accyield,accyield_columns,widths_accyield) #accuracy yield
    results = combineresultsfwf(results,accscs,accscs_columns,widths_accscs) #accuracy scs

    results = results[(results['milk_kg1'].notnull().astype(int) == 1)] #only results for animals in yield file


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


    #Combination of three lactations to make one trait
    results['milk_kgt'] = (results['milk_kg1'] * 0.5 +
                            results['milk_kg2'] * 0.3 +
                            results['milk_kg3'] * 0.2 )

    results['fat_kgt'] = (results['fat_kg1'] * 0.5 +
                            results['fat_kg2'] * 0.3 +
                            results['fat_kg3'] * 0.2 )

    results['prot_kgt'] = (results['prot_kg1'] * 0.5 +
                            results['prot_kg2'] * 0.3 +
                            results['prot_kg3'] * 0.2 )

    results['fat_%t'] = (results['fat_%1'] * 0.5 +
                            results['fat_%2'] * 0.3 +
                            results['fat_%3'] * 0.2 )

    results['prot_%t'] = (results['prot_%1'] * 0.5 +
                            results['prot_%2'] * 0.3 +
                            results['prot_%3'] * 0.2 )

#Creation of a total grade for udder
    results['jugur'] = (results['jugurfesta'] * 0.35 +
                            results['jugurband'] * 0.15 +
                            results['jugurdypt'] * 0.5 )

#Creation of a total grade for tits
    results['spenar'] = (results['spenalengd'] * 0.3 +
                            (200 - results['spenaþykkt']) * 0.3 +
                            results['spenastada'] * 0.4 )

#Creation of a total grade for milking
    results['mjaltir_t'] = (results['mjaltir'] * 0.6 +
                            results['mjaltarod'] * 0.4 )

#Creation of a total grade for bulls with grade for longevity
    results.loc[(results['no_daughters_longevity'] != 0), 'total'] = (
                    results['yieldtotal']*0.36 +
                    results['fertility']*0.10 +
                    results['scs']*0.08 +
                    results['jugur']*0.10 +
                    results['spenar']*0.10 +
                    results['mjaltir_t']*0.08 +
                    results['skap']*0.08 +
                    results['longevity']*0.10)

#Creation of a total grade for cows
    results.loc[(results['no_daughters_longevity'] == 0), 'total'] = (
                    results['yieldtotal']*0.36 +
                    results['fertility']*0.11 +
                    results['scs']*0.09 +
                    results['jugur']*0.11 +
                    results['spenar']*0.13 +
                    results['mjaltir_t']*0.10 +
                    results['skap']*0.10 +
                    results['longevity']*0.0)

#Scaled EBVs written to disc
    if branda == 1:

        branda = results[['id',                                             #1
            'milk_kg1','milk_kg2','milk_kg3','fat_kg1','fat_kg2','fat_kg3', #6
            'prot_kg1','prot_kg2','prot_kg3','fat_%1','fat_%2','fat_%3',    #6
            'prot_%1','prot_%2','prot_%3','scs1','scs2','scs3',             #6
            'milkper', 'fatper', 'protper',                                 #3
            'fer_lact1','fer_lact2','fer_lact3','CR0','ICF','IFL',          #6
            'boldypt', 'utlogur', 'yfirlina', 'malabreidd', 'malahalli', 'malabratti', #6
            'stada_haekla_hlid', 'stada_haekla_aftan', 'klaufhalli', 'jugurfesta', #4
            'jugurband', 'jugurdypt', 'spenalengd', 'spenaþykkt', 'spenastada',  #5
            'mjaltir', 'skap','mjaltarod', 'gaedarod','longevity',              #5
            'milk_kgt','fat_kgt','prot_kgt','fat_%t','prot_%t',                 #5
            'yieldtotal','fertility','scs','jugur','spenar','mjaltir_t','skap', #7
            'total',                                                             #1
            'no_daughters_yield', 'yield_acc', 'no_daughters_SCS', 'SCS_acc', #4
            'offconf','offrankorder','offCR0','offICF1','offICF2','offICF3', #6
            'no_daughters_longevity'                                        #1
            ]].astype(int)                                            #72 columns


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
#
