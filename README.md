# Here is the documentation for the new pipeline in the calculation of EBVs of Icelandic dairy cows

New methods and this new pipeline are in development. 


## Older programs in new pipeline
### prep_tdm.f 
  - Collects data about yield and scs and prepares for DMU
  - Creates radnrkodi for code numbers to be used in pipeline
  - Creates pedigree for all traits for DMU
  - JHE

### solsaman.f
  - ????
  - JHE

### acctdm1.f & acctdm1_f.f
  - Calculates accuracy for yield and scs

## New programs in new pipeline

### 1fromhuppa.py 
  - collects data from Huppa about fertility, conformation traits and rank order and prepares for DMU
  - ÞÞ
### 2branda.py 
  - collects and scales results from DMU SOL files about fertility, conformation traits and rank order. 
  - Collects results about yield, scs, persistancy, accuracy and longevity written by other programs.
  - Creates a datafile to be read into Huppa
  - ÞÞ
  
  
