# optimal_kinase_library
Library of kinase datasets commonly used for modeling tasks. For data/results folders, see [OneDrive](https://adminliveunc-my.sharepoint.com/:f:/r/personal/smgomez_ad_unc_edu/Documents/GomezLabDataShare/GithubData/optimal_kinase_library?csf=1&web=1&e=WdzFzN) (may require permission)

## File Descriptions (WIP)
> Main KIS data sources:
>* Klaeger/Kinomescan
>   * %inhibition: various conc levels
>   * Kd inhibition: dissociation constant/binding affinity, gives a holistic view of drug effect (vs at specific conc levels)
>* Kuster
>   * extra ~800 compounds that are a follow up to original ~200 Klaeger compounds
>* OKL (Optimal Kinase Library)
>   * Peter Sorger's inhibition data

```data/```
* raw data, including ALMANAC, CCLE, Klaeger, LINCS, Kuster, OKL, PDX, and PRISM

```figures/```
* ???

```results/```
* preprocessed and combined datasets
* ```combined_dose_data.fst```
    * pooled data from all sources, with source column for identification/filtering
    * KIS range from 0-1 indicating 0-100% inhibition
* ```combined_kd_data.fst```
    * Kd data (when provided/calculated)
* ```XXX_summarized.fst```
    * summarized data across sources where the mean is used if there is exact overlap between drug AND conc between sources

```src/```
* R scripts for dataset preprocessing

## Notes
* Default to Klaeger + Kuster dataset as our "best" single source, as they are most consistent.
* ```data/PRISM/``` contains ```primary``` and ```secondary``` folders. These correspond to two rounds of expirements/screens. We only use the ```secondary``` screen.