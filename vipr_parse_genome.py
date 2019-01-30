import pandas as pd
from dateutil.parser import parse
import re, os
import shutil


if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='parse metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--input', help="input meta file")
    parser.add_argument('--output', help="output meta")
    parser.add_argument('--regions', help="file to specify regions")
    args = parser.parse_args()

    input_gb_meta = args.input #"data/allEntero-3Oct18.tsv"
    output_gb_meta = args.output #"data/20181003_GenBank.tsv"

    #list countries that need replacing
    # I don't know why Denmark matching doesn't work in lowercase - it should
    country_replace = [['Viet Nam', 'vietnam'], ['Denmark', 'denmark'], [" ", "_"]]

    #Delete the last tab from the header row of the file downloaded from vipr... messes up pandas!!
    #copy - so we can overwrite original
    shutil.copyfile(input_gb_meta, input_gb_meta+"temp")

    from_file = open(input_gb_meta+"temp")
    line = from_file.readline()
    line = line.replace("\t\n", "\n")
    to_file = open(input_gb_meta, mode="w")
    to_file.write(line)
    shutil.copyfileobj(from_file, to_file)
    to_file.close()
    from_file.close()
    os.remove(input_gb_meta+"temp")

    meta = pd.read_csv(input_gb_meta, sep='\t', index_col=False)

    meta.columns = ['strain', 'virus', 'accession', 'seq_len', 'date', 'host', 'genbank_host', 'country', 'moltype']

    #I am not sure where subgenogroup came from?

    meta = meta.reindex(columns = ['strain', 'virus', 'subgenogroup', 'accession', 'seq_len', 'date', 'age', 'host', 'genbank_host', 'country', 'moltype'])

    #replace -N/A- with date format that works in augur
    meta.loc[meta.date == '-N/A-', 'date'] = '20XX-XX-XX'

    #this assumes vipr is consistant in how they do dates...
    newDate = []
    for d in meta.date:
        if 'XX' not in d:
            if d.count('/') == 0:
                augur_date = d+"-XX-XX"
            elif d.count('/') == 1:
                split_date = d.split("/")
                augur_date = split_date[1]+"-"+split_date[0]+"-XX"
            else:
                parsed_date = parse(d)
                augur_date = parsed_date.strftime("%Y-%m-%d")
        else:
            augur_date = d
        newDate.append(augur_date)
    meta.loc[:, 'date'] = newDate

    #Change host 'Unknown' to NA, and '-NA-' to NA
    meta.loc[meta.host == 'Unknown', 'host'] = 'NA'
    meta.loc[meta.genbank_host == '-N/A-', 'genbank_host'] = 'NA'
    meta.loc[meta.country == '-N/A-', 'country'] = 'NA'

    #get age in months
    """
    Age formats possible:
    4M  0.75Y   3.12Y   --> re.search('[0-9.]+[A-Z]', host)    m.group(0)

    age 5   --> re.search('age [0-9]+$', host)

    4 years 1 year  3.4 years   --> re.search('[0-9.]+ year', host)

    4 months    1 months    1 month  --> re.search('[0-9.]+ month', host)
    """
    newAge = []
    for host in meta.genbank_host:
        if re.search('[0-9.]+[A-Z]', host): 
            m = re.search('[0-9.]+[A-Z]', host).group(0)
            if 'Y' in m:
                age_month = round(12*float(m.replace("Y", "")), 2)
            else:
                age_month = float(m.replace("M", ""))

        elif re.search('age [0-9]+$', host): # these are always year
            m = re.search('age [0-9]+$', host).group(0)
            age_month = round(12*float(m.replace("age ","")), 2)

        elif re.search('[0-9.]+ year', host): #these are always year
            m = re.search('[0-9.]+ year', host).group(0)
            age_month = round(12*float(m.replace(" year","")), 2)

        elif re.search('[0-9.]+ month', host):
            m = re.search('[0-9.]+ month', host).group(0)
            age_month = float(m.replace(" month",""))

        else:
            age_month = ''

        newAge.append(age_month)
    meta.loc[:, 'age'] = newAge

    #replace countries and remove spaces
    for i, row in meta.iterrows():
        coun = row.country
        for c,r in country_replace:
            coun = coun.replace(c,r)
        meta.at[i, 'country'] = coun

    #get region if supplied
    if args.regions:
        regions = {}
        with open(args.regions) as f:
            regs = f.readlines()
        regs = regs[1:] #remove first line

        for x in regs:
            pair = x.split()
            if len(pair) is not 0:
                regions[pair[0]] = pair[1]

        meta = meta.reindex(columns = ['strain', 'virus', 'subgenogroup', 'accession', 'seq_len', 'date', 'age', 'host', 'genbank_host', 'country', 'region', 'moltype'])
        newregion = []
        for coun in meta.country:
            coun = coun.lower()
            reg = "NA"
            if coun != "na":
                if coun not in regions:
                    print("No region found for {}! Setting to NA".format(coun))
                else:
                    reg = regions[coun]
            newregion.append(reg)

        meta.loc[:, 'region'] = newregion

    #correct column names to match other meta file
    meta.rename(columns={'seq_len':'seq-len', 'genbank_host':'genbank-host'}, inplace=True)
    meta.to_csv(output_gb_meta, sep='\t', index=False)


