import sys
import numpy as np
import pandas as pd
import requests
from pathlib import Path
import datetime as dt

varlib = {'Temp':0,'Par':0,'SurfacePar':0,'Density':0,'SigmaDensity':0,'Do':0,'Chla':0,'Sal':0,'Turb':0,'Light':0,'Nitrate':0}

locs = pd.read_csv('Water_Quality_Monitoring_Dashboard_-_Monitoring_Sites.csv')

Locators = ["Adm Inlet-1","Adm Inlet-2","Adm Inlet-3","Adm Inlet-4-C14","Adm Inlet-5","Alki#1","Alki#2","Alki#3","Alki#4","CK200P","Colvos Pass","Edmds-1","Edmds-2","Edmds-3","Edmds-4","Edmds-5","Edmds-6","Gedney Is","HNFD01","JSTU01","JSUR01","KSBP01","KSIW02","KSRU03","KSSK02","LSEP01","LSKQ06","LSNT01","LSVV01","LTBC41","LTBC42","LTBC43","LTED04","LTKE03","LTLF04","LTUM03","MSJN02","NSEX01","Orca Buoy","Poss DO-1","Poss DO-2","Poss DO-4","Poss Snd-1","Poss Snd-2","Poss Snd-3-C14","Poss Snd-4","Poss Snd-5","Pt Susan","Pt Wells-1-C14","Pt Wells-2","Pt Wells-3","Pt Wells-4","Pt Wells-5","VO50E","WP#1","WP#2","WP#3","WP#4","Z5-DIFF-CTD","Z6SBE-DIFF-CTD","Z6SBW-DIFF-CTD","Z7SNB-DIFF-CTD","Z7SSA-DIFF-CTD","Z7SSBE-DIFF-CTD","Z7SSBW-DIFF-CTD"]

def data_pugetsound(config):
    '''Example how to use:
    python -m data_pugetsound year var1 var2 var3 var4...'''
    output_dir = Path("Puget_DO_obs")
    output_dir.mkdir(exist_ok=True)
    startdate = dt.datetime(int(config[0]),1,1)
    enddate = dt.datetime(int(config[0]),12,31)
    varlist = config[1:]#['Do','Nitrate']

    for var in varlist:
        varlib[var]+=1

    modurl = f"&StartDate={startdate.month}%2F{startdate.day}%2F{startdate.year}&EndDate={enddate.month}%2F{enddate.day}%2F{enddate.year}&Temp={str(varlib['Temp']>0).lower()}&Par={str(varlib['Par']>0).lower()}&SurfacePar={str(varlib['SurfacePar']>0).lower()}&Density={str(varlib['Density']>0).lower()}&SigmaDensity={str(varlib['SigmaDensity']>0).lower()}&Do={str(varlib['Do']>0).lower()}&Chla={str(varlib['Chla']>0).lower()}&Sal={str(varlib['Sal']>0).lower()}&Turb={str(varlib['Turb']>0).lower()}&Light={str(varlib['Light']>0).lower()}&Nitrate={str(varlib['Nitrate']>0).lower()}"

    for loc in Locators:
        url = f"https://green2.kingcounty.gov/marine/Download/DownloadCTDData?Locator={loc}"+modurl

        filename = loc+"_2018.csv"
        output_file = output_dir / filename
        response = requests.get(url)
        if len(response.content)>3000:
            output_file.write_bytes(response.content)
            print(f"Downloaded station {loc}")

    df =pd.DataFrame(columns=['Locator','Sample_Date','Sample_Depth','UpDown','DO_field','DO_Qual','NO23_field','NO23_Qual','Latitude','Longitude'])

    for csv_file in output_dir.glob("*.csv"):
        dfn = pd.read_csv(csv_file,header=[1],encoding='utf-16', sep=None,engine="python")
        dfn['Latitude'] = locs[locs['Locator']==dfn['Locator'][0]].Latitude.values[0]*np.ones(len(dfn))
        dfn['Longitude'] = locs[locs['Locator']==dfn['Locator'][0]].Latitude.values[0]*np.ones(len(dfn))
        df = pd.concat([df,dfn])

    outfile = f'PugetSound_{startdate.year}.csv'
    df.to_csv(outfile)
    return print(f'Data saved to'+ outfile)


if __name__=="__main__":
    try:
        config = sys.argv[1:]
    except :
        print('Something went wrong')
    data_pugetsound(config)