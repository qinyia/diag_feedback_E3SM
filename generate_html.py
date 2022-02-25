

def generate_html(casedir,webdir=None):

    '''
    Dec 22, 2021: generate html page to visualize all plots.

    input:
    casedir -- the directory where related files are saved.
    webdir -- the directory where you can open the html file.
              Both compy and nersc have such web server. 
              You can also download the final tar file in the 
              current directory to your local computer. 
    '''
    import pandas as pd
    import os
    
    import subprocess
    
    dic = {}
    dic['CRE_globalmean'] = {}
    dic['RadKernel_globalmean'] = {}
    dic['CldRadKernel_globalmean'] = {}
    dic['RadKernel_zonalmean'] = {}
    dic['CldRadKernel_zonalmean'] = {}
    dic['RadKernel_latlon'] = {}
    dic['CldRadKernel_latlon'] = {}
    dic['RadKernel_latlon_dif'] = {}
    dic['CldRadKernel_latlon_dif'] = {}
    dic['tas_latlon'] = {}
    dic['LCF'] = {}
    dic['zm_CLOUD'] = {}
    dic['latlon_CLOUD'] = {}
    dic['webb_decomp'] = {}
    dic['CLOUD_profile'] = {}
    dic['NRMSE_RadKern'] = {}
    dic['cal_regionCor'] = {}

    
    for key in dic.keys():
        dic[key] = pd.DataFrame()
     
        fname = casedir+'csvfile/pd2html_'+key+'.csv'
        if os.path.isfile(fname):
            dic[key] = pd.read_csv(fname,index_col=0)

            # if the dict is empty, don't show the empty table only with columns.
            if len(dic[key]) == 0:
                dic[key] = pd.DataFrame()
    
    # convert into html files 
    html_string = '''
    <!DOCTYPE html>
    <html>
      <head>
      <title>diag_cloud_feedback</title>
      <meta name="viewport" content="width=device-width, initial-scale=1.0">
      <meta charset="utf-8">
      <link rel="stylesheet" href="11.css" type="text/css">
      </head>
      <body>
      <div class="container">
      <div class="row">
      <h1>Cloud Feedback Diagnostics</h1>
      <hr>
      </div>
      <h3>Global Mean Radiative Feedback</h3>
       {table1}
      <h3>Global Mean Radiative Kernel Feedback</h3>
       {table2}
      <h3>Global Mean Cloud Radiative Kernel Feedback</h3>
        {table3}
      <h3>LAT-LON Radiative Kernel Feedback</h3>
        {table13}
      <h3>LAT-LON Cloud Radiative Kernel Feedback</h3>
        {table14}
      <h3>Zonal Mean Radiative Kernel Feedback</h3>
        {table4}
      <h3>Zonal Mean Cloud Radiative Kernel Feedback</h3>
        {table5}
      <h3>Liquid Condensate Fraction</h3>
        {table9}
      <h3>LAT-LON Radiative Kernel Feedback Difference</h3>
        {table6}
      <h3>LAT-LON Cloud Radiative Kernel Feedback Difference</h3>
        {table7}
      <h3>LAT-LON Surface Air Temperature Response Difference</h3>
        {table8}
      <h3>Zonal Mean Cloud Related Variables Difference</h3>
        {table10}
      <h3>LAT-LON Cloud Related Variables Difference</h3>
        {table11}
      <h3>LAT-LON Cloud Feedbacks from Webb Method</h3>
        {table12}
       <h3>Regional Mean Cloud Profiles</h3>
        {table15}
        <h3>NRMSE and COR evolution for RadKern</h3>
        {table16}
        <h3>regional NRMSE and COR evolution for RadKern</h3>
        {table17}
 
      </div>
      </body>
    </html>.
    '''
    
    with open(casedir+'viewer/index.html','w') as f:
        f.write(html_string.format(table1=dic['CRE_globalmean'].to_html(escape=False,index=False,border=0),
                                   table2=dic['RadKernel_globalmean'].to_html(escape=False,index=False,border=0),
                                   table3=dic['CldRadKernel_globalmean'].to_html(escape=False,index=False,border=0),
                                   table4=dic['RadKernel_zonalmean'].to_html(escape=False,index=False,border=0),
                                   table5=dic['CldRadKernel_zonalmean'].to_html(escape=False,index=False,border=0),
                                   table13=dic['RadKernel_latlon'].to_html(escape=False,index=False,border=0),
                                   table14=dic['CldRadKernel_latlon'].to_html(escape=False,index=False,border=0),
                                   table6=dic['RadKernel_latlon_dif'].to_html(escape=False,index=False,border=0),
                                   table7=dic['CldRadKernel_latlon_dif'].to_html(escape=False,index=False,border=0),
                                   table8=dic['tas_latlon'].to_html(escape=False,index=False,border=0),
                                   table9=dic['LCF'].to_html(escape=False,index=False,border=0),
                                   table10=dic['zm_CLOUD'].to_html(escape=False,index=False,border=0),
    
                                   table11=dic['latlon_CLOUD'].to_html(escape=False,index=False,border=0),
                                   table12=dic['webb_decomp'].to_html(escape=False,index=False,border=0),
                                   table15=dic['CLOUD_profile'].to_html(escape=False,index=False,border=0),
                                   table16=dic['NRMSE_RadKern'].to_html(escape=False,index=False,border=0),
                                   table17=dic['cal_regionCor'].to_html(escape=False,index=False,border=0),
                                   )
                                   )
    
    
    print("=============Well done. check viewer/index.html for all plots.==============")
    
    bashCommand = 'sh copy2webdir.sh '+casedir+' '+webdir
    #process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
    #output, error = process.communicate()

    os.system(bashCommand)
    
    print("========You can check your web page here: ",webdir+"====================")

