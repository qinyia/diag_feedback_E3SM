
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
dic['tas_latlon'] = {}
dic['LCF'] = {}
dic['zm_CLOUD'] = {}
dic['latlon_CLOUD'] = {}
dic['webb_decomp'] = {}

for key in dic.keys():
    dic[key] = pd.DataFrame()
 
    fname = 'csvfile/pd2html_'+key+'.csv'
    if os.path.isfile(fname):
        dic[key] = pd.read_csv(fname,index_col=0)

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
  <h3>Zonal Mean Radiative Kernel Feedback</h3>
    {table4}
  <h3>Zonal Mean Cloud Radiative Kernel Feedback</h3>
    {table5}
  <h3>LAT-LON Radiative Kernel Feedback</h3>
    {table6}
  <h3>LAT-LON Cloud Radiative Kernel Feedback</h3>
    {table7}
  <h3>LAT-LON Surface Air Temperature Response</h3>
    {table8}
  <h3>Liquid Condensate Fraction</h3>
    {table9}
  <h3>Zonal Mean Cloud Related Vars</h3>
    {table10}
  <h3>LAT-LON Cloud Related Vars</h3>
    {table11}
  <h3>LAT-LON Cloud Feedbacks from Webb Method</h3>
    {table12}

  </div>
  </body>
</html>.
'''

with open('viewer/index.html','w') as f:
    f.write(html_string.format(table1=dic['CRE_globalmean'].to_html(escape=False,index=False,border=0),
                               table2=dic['RadKernel_globalmean'].to_html(escape=False,index=False,border=0),
                               table3=dic['CldRadKernel_globalmean'].to_html(escape=False,index=False,border=0),
                               table4=dic['RadKernel_zonalmean'].to_html(escape=False,index=False,border=0),
                               table5=dic['CldRadKernel_zonalmean'].to_html(escape=False,index=False,border=0),
                               table6=dic['RadKernel_latlon'].to_html(escape=False,index=False,border=0),
                               table7=dic['CldRadKernel_latlon'].to_html(escape=False,index=False,border=0),
                               table8=dic['tas_latlon'].to_html(escape=False,index=False,border=0),
                               table9=dic['LCF'].to_html(escape=False,index=False,border=0),
                               table10=dic['zm_CLOUD'].to_html(escape=False,index=False,border=0),

                               table11=dic['latlon_CLOUD'].to_html(escape=False,index=False,border=0),
                               table12=dic['webb_decomp'].to_html(escape=False,index=False,border=0),
                               )
                               )


print("=============Well done. check viewer/index.html for all plots.==============")

bashCommand = "tar cvf tar.tar viewer figure"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()

bashCommand = "cp tar.tar /compyfs/www/qiny108/diag_feedback/"
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output, error = process.communicate()




