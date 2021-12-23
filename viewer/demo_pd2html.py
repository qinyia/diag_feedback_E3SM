#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[3]:


df = pd.DataFrame({"Name": ['Anurag', 'Manjeet', 'Shubham', 
                            'Saurabh', 'Ujjawal'],
                     
                   "Address": ['Patna', 'Delhi', 'Coimbatore',
                               'Greater noida', 'Patna'],
                     
                   "ID": [20123, 20124, 20145, 20146, 20147],
                     
                   "Sell": [140000, 300000, 600000, 200000, 600000]})


# In[4]:


df['Sell'] = df['Sell'].apply(lambda x: '<a href="18.jpeg">plot</a>')
df


# In[12]:


html_string = '''
<html>
  <head><title>HTML Pandas Dataframe with CSS</title></head>
  <link rel="stylesheet" type="text/css" href="11.css"/>
  <body>
      <h1>hello</h1>
    {table}
      <h1>world</h1> 
    {table1}
  </body>
</html>.
'''


# In[13]:


with open('myhtml.html', 'w') as f:
    f.write(html_string.format(table=df.to_html(escape=False, index=False), 
                               table1=df.to_html(escape=False, index=False)))


# In[ ]:




