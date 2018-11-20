'''
To do the web scrapping of FRB website for their transcript and Beige books
for the FRB textual project
'''

from collections import Counter, defaultdict
import re
from nltk.tokenize import word_tokenize
import pandas as pd
import numpy as np
import scipy
from decimal import getcontext
from decimal import Decimal as Dec
from bs4 import BeautifulSoup
import requests,re,codecs
import urllib3

def download_file(download_url, file_name):
    http = urllib3.PoolManager()
    response = http.request('GET', download_url)
    with open(file_name, 'wb') as f:
        f.write(response.data)
    response.release_conn()

def refine_text(text_raw, pattern_head=None, pattern_tail=None):
    try:
        position_head = re.search(pattern_head, text_raw).span()[1]
    except:
        position_head = 0
    try:
        position_tail = re.search(pattern_tail, text_raw).span()[0]
    except:
        position_tail = len(text_raw)
    text = text_raw[position_head:position_tail]
    return text

def get_transcript(transcript_link, beige_book_date, meeting_date, year):
    try:
        file_name = "Transcript" + transcript_link[-23:]
        download_file(transcript_link, file_name)
        print('%d, %d Transcript finished' % (year, meeting_date))
    except:
        print("Transcript error %d" % meeting_date)

def get_old_beige_book_text(page_link, div_id, beige_book_date, year):
    beige_book_page = requests.get(page_link)
    beige_book_page_soup = BeautifulSoup(beige_book_page.content, 'html.parser')
    try:
        text_raw = beige_book_page_soup.find_all('table')[0].find_all('tr')[1].find_all('td')[2].text
    except:
        text = ""
        print("Beige book error %d, %s" % (beige_book_date, div_id))
        return text
    pattern_head = re.compile(r'(is not a commentary on the views of Federal Reserve officials\.)', re.IGNORECASE)
    pattern_tail = re.compile(r'(Return to top)', re.IGNORECASE)
    text = refine_text(text_raw, pattern_head, pattern_tail)
    return text
    
def get_new_beige_book_text(page_link, div_id, beige_book_date, year):
    beige_book_page = requests.get(page_link)
    beige_book_page_soup = BeautifulSoup(beige_book_page.content, 'html.parser')
    beige_book_page_content = beige_book_page_content = beige_book_page_soup.find(id="leftText")
    try:   
        text_raw = beige_book_page_content.find(id=div_id).text
    except:
        text = ""
        print("Beige book error %d, %s" % (beige_book_date, div_id))
        return text
    pattern_head = re.compile(r'(is not a commentary on the views of Federal Reserve officials\.)', re.IGNORECASE)
    pattern_tail = re.compile(r'(Return to top)', re.IGNORECASE)
    text = refine_text(text_raw, pattern_head, pattern_tail)
    return text

def get_beige_book(beige_book_link, df_beige_books, beige_book_date, meeting_date, year):
    div_id_list = ['div_summary','div_boston','div_new_york','div_philadelphia',\
                   'div_cleveland','div_richmond','div_atlanta','div_chicago',\
                   'div_st_louis','div_minneapolis','div_kansas_city',\
                   'div_dallas','div_san_francisco']
    beige_book_dict = {}
    beige_book_dict['Year'], beige_book_dict['Meeting_Date'], beige_book_dict['Beige_Book_Date']\
     = [year], [meeting_date], [beige_book_date]
    if year < 2011:
        beige_book_dict['div_summary'] = [get_old_beige_book_text(beige_book_link, 'div_summary', beige_book_date, year)]
        for i in range(1,13):
            district_link = beige_book_link[0:-11] + str(i) + '.htm'
            beige_book_dict[div_id_list[i]] = [get_old_beige_book_text(district_link, div_id_list[i], beige_book_date, year)]
    else:
        for i in range(0,13):
            beige_book_dict[div_id_list[i]] = [get_new_beige_book_text(beige_book_link, div_id_list[i], beige_book_date, year)]
    beige_book_df = pd.DataFrame(beige_book_dict)
#    df_beige_books = df_beige_books.append(beige_book_df)
    print('%d, %d Beige Book finished' % (year, beige_book_date))
    return df_beige_books.append(beige_book_df)
    # signal printed in subfunctions
    
def get_statement(statement_link, df_statements, beige_book_date, meeting_date, year):
    try:
        statement_page = requests.get(statement_link)
        statement_page_soup = BeautifulSoup(statement_page.content, 'html.parser')        
        if year < 2006:
            text_raw = statement_page_soup.find_all('table')[-1].text
        else:
            text_raw = statement_page_soup.find(id="article").\
             find(class_="col-xs-12 col-sm-8 col-md-8").text
        
        pattern_tail = re.compile(r'(Voting for the FOMC monetary policy action were|\
                                [0-9]{4} Monetary Policy)', re.IGNORECASE)
        text = refine_text(text_raw, None, pattern_tail)
        statement_df = pd.DataFrame({'Year':[year], 'Meeting_Date':[meeting_date],\
                                     'Beige_Book_Date':[beige_book_date], 'Statement':[text]})
    except:
        statement_df = pd.DataFrame({'Year':[year], 'Meeting_Date':[meeting_date],\
                                     'Beige_Book_Date':[beige_book_date], 'Statement':[""]})
        print("Statement error %d" % meeting_date)
    print('%d, %d Statement finished' % (year, meeting_date))
#    df_statements = df_statements.append(statement_df)
    return df_statements.append(statement_df)

def get_materials(df_links, df_beige_books, df_statements):
    for row in df_links.itertuples():
        year, meeting_date	, beige_book_date, beige_book_link, statement_link, transcript_link\
         = row.Year, row.Meeting_Date, row.Beige_Book_Date, row.Beige_Book_Link, row.Statement_Link, row.Transcript_Link
        get_transcript(transcript_link, beige_book_date, meeting_date, year)
#        df_beige_books = get_beige_book(beige_book_link, df_beige_books, beige_book_date, meeting_date, year)
        df_statements = get_statement(statement_link, df_statements, beige_book_date, meeting_date, year)
    return df_beige_books, df_statements

def get_links(year_start, year_end):
    df_links = pd.DataFrame({'Year':[], 'Meeting_Date':[], 'Beige_Book_Date':[],\
                             'Beige_Book_Link':[], 'Statement_Link':[], 'Transcript_Link':[]})
    for year in range(year_start, year_end+1):
        year_link = "https://www.federalreserve.gov/monetarypolicy/fomchistorical" + str(year) +".htm"
        ###Yearly page search###
        year_page = requests.get(year_link)
        year_page_soup = BeautifulSoup(year_page.content, 'html.parser')
        year_page_content = year_page_soup.find(id="article")
        meeting_contents = year_page_content.find_all(class_="panel panel-default")
        inyear_count = 0
        for meeting_content in meeting_contents:
            try:
                beige_book_link = "https://www.federalreserve.gov" + \
                 meeting_content.find('a', href=re.compile(r'.+beigebook.+\.htm')).get('href')
                beige_book_date = re.findall(r'[0-9]{6}', beige_book_link)[0]
            except:
                beige_book_link = ""
                beige_book_date = ""
            try:
                statement_link = "https://www.federalreserve.gov" + \
                 meeting_content.find('a', text='Statement').get('href')
            except:
                statement_link = ""
            try:
                transcript_link = "https://www.federalreserve.gov" + \
                 meeting_content.find('a', href=re.compile(r'.+meeting\.pdf')).get('href')
                meeting_date = transcript_link[-19:-11]
            except:
                meeting_date = ""
                transcript_link = ""
            link_df = pd.DataFrame({'Year':[year], 'Meeting_Date':[meeting_date], 'Beige_Book_Date':[beige_book_date],\
                                    'Beige_Book_Link':[beige_book_link], 'Statement_Link':[statement_link], 'Transcript_Link':[transcript_link]})
            df_links = df_links.append(link_df)
            inyear_count += 1
            print('finish %d' % inyear_count)
        print("%d finished" % year)
    return df_links

if __name__ == "__main__":
    ###Get the Links###
    df_links = get_links(1970, 1992)
    df_links.to_csv("FOMC_Links.csv", index=False)
    
    ###Do the scrapping###
    df_beige_books = pd.DataFrame({'Year':[], 'Meeting_Date':[], 'Beige_Book_Date':[],\
                                   'div_summary':[],'div_boston':[],\
                                   'div_new_york':[],'div_philadelphia':[],\
                                   'div_cleveland':[],'div_richmond':[],\
                                   'div_atlanta':[],'div_chicago':[],\
                                   'div_st_louis':[],'div_minneapolis':[],\
                                   'div_kansas_city':[],'div_dallas':[],\
                                   'div_san_francisco':[]})
    df_statements = pd.DataFrame({'Year':[], 'Meeting_Date':[],\
                                  'Beige_Book_Date':[], 'Statement':[]})
    df_links = pd.read_csv("FOMC_Links.csv", encoding='utf-8')
    df_beige_books, df_statements = get_materials(df_links, df_beige_books, df_statements)
    df_beige_books.to_csv("FOMC_Beige_Books.csv", index=False)
    df_statements.to_csv("FOMC_Statements.csv", index=False)
            
