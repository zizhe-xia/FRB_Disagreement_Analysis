'''
Convert the PDF files to text
'''

import os
import sys, getopt
from collections import Counter, defaultdict
import re
from nltk.tokenize import word_tokenize
import pandas as pd
import numpy as np
import pickle

def remove_bad_characters(text):
    return text.replace('‘',"'").replace('’',"'").replace('“','"').replace('”','"')\
               .replace('！','!').replace('—',' ').replace('-',' ')
#    return text

def refine_text_front(text_raw, pattern_head=None, pattern_tail=None):
    try:
        position_front, position_head = re.search(pattern_head, text_raw).span()
    except:
        position_front, position_head = 0,0
    try:
        position_tail = re.search(pattern_tail, text_raw).span()[0]
    except:
        position_tail = len(text_raw)
    dirty_text = re.sub(r'([0-9]{1,2}/[0-9]{1,2}/[0-9]{2} [0-9]+|\[(Laughter|No response)\])', ' ', text_raw[position_head:position_tail])
    text = remove_bad_characters(dirty_text)
    front_matter = text_raw[0:position_front]
    return text, front_matter

def label_transcripts(text_dir, transcript_names, df_labelled_transcripts):
    os.chdir(text_dir)
    for file_name in transcript_names:
        file_name = file_name +".txt"
        meeting_date = int(file_name[-23:-15])
        with open (file_name, "r", encoding="utf-8") as f:
            text_raw = f.read()
        ###Get rid of front and end matter###
        pattern_head = re.compile(r'Transcript of Federal Open Market Committee Meeting of', re.IGNORECASE)
        pattern_tail = re.compile(r'END OF MEETING')
        text, front_matter = refine_text_front(text_raw, pattern_head, pattern_tail)
        ###Parsing the text###
        speaker_list, content_list = [],[]
        speaker_list.append('FRONT_MATTER')
        content_list.append(front_matter)
        pattern_speaker = re.compile(r'(MR\.|MS\.|MRS\.|CHAIRMAN|VICE CHAIRMAN) [A-Z]+\. ')
        match_iterator = re.finditer(pattern_speaker, text)
        first_match = next(match_iterator)
        speech_end = first_match.span()[0]
        speech_begin = first_match.span()[1]
        speaker_list.append(text[speech_end:speech_begin])
        for match in match_iterator:
            speech_end = match.span()[0]
            content_list.append(re.sub(r'(\n|/)', ' ', text[speech_begin:speech_end]))
            speech_begin = match.span()[1]
            speaker_list.append(text[speech_end:speech_begin])
        else:
            content_list.append(re.sub(r'(\n|/)', ' ', text[speech_begin:speech_end]))
        ###Store the parsed text###
        labelled_transcript_df = pd.DataFrame({'Meeting_Date':[int(meeting_date)]*len(speaker_list),\
                                               'Speaker':speaker_list, 'Content':content_list})
        df_labelled_transcripts = df_labelled_transcripts.append(labelled_transcript_df)
    return df_labelled_transcripts

if __name__ == "__main__":
    data_dir = "C:/Users/xiazizhe.HKUPC2.006/Downloads"
    text_dir = "C:/Users/xiazizhe.HKUPC2.006/Downloads"
    os.chdir(data_dir)
    transcript_names = pd.read_csv("FOMC_Links.csv", encoding='utf-8').Transcript_Name.tolist()
    df_labelled_transcripts = pd.DataFrame({'Meeting_Date':[], 'Speaker':[], 'Content':[]})
    df_labelled_transcripts = label_transcripts(text_dir, transcript_names, df_labelled_transcripts)
    os.chdir(data_dir)
    df_labelled_transcripts.to_csv("FOMC_Transcripts.txt", sep="|", index=False)
    df_labelled_transcripts.to_csv("FOMC_Transcripts.csv", index=False)
    df = pd.read_csv("FOMC_Transcripts.txt", sep="|", encoding='utf-8')
    pickle.dump(df, open('FOMC_transcripts_df.pickle', 'wb'))
    