#!/Library/Frameworks/Python.framework/Versions/3.5/bin/python3

import csv
import re
import string
from os import listdir
from bs4 import BeautifulSoup
from nltk import word_tokenize
from nltk.corpus import stopwords
from nltk.stem.snowball import SnowballStemmer

# working directory
wd = "/Users/brycedietrich/Dropbox/workshops/text_as_data/prelim/"

# list files
articles = []
with open(wd + 'nyt_results.csv', 'r') as nyt_csv:
    nyt_data = csv.reader(nyt_csv, delimiter=',')
    for nyt_row in nyt_data:
        articles.append(nyt_row[3])

# preprocessing
stopwords = stopwords.words('english')
stemmer = SnowballStemmer('english')

# positive words
with open(wd + 'pos_words.txt') as f:
    positive_words = f.read().splitlines()

# negative words
with open(wd + 'neg_words.txt') as f:
    negative_words = f.read().splitlines()

# stem dictionaries
positive_words = [stemmer.stem(word) for word in positive_words]
negative_words = [stemmer.stem(word) for word in negative_words]

article_number = 0
for article in articles:
    article_number += 1
    print(article_number)

    # remove punctuation and capitalization
    text = re.sub('\W', ' ', article.lower())

    # get unigrams
    words = word_tokenize(text)

    # remove stop words
    clean_text = filter(lambda x: x not in stopwords, words)
    clean_text = list(clean_text)

    # apply stemmer
    snowball_words = [stemmer.stem(word) for word in clean_text]

    # word count
    word_count = len(snowball_words)

    # count positive
    pos_count = 0
    for word in positive_words:
        if word in snowball_words:
            pos_count += 1

    # count negative
    neg_count = 0
    for word in negative_words:
        if word in snowball_words:
            neg_count += 1

    # append the csv
    with open(wd + 'dictionary_results.csv', 'a') as my_csv:
        # create row
        row = [article_number, pos_count, neg_count, word_count]

        # write the row
        data_writer = csv.writer(my_csv)
        data_writer.writerow(row)
