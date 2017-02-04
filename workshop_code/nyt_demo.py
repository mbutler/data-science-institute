import requests
import csv
import time

# intitate while loop
keep_going = True
page_number = 0

while keep_going:
    # pause
    time.sleep(5)
    
    # increase page number
    page_number += 1
    print(page_number)

    # get url response
    url = 'http://api.nytimes.com/svc/search/v2/articlesearch.json?fq=obamacare&page=' + str(page_number) + '&begin_date=20160101&end_date=20161231&api-key=2b621b40460346d889170e2acb66f636'
    response = requests.get(url)
    data = response.json()

    # check whether there are documents
    if len(data['response']['docs']) == 0:
        keep_going = False
    else:
        # if you find documents then cycle through them
        start_doc = 0
        end_doc = len(data['response']['docs'])
        for doc_number in range(start_doc, end_doc):
            temp_headline = data['response']['docs'][doc_number]['headline']
            temp_byline = data['response']['docs'][doc_number]['byline']

            # check whether there is a headline
            if len(temp_headline) == 0:
                temp_headline = 'NA'
            else:
                temp_headline = data['response']['docs'][doc_number]['headline']['main']

            # check whether byline is None
            if temp_byline == None:
                temp_byline == 'NA'
            else:
                # check whether there is a byline
                if len(temp_byline) == 0:
                    temp_byline = 'NA'
                else:
                    temp_byline = data['response']['docs'][doc_number]['byline']['original']

            # create a row
            row = [data['response']['docs'][doc_number]['pub_date'], temp_headline, temp_byline, data['response']['docs'][doc_number]['lead_paragraph']]

            # append the csv
            with open('/Users/brycedietrich/Dropbox/workshops/text_as_data/prelim/nyt_results.csv', 'a') as my_csv:
                # write the row
                data_writer = csv.writer(my_csv)
                data_writer.writerow(row)
