from __future__ import division
import HTMLParser
import urllib
import re
reload(HTMLParser)

class TitleParser(HTMLParser.HTMLParser):
    def __init__(self):
        self.titles = list()
        self.in_title = False
        HTMLParser.HTMLParser.__init__(self)

    def get_titles(self, url):
        page = urllib.urlopen(url)
        raw = page.read()
        # Preprocessing to get rid of url's that aren't in quotes because
        # some of the characters in them screw up the parsera
        processed = re.sub('href=http://.*?>', 'href="">', raw)

        
        # Ridiculous hack to fix a bug I don't understand
        first_value_tag = processed.index('<value')
        processed=processed[first_value_tag:]

        self.feed(processed)
        return self.titles

    def handle_starttag(self, tag, attrs):
        if tag == 'value':
            self.in_title = True

    def handle_endtag(self, tag):
        if tag == 'value':
            self.in_title = False
    
    def handle_data(self, data):
        if self.in_title:
            self.titles.append(data)

