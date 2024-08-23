import re
import html
from urllib import parse
import requests
import sys

GOOGLE_TRANSLATE_URL = 'http://translate.google.cn/m?q=%s&tl=%s&sl=%s'

def translate(text, to_language="auto", text_language="auto"):

    text = parse.quote(text)
    url = GOOGLE_TRANSLATE_URL % (text,to_language,text_language)
    response = requests.get(url)
    data = response.text
    expr = r'(?s)class="(?:t0|result-container)">(.*?)<'
    result = re.findall(expr, data)
    if (len(result) == 0):
        return ""

    return html.unescape(result[0])

data=sys.argv[1]
datamode=sys.argv[2]
if datamode=="s":
    data_tran=translate(data, "zh-CN","en")
    print(data_tran)    
    # print(translate(data, "zh-CN","en")) #英语转汉语
elif datamode=="m":
    dat=open(file=data,mode="r+")
    datw=open(file=data+".tran",mode="w+")
    for i in dat:
        # print(translate(i, "zh-CN","en"))
        # print(translate(i, "zh-CN","en")) #英语转汉语
        seq=translate(i, "zh-CN","en").strip('"')
        print(seq)
        datw.write(seq+"\n")
        
else:
    print("Error: Please input the right mode.")


