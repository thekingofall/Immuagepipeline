from email import encoders, message, message_from_binary_file
from email.header import Header
from email.mime.text import MIMEText
from email.utils import parseaddr, formataddr
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
import os,sys

import smtplib
from time import time

from pkg_resources import to_filename

keeptime=sys.argv[1]
attachmessage=sys.argv[2]
projectname=sys.argv[3]

def _format_addr(s):       #格式化一个邮件地址
    name, addr = parseaddr(s)
    return formataddr((Header(name, 'utf-8').encode(), addr))

from_addr = '3308450235@qq.com'    #发送者邮箱地址
password = 'lpyllcjtuehpcjah'                   #发送者邮箱密码,不告诉你密码=。=
to_addr = '1732367346@qq.com'           #接收者邮箱地址
smtp_server = 'smtp.qq.com'          #发送者所在的邮箱供应商的MTA地址


#from_addr = input('From: ')
#password = input('Password: ')
# to_addr = input('To: ')
#smtp_server = input('SMTP server: ')
msg= MIMEMultipart()
msg['From'] = _format_addr('Immuclock <%s>' % from_addr)     #邮件头部，发送者信息
msg['To'] = _format_addr('Immuclock.rep <%s>' % to_addr)           #接收者信息
msg['Subject'] = Header(str(projectname)+"   end", 'utf-8').encode()       #邮件主题

msg.attach(MIMEText('The project  '+projectname+ ' about your immuclock keep  '+keeptime+" h", 'plain', 'utf-8'))

if os.path.exists(attachmessage):
    with open(attachmessage,'rb')  as fhandle:
        mime = MIMEBase('html','html',filename='all.html')
        mime.add_header('Content-Disposition', 'attachment', filename=('gbk','','Immuclock.html'))
        mime.add_header('Content-ID', '<0>')
        mime.add_header('X-Attachment-Id', '0')
        # 把附件的内容读进来:
        mime.set_payload(fhandle.read())
        # 用Base64编码:
        encoders.encode_base64(mime)
        # 添加到MIMEMultipart:
        msg.attach(mime)


try:
    server = smtplib.SMTP(smtp_server, 25)
    server.set_debuglevel(1)
    server.login(from_addr, password)
    server.sendmail(from_addr, [to_addr], msg.as_string())
    server.quit()
    print('发送成功')
except smtplib.SMTPException as e:
    print('发送失败')
    print(e)
