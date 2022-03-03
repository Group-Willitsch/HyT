import numpy as np
from string import Template
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
import smtplib


MY_ADDRESS = 'labwatch-pc@unibas.ch'
PASSWORD   = 'UltraCold#'

def get_contacts():
    emails = ['thomas.kierspel@unibas.ch']
    return emails

def send_email(text):
    emails = get_contacts() # read contacts

    # set up the SMTP server
    s = smtplib.SMTP(host='smtp-ext.unibas.ch', port=587)
    s.ehlo()
    s.starttls()
    s.login(MY_ADDRESS, PASSWORD)


    # create a message
    msg = MIMEMultipart()
    email = emails[0]

    # setup the parameters of the message
    msg['From']   = MY_ADDRESS
    msg['To']     = email
    msg['Subject']= "Pressure note"
    if len(text)<=1:
        msg.attach(MIMEText('This is serious', 'plain'))
    else:
        msg.attach(MIMEText(str(text), 'plain'))

    #send the message via the server set up earlier.
    s.send_message(msg)
    del msg
    s.quit()



if __name__ == '__main__':
    send_email('WOW is ths really working?')
