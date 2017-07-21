%insert("python"){

import os
_pycadir = os.path.join(os.path.expanduser('~'), '.PyCA')
_statfile = os.path.join(_pycadir, 'reportedStatistics.txt')

def _choosePrompt(prompt, defaultYes):
    c = "_"
    while True:
        c = raw_input(prompt + ' ')
        if c in ['', 'y', 'Y', 'n', 'N']:  # valid choices
            break

    if c is '':
        return defaultYes
    elif c in ['y', 'Y']:
        return True
    elif c in ['n', 'N']:
        return False
    else:
        raise Exception("Somehow got a nonsensical choice from prompt")

if not os.path.exists(_statfile):
    print """Welcome to PyCA.

Since this is your first time loading the PyCA.Core module, we'd like to ask if
you would be comfortable providing some feedback.  Your responses to the
following four brief questions help us taylor development and these
statistics are also very useful in helping us fund the project.  We appreciate
your responses, but if you do not feel comfortable sharing this information
with us, you may elect to skip this step.  Also note that any responses may be
left blank and that you will not see this message again.

Thank you,
The PyCA team

"""
    msg = None  # default to no message
    if _choosePrompt("Answer four brief questions? (Y/n)", True):
        name = raw_input("Your Name: ")
        inst = raw_input("Institution: ")

        import urllib2
        ret = urllib2.urlopen('https://enabledns.com/ip')
        ip = ret.read()
        if not _choosePrompt("Report IP Address " + ip + " (Y/n)", True):
            ip = ''

        use = raw_input("Briefly describe what you plan to use PyCA for: ")

        msg = "name: %s\ninst: %s\n ip: %s\nuse: %s" % (name, inst, ip, use)

    try:
        os.makedirs(_pycadir)
    except OSError as exc: # Python >2.5
        pass

    if msg is not None:
        print "Thank you for helping PyCA development"

        # save this information locally in _statfile
        with open(_statfile, 'w') as f:
            f.write(msg)

        try:
            import urllib, urllib2
            url = 'https://docs.google.com/forms/d/1Hiz9cjzr3eeN-DnbhNGiiEYzhArGgPnCm1GzsB0zbtI/formResponse'
            values = {'entry.988845752': name,
                      'entry.1776406582': inst,
                      'entry.1681972327': ip,
                      'entry.1178503519': use}
            data = urllib.urlencode(values)
            req = urllib2.Request(url, data)
            response = urllib2.urlopen(req)
        except Exception, detail:
            print "Error sending", detail


        print "Thank you!"
    else:  # still write a stub so we know not to ask again
        print "Saving empty file %s" % (_statfile)
        with open(_statfile, 'w') as f:
            f.write('')

}
