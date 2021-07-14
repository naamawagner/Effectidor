#!/groups/pupko/modules/python-anaconda3.6.5/bin/python

import os

# constants to use when sending e-mails using the server admin's email address.
ADMIN_EMAIL = 'naamawagner@mail.tau.ac.il' #'orenavram@gmail.com' #'shiranos@gmail.com' 'evolseq@mail.tau.ac.il'
SMTP_SERVER = 'mxout.tau.ac.il'

QUEUE_NAME = 'pupkoweb'
#PIPELINE_NAME = 'ModelTeller'

# general paths
SERVERS_RESULTS_DIR = '/bioseq/data/results'
SERVERS_LOGS_DIR = '/bioseq/data/logs'

WEBSERVER_NAME = 'effectidor'
WEBSERVER_URL = 'https://effectidor.tau.ac.il'
#MODELTELLER_LOG = '/bioseq/modelteller/MODELTELLER_runs.log'
#APACHE_ERROR_LOG = '/var/log/httpd/modelteller.error_log'

RELOAD_INTERVAL = 30
RELOAD_TAGS = f'<META HTTP-EQUIV="REFRESH" CONTENT="{RELOAD_INTERVAL}"/>'

HR_STYLE = 'style="height:1px;border:none;color:#333;background-color:#333;"'

EFFECTIROT_RESULTS_DIR = os.path.join(SERVERS_RESULTS_DIR, 'effectidor')
EFFECTIDOR_LOGS_DIR = os.path.join(SERVERS_LOGS_DIR, 'effectidor')
EFFECTIDOR_RESULTS_URL = os.path.join(WEBSERVER_URL, 'results')
EFFECTIDOR_HTML_DIR = '/data/www/html/effectidor'
EFFECTIDOR_EXEC = '/groups/pupko/naamawagner/T3Es_webserver/scripts/debug'

MAIN_SCRIPT = os.path.join(EFFECTIDOR_EXEC, 'main_T3Es.py')

RESULT_MSG = 'Unresolved error'


CONTAINER_WIDTH = 'width: 80%'
CONTAINER_NO_MARGIN = 'margin: 0 auto'
CONTAINER_FONT = 'font-size: 20px'

CONTAINER_STYLE = f'{CONTAINER_WIDTH}; {CONTAINER_NO_MARGIN}; {CONTAINER_FONT}'

PROCESSING_MSG = f'<i>{WEBSERVER_NAME.upper()}</i> is now processing your request. This page will be automatically ' \
    f'updated every few seconds (until the job is done). You can also reload it manually. Once the job has finished, ' \
    f'several links to the output files will appear below. '

PROGRESS_BAR_ANCHOR = '''<!--progress_bar_anchor-->'''
PROGRESS_BAR_TAG = '''<div class="progress">
        <div class="progress-bar progress-bar-striped active" role="progressbar" style="width:100%">
        </div>
    </div>'''