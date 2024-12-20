#!/groups/pupko/modules/python-anaconda3.6.5/bin/python

import os
import sys
import cgi
import cgitb
import subprocess
from time import time, ctime
from random import randint

if os.path.exists('/bioseq/effectidor'): #remote run
    sys.path.append('/bioseq/effectidor/auxiliaries/')
    sys.path.append('/bioseq/bioSequence_scripts_and_constants/')

import effectidor_CONSTANTS as CONSTS # from /bioseq/effectidor/auxiliaries/
from directory_creator import create_dir # from /bioseq/bioSequence_scripts_and_constants/
from email_sender import send_email # from /bioseq/bioSequence_scripts_and_constants/

def print_hello_world(output_path = '', run_number = 'NO_RUN_NUMBER'):

    hello_world_html = """
    <html>
        <body>
            <h2>Hello World! """ + run_number + """</h2>
        </body>
    </html>
    """
    if not output_path:
        print(hello_world_html)
    else:
        with open(output_path, 'w') as f:
            f.write(hello_world_html)

def write_to_debug_file(cgi_debug_path, content):
    with open(cgi_debug_path, 'a') as f:
        f.write(f'{ctime()}: {content}\n')

def write_html_prefix(output_path, run_number):
    with open(output_path, 'w') as f:
        f.write(f'''<html><head>

    <meta http-equiv="cache-control" content="no-cache, must-revalidate, post-check=0, pre-check=0" />
    <meta http-equiv="cache-control" content="max-age=0" />
    <meta http-equiv="expires" content="0" />
    <meta http-equiv="expires" content="Tue, 01 Jan 1980 1:00:00 GMT" />
    <meta http-equiv="pragma" content="no-cache" />
    {CONSTS.RELOAD_TAGS}

    <title>Effectidor Job #{run_number}</title>
    <link rel="shortcut icon" type="image/x-icon" href="{CONSTS.WEBSERVER_URL}/pics/logo.gif" />

    <meta charset="utf-8">
    <!--<meta name="viewport" content="width=device-width, initial-scale=1">-->
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/css/bootstrap.min.css">
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.3.1/jquery.min.js"></script>
    <script src="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.7/js/bootstrap.min.js"></script>
    <link rel="stylesheet" href="https://gitcdn.github.io/bootstrap-toggle/2.2.2/css/bootstrap-toggle.min.css">

    <link rel="stylesheet" href="{CONSTS.WEBSERVER_URL}/css/general.css">

    </head><body>
    <nav role="navigation" class="navbar navbar-fixed-top">
        <div class="jumbotron" id="jumbo">
            <div class="container">            
                <div class="row" id="title-row" align="center">
                    <div class="col-md-1">
                    </div>
                    <div class="col-md-10" align="center">
                        <!--<span id="server-title">Effectidor</span><br>-->
                        <img src="{CONSTS.WEBSERVER_URL}/pics/logo.gif" id="nav_bar_image" style="height: 120px;"><br>
                        <span id="sub-title">A machine learning based type 3 effectors predictor</span>
                        <br><br>
                    </div>
                </div>
            </div>       
        </div>
    </nav>
    <div id="behind-nav-bar-results">
    </div>
    <br><div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify"> 
    <H1 align=center>Job Status - <FONT color='red'>RUNNING</FONT></h1>
    <br>Effectidor is now processing your request. This page will be automatically updated every {CONSTS.RELOAD_INTERVAL} seconds (until the job is done). You can also reload it manually. Once the job has finished, the output will appear below. A link to this page was sent to your email in case you wish to view these results at a later time without recalculating them. Please note that the results will be kept in the server for 3 months.
    <br><br></div>''')
        f.flush()


def upload_file(form, form_key_name, file_path, cgi_debug_path):
    write_to_debug_file(cgi_debug_path, f'{"#"*80}\nuploading file\n')
    filename = form[form_key_name].filename
    write_to_debug_file(cgi_debug_path, f'file name is:\n{filename}\n')
    content = form[form_key_name].value
    write_to_debug_file(cgi_debug_path, f'{filename} first 100 chars are: {content[:100]}\n')
    with open(file_path, 'wb') as f:
        f.write(content)


def write_running_parameters_to_html(cgi_debug_path, output_path, job_title, ORFs_name, effectors_name, T3Es_name,
                                     host_name, non_T3SS_name, gff_name, genome_f_name, PIP, hrp, mxiE, exs, tts,
                                     T3signal):
    with open(output_path, 'a') as f:
        write_to_debug_file(cgi_debug_path, f'24\n') 
        # regular params row
        f.write(f"""<div class="container" style="{CONSTS.CONTAINER_STYLE}" align="justify">""")
        write_to_debug_file(cgi_debug_path, f'25\n')
        if job_title != '':
            write_to_debug_file(cgi_debug_path, f'26\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>Job title: </b>{job_title}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'27\n')
        write_to_debug_file(cgi_debug_path, f'28\n')
        f.write('<div class="row" style="font-size: 20px;">')
        f.write('<div class="col-md-12">')
        f.write(f'<b>ORFs input file: </b>{ORFs_name}')
        f.write('</div></div>')
        write_to_debug_file(cgi_debug_path, f'29\n')
        if effectors_name != '':
            write_to_debug_file(cgi_debug_path, f'30\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>effectors input file: </b>{effectors_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'31\n')
            
        if T3Es_name !='':
            write_to_debug_file(cgi_debug_path, f'32\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>T3Es of other bacteria input file: </b>{T3Es_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'33\n')
            
        if host_name != '':
            write_to_debug_file(cgi_debug_path, f'34\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>host input file: </b>{host_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'35\n')
            
        if non_T3SS_name != '':
            write_to_debug_file(cgi_debug_path, f'36\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>bacterial proteomes without T3SS archive: </b>{non_T3SS_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'37\n')
        
        if gff_name:
            write_to_debug_file(cgi_debug_path, f'40\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>GFF file:</b>{gff_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'41\n')
            
        if genome_f_name:
            write_to_debug_file(cgi_debug_path, f'42\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>full genome file:</b>{genome_f_name}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'43\n')
            
        if PIP or hrp or mxiE or exs or tts:
            write_to_debug_file(cgi_debug_path, f'44\n')
            to_print = []
            if PIP:
                write_to_debug_file(cgi_debug_path, f'45\n')
                to_print.append('PIP-box')
            if hrp:
                write_to_debug_file(cgi_debug_path, f'46\n')
                to_print.append('hrp-box')
            if mxiE:
                write_to_debug_file(cgi_debug_path, f'47\n')
                to_print.append('mxiE-box')
            if exs:
                write_to_debug_file(cgi_debug_path, f'48\n')
                to_print.append('exs-box')
            if tts:
                write_to_debug_file(cgi_debug_path, f'49\n')
                to_print.append('tts-box')
            write_to_debug_file(cgi_debug_path, f'50\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>Cis regulatory elements:</b>{",".join(to_print)}')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'51\n')
        if T3signal:
            write_to_debug_file(cgi_debug_path, f'52\n')
            f.write('<div class="row" style="font-size: 20px;">')
            f.write('<div class="col-md-12">')
            f.write(f'<b>Including Type III secretion signal prediction</b>')
            f.write('</div></div>')
            write_to_debug_file(cgi_debug_path, f'53\n')
        f.write('</div><br>')
        write_to_debug_file(cgi_debug_path, f'57\n')


def write_cmds_file(cmds_file, run_number, parameters):
    # the queue does not like very long commands so I use a dummy delimiter (!@#) to break the commands for q_submitter
    new_line_delimiter = ';!@#'
    # the code contains features that are exclusive to Python3.6 (or higher)!
    required_modules = ' '.join(['MMseqs2/May2024'])
    with open(cmds_file, 'w') as f:
        f.write(f'module load {required_modules}')
        f.write(new_line_delimiter)
        cmd = " ".join(["python", f"{CONSTS.MAIN_SCRIPT}", parameters])+f'\teffectidor{run_number}\n'
        f.write(cmd)

def run_cgi():

    # prints detailed error report on BROWSER when backend crashes
    # This line MUST appear (as is) BEFORE any error occurs to get a report about the exception!! otherwise you'll get a non-informatvie message like "internal server error"
    cgitb.enable()

    # print_hello_world() # for debugging
    form = cgi.FieldStorage()  # extract POSTed object

    # random_chars = "".join(choice(string.ascii_letters + string.digits) for x in range(20))
    run_number = str(round(time())) + str(randint(10 ** 19, 10 ** 20 - 1))  # adding 20 random digits to prevent users see data that are not their's
    if False:
        run_number = 'debug'  # str(round(time())) + str(randint(1000,9999)) # adding 4 random figures to prevent users see data that are not their's

    results_url = os.path.join(CONSTS.EFFECTIDOR_RESULTS_URL, run_number)
    output_url = os.path.join(results_url, 'output.html')

    wd = os.path.join(CONSTS.EFFECTIROT_RESULTS_DIR, run_number)
    create_dir(wd)
    output_path = os.path.join(wd, 'output.html')
    cgi_debug_path = os.path.join(wd, 'cgi_debug.txt')
    #print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
    # print_hello_world(wd+'/test.txt') # comment out for debugging
    # print_hello_world(output_html_path, run_number) # comment out for debugging

    write_html_prefix(output_path, run_number)  # html's prefix must be written BEFORE redirecting...

    print(f'Location: {output_url}')  # Redirects to the results url. MUST appear before any other print.
    print('Content-Type: text/html\n')  # For more details see https://www.w3.org/International/articles/http-charset/index#scripting
    sys.stdout.flush()  # must be flushed immediately!!!

    # email field should ALWAYS exist in the form (even if it's empty!!)
    # if it's not there, someone sent a request not via the website so they should be blocked.
    # confirm_email is hidden field that only spammer bots might fill in...
    if 'email' not in form or ('confirm_email' in form and form['confirm_email'].value != ''):
        exit()

    # uncomment to send the admin a notification email EVERY time there's a new request
    send_email(smtp_server=CONSTS.SMTP_SERVER, sender=CONSTS.ADMIN_EMAIL,
               receiver=CONSTS.OWNER_EMAIL, subject=f'Effectidor - A new job has been submitted: {run_number}',
               content=f"{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, 'cgi_debug.txt')}\n{os.path.join(CONSTS.WEBSERVER_URL, 'results', run_number, 'output.html')}")

    try:
        if form['email'].value != '':
            write_to_debug_file(cgi_debug_path, f"{form['email'].value.strip()}\n\n")
        if len(form['ORFs'].value) < 100:
            raise ValueError('ORFs file is too small')

        with open(cgi_debug_path, 'a') as f:
            # for debugging
            f.write(f'{"#"*50}\n{ctime()}: A new CGI request has been recieved!\n')
            sorted_form_keys = sorted(form.keys())
            f.write(f'These are the keys that the CGI received:\n{"; ".join(sorted_form_keys)}\n\n')
            f.write('Form values are:\n')
            for key in sorted_form_keys:
                if 'ORFs' in key or 'effectors' in key or 'host' in key or 'no_T3SS' in key or 'gff' in key or 'genome' in key or 'T3Es' in key:
                    # avoid writing the whole file
                    f.write(f'100 first characters of {key} = {form[key].value[:100]}\n')
                else:
                    f.write(f'{key} = {form[key]}\n')
            f.write('\n\n')

        # extract form's values:
        user_email = form['email'].value.strip()

        job_title = ''
        if form['job_title'].value != '':
            job_title = form['job_title'].value.strip()
        write_to_debug_file(cgi_debug_path, f'1\n')
        # This is hidden field that only spammer bots might fill in...
        confirm_email_add = form['confirm_email'].value  # if it is containing a value it is a spammer.
        write_to_debug_file(cgi_debug_path, f'2\n')
        # human readable parameters for results page and confirmation email
        ORFs_name = form['ORFs'].filename
        write_to_debug_file(cgi_debug_path, f'3\n')
        effectors_name = form['effectors'].filename
        write_to_debug_file(cgi_debug_path, f'4\n')
        T3Es_name = form['T3Es'].filename
        write_to_debug_file(cgi_debug_path, f'5\n')
        host_name = form['host'].filename
        write_to_debug_file(cgi_debug_path, f'6\n')
        non_T3SS_name = form['no_T3SS'].filename
        write_to_debug_file(cgi_debug_path, f'7\n')
        #full_genome = form['full_genome'].value.strip()
        write_to_debug_file(cgi_debug_path, f'8\n')
        gff_name = form['gff'].filename
        write_to_debug_file(cgi_debug_path, f'9\n')
        genome_f_name = form['genome'].filename
        write_to_debug_file(cgi_debug_path, f'10\n')
        effectors_homology = form['homology_search'].value.strip()
        translocation_signal = form['translocation_signal'].value.strip()
        identity = form['identity_cutoff'].value.strip()
        coverage = form['coverage_cutoff'].value.strip()
        if ORFs_name.endswith('zip'): #ZIP archive
            ORFs_path = os.path.join(wd, 'ORFs.zip')
        else: #FASTA
            ORFs_path = os.path.join(wd, 'ORFs.fasta')
        upload_file(form, 'ORFs', ORFs_path, cgi_debug_path)
        write_to_debug_file(cgi_debug_path, f'ORFs file was saved to disk successfully\n\n')
        
        parameters = f'{ORFs_path} {wd} --html_path {output_path} -q {CONSTS.QUEUE_NAME}'
        if form['effectors'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'11\n')
            if effectors_name.endswith('zip'): # ZIP archive
                effectors_path = os.path.join(wd, 'effectors.zip')
            else: # FASTA
                effectors_path = os.path.join(wd, 'effectors.fasta')
            upload_file(form, 'effectors', effectors_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path, f'effectors file was saved to disk successfully\n\n')
            
            parameters += f' --input_effectors_path {effectors_path}'
            
            if effectors_homology == 'yes':
                parameters += f' --homology_search'
        if translocation_signal == 'yes':
            write_to_debug_file(cgi_debug_path, 'include translocation signal')
            parameters += ' --translocation_signal'
            
        if form['T3Es'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'12\n')
            T3Es_path = os.path.join(wd, 'user_T3Es.fasta')
            upload_file(form, 'T3Es', T3Es_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path, f'user_T3Es file was saved to disk successfully\n\n')
            
            parameters += f' --input_T3Es_path {T3Es_path}'
            
        if form['host'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'13\n')
            host_path = os.path.join(wd, 'host.zip')
            upload_file(form, 'host', host_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path, f'host file was saved to disk successfully\n\n')
            
            parameters += f' --host {host_path}'
            
        if form['no_T3SS'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'14\n')
            no_T3SS_path = os.path.join(wd, 'non_T3SS.zip')
            upload_file(form, 'no_T3SS', no_T3SS_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path, f'no_T3SS file was saved to disk successfully\n\n')
        
            parameters += f' --no_T3SS {no_T3SS_path}'

        parameters += f' --identity_cutoff {identity}'
        parameters += f' --coverage_cutoff {coverage}'
        PIP, hrp, mxiE, exs, tts = False, False, False, False, False
        if form['gff'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'16\n')
            if gff_name.endswith('.zip'):
                gff_path = os.path.join(wd, 'genome_features.zip')
            else:
                gff_path = os.path.join(wd, 'genome.gff3')
            upload_file(form,'gff', gff_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path,'gff file was saved to disc successfully\n\n')
            
            parameters += f' --gff_path {gff_path}'
        
            write_to_debug_file(cgi_debug_path, f'17\n')
            
        if form['genome'].value: # not empty string / empy bytes - the file was supplied by the user
            write_to_debug_file(cgi_debug_path, f'15\n')
            if genome_f_name.endswith('.zip'):
                genome_path = os.path.join(wd, 'genome_sequence.zip')
            else:
                genome_path = os.path.join(wd, 'genome.fasta')
            upload_file(form, 'genome', genome_path, cgi_debug_path)
            write_to_debug_file(cgi_debug_path, 'genome file was saved to disc successfully\n\n')
            
            parameters += f' --genome_path {genome_path}'
            if 'PIP-box' in form:
                write_to_debug_file(cgi_debug_path, f'18\n')
                parameters += ' --PIP'
                PIP = True
            if 'hrp-box' in form:
                write_to_debug_file(cgi_debug_path, f'19\n')
                parameters += ' --hrp'
                hrp = True
            if 'mxiE-box' in form:
                write_to_debug_file(cgi_debug_path, f'20\n')
                parameters += ' --mxiE'
                mxiE = True
            if 'exs-box' in form:
                write_to_debug_file(cgi_debug_path, f'21\n')
                parameters += ' --exs'
                exs = True
            if 'tts-box' in form:
                write_to_debug_file(cgi_debug_path, f'22\n')
                parameters += ' --tts'
                tts = True
        T3signal = False
        if translocation_signal == 'yes':
            T3signal = True
        write_to_debug_file(cgi_debug_path, f'23\n')   
        write_running_parameters_to_html(cgi_debug_path, output_path, job_title, ORFs_name, effectors_name, T3Es_name,
                                         host_name, non_T3SS_name, gff_name, genome_f_name, PIP, hrp, mxiE, exs, tts,
                                         T3signal)
        write_to_debug_file(cgi_debug_path, f'{ctime()}: Running parameters were written to html successfully.\n')

        cmds_file = os.path.join(wd, 'qsub.cmds')
        write_cmds_file(cmds_file, run_number, parameters)

        log_file = cmds_file.replace('cmds', 'log')
        # complex command with more than one operation (module load + python q_submitter.py)
        # submission_cmd = 'ssh bioseq@powerlogin "module load python/anaconda_python-3.6.4; /bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_file} {wd} -q {CONSTS.QUEUE_NAME} --verbose > {log_file}"'

        # simple command when using shebang header
        submission_cmd = f'ssh bioseq@power9login /bioseq/bioSequence_scripts_and_constants/q_submitter_power.py {cmds_file} {wd} -q {CONSTS.QUEUE_NAME} --verbose > {log_file}'

        write_to_debug_file(cgi_debug_path, f'\nSSHing and SUBMITting the JOB to the QUEUE:\n{submission_cmd}\n')

        subprocess.call(submission_cmd, shell=True)

        if user_email != '':
            with open(os.path.join(wd, 'user_email.txt'), 'w') as f_email:
                f_email.write(f'{user_email}\n')

            notification_content = f"Your submission configuration is:\n\n"
            if job_title:
                notification_content += f'Job title: {job_title}\n'
            notification_content += f'ORFs file: {ORFs_name}\n' 
            if effectors_name:
                notification_content += f'known effectors file: {effectors_name}\n'
            if T3Es_name:
                notification_content += f'type 3 effectors of other bacteria file: {T3Es_name}\n'
            if host_name:
                notification_content += f'host proteome file: {host_name}\n'
            if non_T3SS_name:
                notification_content += f'proteomes with no T3SS archive: {non_T3SS_name}\n'
            if gff_name:
                notification_content += f'GFF file: {gff_name}\n'
            if genome_f_name:
                notification_content += f'genome file: {genome_f_name}\n'
            notification_content += f'Minimal identity percentage of orthologs and paralogs: {identity}\n'
            notification_content += f'Minimal coverage percentage of orthologs and paralogs: {coverage}\n'
            notification_content += f'\nYou can track the progress of your job at:\n{os.path.join(CONSTS.WEBSERVER_URL, "results", run_number, "output.html")}\n\n'

            # Send the user a notification email for their submission
            try:
                send_email(smtp_server=CONSTS.SMTP_SERVER,
                           sender=CONSTS.ADMIN_EMAIL,
                           receiver=f'{user_email}',
                           subject=f'Effectidor - Your job has been submitted!{" (Job name: "+str(job_title) if job_title else ""})',
                           content=notification_content)
            except:
                write_to_debug_file(cgi_debug_path,'\nInvalid Email - Email was not sent to the user\n\n')

        write_to_debug_file(cgi_debug_path, f'\n\n{"#"*50}\nCGI finished running!\n{"#"*50}\n')

    except Exception as e:
        msg = 'CGI crashed before the job was submitted :('
        with open(output_path) as f:
            html_content = f.read()
        html_content = html_content.replace('RUNNING', 'FAILED')
        html_content += f'<br><br><br><center><h2><font color="red">{msg}</font><br><br>Please make sure all files are in the right format and try to re-run your job or <a href="mailto:{CONSTS.ADMIN_EMAIL}?subject={CONSTS.PIPELINE_NAME}%20Run%20Number%20{run_number}">contact us</a> for further information</h2></center><br><br>\n</body>\n</html>\n'
        with open(output_path, 'w') as f:
            f.write(html_content)

        exc_type, exc_obj, exc_tb = sys.exc_info()
        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
        write_to_debug_file(cgi_debug_path, f'\n{"$"*50}\n\n{msg}\n\n{fname}: {exc_type}, at line: {exc_tb.tb_lineno}\n\n{"$"*60}')

        # logger.info(f'Waiting {2*CONSTS.RELOAD_INTERVAL} seconds to remove html refreshing headers...')
        # Must be after flushing all previous data. Otherwise it might refresh during the writing.. :(
        from time import sleep

        sleep(2 * CONSTS.RELOAD_INTERVAL)
        with open(output_path) as f:
            html_content = f.read()
        html_content = html_content.replace(CONSTS.RELOAD_TAGS, '')
        with open(output_path, 'w') as f:
            f.write(html_content)

    # logging submission
    with open(CONSTS.SUBMISSIONS_LOG, 'a') as f:
        f.write(f'{user_email}\t{run_number}\t{ctime()}\n')

if __name__ == '__main__':
    run_cgi()