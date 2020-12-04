import sys, os
from array import array

dirs = ['jobs/samples', 'jobs/analisi', 'jobs/makefiles', 'jobs/jobfiles', 'jobs/outfiles', 'jobs/completed', 'jobs/executables']
for adir in dirs:
    if not os.path.isdir(adir): os.makedirs(adir)

data_backgrounds = []
signals = []

infile = open(sys.argv[1], 'r')
for line in infile.readlines():
    if '#' in line: continue
    elif 'HeavyNeutrino' in line: signals.append(line)
    else: data_backgrounds.append(line)


def makeCondorFile( script,runscriptname, label):
    submit_file_name = 'condor_'+label+'.sub'
    submit_file = open(submit_file_name, 'w')
    submit_file.write('Universe     = vanilla \n')
    submit_file.write('Executable   = '+script+' \n')
    submit_file.write('should_transfer_files = YES  \n')
    submit_file.write('Log          = '+'jobs/outfiles/output_'+label+'.log \n')
    submit_file.write('Output       = '+'jobs/outfiles/output_'+label+'.out \n')
    submit_file.write('Error        = '+'jobs/outfiles/output_'+label+'.err \n')
    submit_file.write('Queue \n')
    submit_file.close()
    return submit_file_name


for sgn in signals:
    ## Create list of samples
    label = sgn.split()[0]
    if os.path.isfile('jobs/completed/done_'+label): continue
    outlist = open('jobs/samples/samples_'+label+'.txt', 'w+')
    outlist.write(data_backgrounds[0])
    outlist.write(sgn)
    for ibkg in data_backgrounds[1:]: outlist.write(ibkg)
    #for i in range(1, len(data_backgrounds)): outlist.write(data_backgrounds[i])
    ## Create analisi.C
    with open('analisi_TEMPLATE.C', 'r') as inanalisi: analisidata = inanalisi.read()
    analisidata = analisidata.replace('TEMPLABEL', label)
    with open('jobs/analisi/analisi_'+label+'.C', 'w+') as outanalisi: outanalisi.write(analisidata)
    ## Create makefile
    with open('makeFiles/makeAnalisi_TEMPLATE', 'r') as inmake: makedata = inmake.read()
    makedata = makedata.replace('TEMPLABEL', label)
    with open('jobs/makefiles/make_'+label, 'w+') as outmake: outmake.write(makedata)
    ## Create job file
    with open('limit_job_TEMPLATE.sh', 'r') as injob: jobdata = injob.read()  # for sh/bash
    #with open('limit_job_TEMPLATE.csh', 'r') as injob: jobdata = injob.read()  # for (t)csh
    jobdata = jobdata.replace('TEMPLABEL', label).replace('TEMPDIR', os.getcwd())
    with open('jobs/jobfiles/limit_job_'+label+'.sh', 'w+') as outjob: outjob.write(jobdata)  # for sh/bash
    #with open('jobs/jobfiles/limit_job_'+label+'.csh', 'w+') as outjob: outjob.write(jobdata)  # for (t)csh
    os.chmod('jobs/jobfiles/limit_job_'+label+'.sh', 0o777)  # for sh/bash
    #os.chmod('jobs/jobfiles/limit_job_'+label+'.csh', 0o777)  # for (t)csh
    ## Submit job
    runscript_name = 'jobs/executables/analisi_'+label
    script = 'jobs/jobfiles/limit_job_'+label+'.sh'
    submit_file_name = makeCondorFile(script, runscript_name,label)
    os.system('condor_submit '+submit_file_name)
