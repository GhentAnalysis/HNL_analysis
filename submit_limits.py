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

for sgn in signals:
    ## Create list of samples
    label = sgn.split()[0]
    if os.path.isfile('jobs/completed/done_'+label): continue
    outlist = open('jobs/samples/samples_'+label+'.txt', 'w+')
    outlist.write(data_backgrounds[0])
    outlist.write(sgn)
    for i in range(1, len(data_backgrounds)): outlist.write(data_backgrounds[i])
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
    os.system('qsub jobs/jobfiles/limit_job_'+label+'.sh -q express -o jobs/outfiles/output_'+label+'.stdout -e jobs/outfiles/output_'+label+'.stdout > jobs/outfiles/submitCheck_'+label+'.out 2>> jobs/outfiles/submitCheck_'+label+'.out')  # for sh/bash
    #os.system('qsub jobs/jobfiles/limit_job_'+label+'.csh -q express -o jobs/outfiles/output_'+label+'.stdout -e jobs/outfiles/output_'+label+'.stderr > jobs/outfiles/submitCheck_'+label+'.out 2>> jobs/outfiles/submitCheck_'+label+'.out')  # for (t)csh
