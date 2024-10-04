import os
import subprocess
from threading import Lock
import pandas as pd

###########################################################################
#                                                                         #
#                INPUT PARAMETERS                                         #
#                                                                         #
###########################################################################

#Liner
OD                 = [ 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0, 308.0]
WT                 = [  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0,  12.0]
elem_length        = [   2.5,   2.5,   2.5,   2.5,   2.5,   2.5,   2.5,   2.5,  10.0,  10.0,   2.5,   2.5]
no_el_thru_thkness = [    20,    20,    20,    20,    20,    20,    20,    20,    21,    21,     5,     5]
refLength          = [ 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 250.0, 500.0, 250.0]

#Die
dia_out = [ 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0, 266.0]
angle   = [  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0,  15.0]

#Backing pipe
pipeFlag = [  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True,  True] # True or False
pipeOD   = [ 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9, 323.9]
pipeWT   = [  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5,  17.5]

#Material
anisotropy        = [ 1.0,   1.0,   0.3,   1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0]
stressStrainCurve = [ 'E', 'EPP', 'EPP', 'EPP', 'EP', 'EP', 'EP', 'EP', 'EP', 'EP', 'EP', 'EP'] # E for Elastic or EPP for Elastic-Perfectly-Plastic or EP for Elastic-Plastic
yieldCriterion    = [ 'M',   'M',   'M',   'M',  'M', 'DP', 'DP', 'DP', 'DP',  'M',  'M',  'M'] # M for Mises or DP for Drucker-Prager
kFactor           = [ 1.0,   1.0,   1.0,   1.0,  1.0,  1.0,  0.8,  1.0,  1.0,  1.0,  1.0,  1.0]
angleFriction     = [13.7,  13.7,  13.7,  13.7, 13.7, 13.7, 30.0, 13.7, 13.7, 13.7, 13.7, 13.7]
referenceStress   = [16.1,  16.1,  16.1,  48.3, 16.1, 16.1, 16.1, 16.1, 16.1, 16.1, 16.1, 16.1]
relaxation        = [   0,    0,      0,     0,    0,    0,    0,    0,    0,    1,    2,    2] # 0 for no relaxation, 1 for viscoelastic and 2 for two-layer viscoplastoc

#PE and Steel Friction
friction = [ 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

#Load details
temperature = [ 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0]
run_job     = [ True, True, True, True, True, True, True, True, True, True, True, True] # True or False

#Boundary conditions
BC = ['base','base','base','base','base','base','base','alt','base','alt','alt','base'] # base or alt

#Element type
elementType = ['axis','axis','axis','axis','axis','axis','axis','axis','shell','shell','axis','axis'] # axis or shell

###########################################################################
#                                                                         #
#                SCRIPT                                                   #
#                                                                         #
###########################################################################

for i in range(11,12):#len(OD)):y

    filename = '%iOD_%idegC_%ssSC_%syC_%srefS_%saniS_%skFactor_%saF_%sbP_%sBC_%s' % (OD[i], temperature[i], stressStrainCurve[i], yieldCriterion[i], str(referenceStress[i]).replace('.','d'), str(anisotropy[i]).replace('.','d'), str(kFactor[i]).replace('.','d'), str(angleFriction[i]).replace('.','d'), pipeFlag[i], BC[i], elementType[i])
    
    lock = Lock()
    lock.acquire()
    proc = subprocess.Popen("abq2019 cae nogui=create_model_library.py", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    #proc = subprocess.Popen("abq2019 cae script=create_model_library.py", shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    info = f'{OD[i]}\\{WT[i]}\\{elem_length[i]}\\{no_el_thru_thkness[i]}\\{refLength[i]}\\{dia_out[i]}\\{angle[i]}\\{pipeFlag[i]}\\{pipeOD[i]}\\{pipeWT[i]}\\{anisotropy[i]}\\{stressStrainCurve[i]}\\{yieldCriterion[i]}\\{kFactor[i]}\\{angleFriction[i]}\\{referenceStress[i]}\\{relaxation[i]}\\{friction[i]}\\{temperature[i]}\\{BC[i]}\\{filename}\\{elementType[i]}'
    utf_info = info.encode()
    proc.stdin.write(utf_info)
    proc.stdin.close()
    proc.stdout.read(1)
    lock.release()

    if elementType[i] == 'axis':
        #MPC
        df = pd.DataFrame()
        flag = False
        j = 0
        with open('%s.inp' % filename) as f:
            for line in f:
                if '*Element' in line:
                    break
                if flag == True:
                    df.loc[j,'nodeNo'] = float(line.split('\n')[0].split(',')[0])
                    df.loc[j,'x'] = float(line.split('\n')[0].split(',')[1])
                    df.loc[j,'y'] = float(line.split('\n')[0].split(',')[2])
                    j = j + 1
                if '*Node' in line:
                    flag = True
        df = df.sort_values(by = ['x','y'])
        dfInner = df[df['y'] == 0.5* (OD[i] - 2 * WT[i])]
        dfOuter = df[df['y'] == 0.5 * OD[i]]
        dfMiddle = df[~(df['y'] == 0.5* (OD[i] - 2 * WT[i])) & ~(df['y'] == 0.5 * OD[i])]
        writeLine='*Part, name=Liner\n*MPC\n'
        for k, row in dfMiddle.iterrows():
            writeLine = writeLine + ' SLIDER, %i, %i, %i\n' % (row['nodeNo'],dfInner[dfInner['x'] == row['x']].iloc[0,0],dfOuter[dfOuter['x'] == row['x']].iloc[0,0])

        with open('%s.inp' % filename) as f:
            writeLineNew=f.read().replace('*Part, name=Liner', writeLine[:-1])

        with open('%s.inp' % filename, 'w') as f:
            f.write(writeLineNew)

    #Run job if selected
    if run_job[i] == True:
        os.system(f"abq2019 int ask_delete=off job={filename} cpus=2")
        ext_to_delete = ['.com', '.dat', '.msg', '.prt', '.sim', '.log']
        #for ext in ext_to_delete:
        #    os.system('del {}{}'.format(filename, ext))