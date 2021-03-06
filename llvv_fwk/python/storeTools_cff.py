import FWCore.ParameterSet.Config as cms
import os,sys
import getopt
import commands
import re
import json

def natural_sort(l):
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)

"""
lists the files available in castor
"""
def fillFromStore(dir,ffile=0,step=-1,generatePfn=True):

    localdataset=cms.untracked.vstring()
    if(len(dir)==0) : return localdataset
 
    #check if it is a directory (check if it is castor, eos or local)
    prefix='singlefile'
    lsout=[dir]
    if(dir.find('path=')>=0) :
        print 'Using dbs to query %s'%(dir)
        prefix='eoscms'
        lsout=commands.getstatusoutput('dbs lsf --' + dir)[1].split()

    elif(dir.find('castor')>=0) :
        prefix='rfio'
        lscommand ='rfdir ' + dir + ' | awk \'{print $9}\''
        lsout = commands.getstatusoutput(lscommand)[1].split()

    elif(dir.find('/store/')==0):
        prefix='eoscms'
        lscommand = 'cmsLs -R ' + dir + ' | grep root | awk \'{print $5}\''
        lsouttmp = commands.getstatusoutput(lscommand)[1].split()

        #this is needed as cmsLs lists twice files staged from castor
        #(only needed during transition to EOS, can remove now)
        lsout=[]
        nduplicate=0
        for l in lsouttmp :
            if len(l)==0 : continue
            if l in lsout :
                nduplicate += 1
                continue
            #filter out CMG trees and histograms
            basename = os.path.basename(l)
            if(basename.find('tree_')==0) : continue
            if(basename.find('histogram')==0): continue
            lsout.append(l)
#        print 'Discarded ' + str(nduplicate)  + ' files duplicated in cmsLs output'
       
    elif(dir.find('.root')<0):
        prefix='file'
        lscommand='ls ' + dir
        lsout = commands.getstatusoutput(lscommand)[1].split()

    #check for the files needed (first file, firstfile+step)
    ifile=0
    for line in lsout :
        if(type(line) is not str ) : continue
        if(len(line)==0) : continue
        if(line.find('root')<0) : continue
        if(ifile<ffile): continue

        if( (step<0) or  (step>0 and ifile<ffile+step) ):                
            sline=''
            if(prefix=='eoscms') :
                if(generatePfn) : sline=commands.getstatusoutput('cmsPfn ' + line )[1]
                else            : sline=line
            elif(prefix=='singlefile') :
                sline='file://' + line
            else :
                sline=str(prefix+'://' + dir + '/' + line.split()[0])
                if(len(sline)==0): continue
            sline = sline.replace('?svcClass=default', '')

            localdataset.extend( [ sline ] )
    ifile=ifile+1

    return natural_sort(localdataset)

"""
check that a file exist and is not corrupted
"""
def checkInputFile(url):
    if(url.startswith('/store')==True):
       url= 'root://eoscms//eos/cms'+url
    command_out = commands.getstatusoutput("root -l -b -q " + url)
    if(command_out[1].find("Error")>=0 or command_out[1].find("probably not closed")>=0 or command_out[1].find("Corrupted")>=0):return False
    return True


"""
check store for duplicates
"""
def checkStoreForDuplicates(outdir):
    ls_cms = "ls " + outdir
    isEOS=False
    isCastor=False
    if(outdir.find('/store/')==0) :
        isEOS=True
        splitOnString=','
        ls_cms='cmsLs ' + outdir + ' | grep root | awk \'{print $5}\''	
    elif(outdir.find('/store/')==0) :
        isCastor = True
        ls_cms = "rfdir " + outdir + + ' | grep root'
    nOutFile = 0
    lsCmd_out = commands.getstatusoutput(ls_cms)
    jobNumbers = []
    duplicatedJobs = []
    origFiles=[]
    duplicatedFiles=[]
    if lsCmd_out[0] == 0:
        lsCmd_outLines = lsCmd_out[1].split("\n")
        if len(lsCmd_outLines) != 0:
            for fileLine in lsCmd_outLines:
                if not "root" in fileLine: continue
                fileName=fileLine
                if(isCastor) : fileName = fileLine.split()[8]

	        if(checkInputFile(fileName)==True):
		    jobNumber=-1
		    try:
			fileBaseName=os.path.basename(fileName)
			jobNumber=int(fileBaseName.split("_")[1])
		    except:
			continue

		    if jobNumber in jobNumbers:
			if not jobNumber in duplicatedJobs:  duplicatedJobs.append(jobNumber)
			duplicatedFiles.append(fileName)
		    else :
			jobNumbers.append(jobNumber)
			origFiles.append(fileName)
			nOutFile += 1
   	        else:
		    print("   #corrupted file found : " + fileName)
		    duplicatedFiles.append(fileName)
    return natural_sort(duplicatedFiles)


"""
clean up for duplicats in the storage area
"""
def removeDuplicates(dir):
    duplicatedFiles=checkStoreForDuplicates(dir)
    print 'Removing ' + str(len(duplicatedFiles)) + ' duplicated files in ' + dir
    isNCG=False
    isEOS=False
    isCastor=False
    if(dir.find('/lustre/ncg.ingrid.pt/')==0) : isNCG=True
    if(dir.find('/store/')==0) : isEOS=True
    if(dir.find('castor')>=0) : isCastor=True
    print "removing duplicatessssssssssssssssssssss"
    for f in duplicatedFiles :
        print f
        if(isNCG) :
            
            commands.getstatusoutput('') # REPLACE
        elif(isEOS) : commands.getstatusoutput('cmsRm ' + f)
        elif(isCastor) : commands.getstatusoutput('rfrm ' +dir + '/' + f)
        else : commands.getstatusoutput('rm ' +dir + '/' + f)

            
"""
wrapper to read the configuration from comand line
args: castor_directory output_file first_file step
"""
def configureSourceFromCommandLine() :
    storeDir=''
    outputFile='Events.root'
    ffile=0
    step=-1
    try:
        if(len(sys.argv)>2 ):
            if(sys.argv[2].find('/')>=0 or sys.argv[2].find('.root')>0) :
                storeDir=sys.argv[2]
                if(len(sys.argv)>3 ):
                    if(sys.argv[3].find('.root')>0):
                        outputFile=sys.argv[3]
                    if(len(sys.argv)>4 ):
                        if(sys.argv[4].isdigit()) : ffile=int(sys.argv[4])
                        if(len(sys.argv)>5 ):
                            if(sys.argv[5].isdigit()) : step=int(sys.argv[5])
    except:
        print '[storeTools_cff] Could not configure from command line, will return default values'

    return outputFile, fillFromStore(storeDir,ffile,step)



def addPrefixSuffixToFileList(Prefix, fileList, Suffix):
   outList = []
   for s in fileList:
      outList.append(Prefix+s+Suffix)
   return natural_sort(outList)







def keepOnlyFilesFromGoodRun(fileList, jsonPath):
   if(jsonPath == ''): 
     print "No golden json file"
     return fileList

   goodLumis = {}
   print "found golden json file"
   jsonFile = open(jsonPath,'r')
   runList=json.load(jsonFile,encoding='utf-8').items()
   runList.sort()
   for run in runList :
      #print run
      goodLumis[int(run[0])] = []
      for lumi in run[1] : 
         #print lumi
         goodLumis[int(run[0])].append(lumi)

   outFileList = []
   for F in fileList:
      try:
         #FAST technique only based on run number and file name parsing
         if '/00000/' in F:
            Fsplit = F.split('/00000/')[0].split('/')
            run = int(Fsplit[-2])*1000+int(Fsplit[-1])
            #print run 
            if(run in goodLumis): outFileList.extend([F])
         else:
            #SLOW technique
            ## To be fixed: new way to store the data is different compared to the one used for 74X
            #check if the file contains at least one good lumi section using DAS_CLIENT --> commented out because VERY SLOW!
            print 'das_client.py --limit=0 --query "lumi file='+f+' |  grep lumi.run_number,lumi.number"'      
            containsGoodLumi = False
            for run in commands.getstatusoutput('das_client.py --limit=0 --query "run file='+f+'"')[1].replace('[','').replace(']','').split(','):
               if(not IsGoodRun(int(run))):continue
               #comment the part bellow because it makes it very very slow
               #for lumi in commands.getstatusoutput('das_client.py --limit=0 --query "lumi file='+f+'"')[1].replace('[','').replace(']','').split(','):
                  #if(not IsGoodLumi(run, lumi)):continue
               outFileList.extend([F])
      except:
         #if something went wrong, always accept the file
         outFileList.extend([F])

   return outFileList


"""
Takes a directory on eos (starting from /store/...) and returns a list of all files with 'prepend' prepended
"""
def getLslist(directory, mask='', prepend='root://eoscms//eos/cms', local='False'):
    from subprocess import Popen, PIPE
    print 'looking into: '+directory+'...'

    if(directory.find('/lustre/ncg.ingrid.pt/')==0) : 
        isNCG=True
        isLocal=True

    eos_cmd = '/afs/cern.ch/project/eos/installation/0.2.41/bin/eos.select'

    data=""
    if local:
        data = Popen(['ls', directory],stdout=PIPE)
    else:
        data = Popen([eos_cmd, 'ls', '/eos/cms/'+directory],stdout=PIPE)
    
    out,err = data.communicate()

    full_list = []

    ## if input file was single root file:
    if directory.endswith('.root'):
        if len(out.split('\n')[0]) > 0:
            return [prepend + directory]

    ## instead of only the file name append the string to open the file in ROOT
    for line in out.split('\n'):
        if len(line.split()) == 0: continue
        full_list.append(prepend + directory + '/' + line)

    ## strip the list of files if required
    if mask != '':
        stripped_list = [x for x in full_list if mask in x]
        return stripped_list

    print full_list
    ## return 
    return full_list
