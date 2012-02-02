from optparse import OptionParser
import os
import sys

parser = OptionParser()
parser.add_option("-c", "--command", dest="command", help="Input Commands")
parser.add_option("-i", "--inputdirectory", dest="inputdir", help="Location of Data")
parser.add_option("-s", "--outputdirectory", dest="outputdir", help="Final Location of Data")
parser.add_option("-l", "--logdirectory", dest="logdir", default="BatchSubmission_DATE", help="Location of batch logs")
parser.add_option("-q", "--queue", dest="queue", default="1nh", help="Batch queue to submit jobs to")
parser.add_option("--debug", action="store_true", dest="debug", default=False, help="Debug Mode")
parser.add_option("--dryrun", action="store_true", dest="dryrun", default=False, help="Do not submit to batch queue.")
(options, args) = parser.parse_args()

HomeDir = os.getcwd()+"/"
if not os.path.exists(options.logdir): os.mkdir(options.logdir)
InputFiles = os.popen("nsfind "+options.inputdir+" | grep .root").readlines()
print "Submitting %i jobs for directory %s" %(len(InputFiles),options.inputdir)
for infile in InputFiles:
	logfile = options.logdir+"/"+infile.strip("\n")[infile.rfind("/"):]+".log"
	bsubcommand = "bsub -q "+options.queue+" -o "+logfile+" batch_submit.sh "+HomeDir+" \""+options.command+"\" "+infile.strip("\n")+" "+options.outputdir
	bsubcommand = bsubcommand.replace("//","/")
	if options.debug: print bsubcommand
	if not options.dryrun: os.popen(bsubcommand)
print "%i jobs submitted" %len(InputFiles)
print "Done!"
